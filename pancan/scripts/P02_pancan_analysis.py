"""Pan-cancer analysis following plan/pancan_analysis_plan_2026-11-25.md (SHA-256 c440f092...fa3b) and its deviations log.
Run: python3 -I P02_pancan_analysis.py <processed_pancan dir> <ref dir> <breast processed dir> <results dir>
ref dir holds: genes_tss_gencode_v36.tsv, HM450.hg38.mask.tsv.gz, probe_annotation_shared.tsv (breast), DataS1_samples.xlsx,
ABSOLUTE_pancan_purity.txt, cohort_all3_primary.tsv
"""
import glob, json, os, sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore', message='Mean of empty slice')
warnings.filterwarnings('ignore', message='Degrees of freedom')

P, REF, BP, R = sys.argv[1:5]
os.makedirs(R, exist_ok=True)
rng = np.random.default_rng(20261125)
ELIGIBLE = ['KIRP', 'PRAD', 'COAD', 'LUAD', 'STAD', 'ESCA', 'LIHC', 'LUSC', 'KIRC']

# ---------------- shared helpers (identical maths to the breast code) ----------------
def cv_r2(Xs, y):
    """Leave-one-tumour-out R2 for a batch of genes. Xs: (G,n,k) incl. intercept; y: (G,n)."""
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return 1 - (e ** 2).sum(1) / ((y - y.mean(1, keepdims=True)) ** 2).sum(1)

def press(Xs, y):
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return (e ** 2).sum(1)

def zs(X): return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)
def design(cols, G, n): return np.stack([np.ones((G, n))] + [np.broadcast_to(c, (G, n)) for c in cols], axis=2)

def models(M, A, E, cov=()):
    G, n = E.shape; cov = list(cov)
    return (cv_r2(design([zs(M)] + cov, G, n), E), cv_r2(design([zs(A)] + cov, G, n), E),
            cv_r2(design([zs(M), zs(A)] + cov, G, n), E))

def summ(m, a, ma, sel=None):
    if sel is not None: m, a, ma = m[sel], a[sel], ma[sel]
    return {'D_median_A_minus_M': np.median(a - m), 'G_median_MA_minus_A': np.median(ma - a),
            'median_cvR2_M': np.median(m), 'median_cvR2_A': np.median(a), 'median_cvR2_MA': np.median(ma),
            'median_MA_minus_M': np.median(ma - m), 'frac_A_gt_M': np.mean(a > m)}

def jackknife(M, A, E, cov=(), sel=None):
    n = E.shape[1]; full = summ(*models(M, A, E, cov), sel); reps = []
    for i in range(n):
        k = np.arange(n) != i
        reps.append(summ(*models(M[:, k], A[:, k], E[:, k], [c[k] for c in cov]), sel))
    out = {}
    for key in full:
        v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
        out[key] = {'estimate': float(full[key]), 'SE': float(se), 'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
    return out

def spearman_rows(X, Y):
    rx = stats.rankdata(X, axis=1); ry = stats.rankdata(Y, axis=1)
    rx -= rx.mean(1, keepdims=True); ry -= ry.mean(1, keepdims=True)
    return (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))

def dersimonian_laird(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float); w = 1 / se ** 2
    fe = np.sum(w * est) / np.sum(w); Q = np.sum(w * (est - fe) ** 2); k = len(est)
    tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    wr = 1 / (se ** 2 + tau2); re = np.sum(wr * est) / np.sum(wr); sre = np.sqrt(1 / np.sum(wr))
    I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return {'estimate': float(re), 'SE': float(sre), 'CI95': [float(re - 1.96 * sre), float(re + 1.96 * sre)],
            'tau2': float(tau2), 'Q': float(Q), 'I2': float(I2), 'k': int(k)}

# ---------------- reference data ----------------
genes = pd.read_csv(f'{REF}/genes_tss_gencode_v36.tsv', sep='\t')
genes['w0'] = np.where(genes.strand == '+', genes.tss - 1000, genes.tss - 500)
genes['w1'] = np.where(genes.strand == '+', genes.tss + 500, genes.tss + 1000)
genes = genes.drop_duplicates('gene_name').set_index('gene_name')
mk = pd.read_csv(f'{REF}/HM450.hg38.mask.tsv.gz', sep='\t')
MASK = set(mk.loc[mk.M_general.astype(str).str.upper() == 'TRUE', 'Probe_ID'])
SHARED = set(pd.read_csv(f'{REF}/probe_annotation_shared.tsv', sep='\t', usecols=['probeID']).probeID)
pann = pd.read_csv(f'{P}/probe_annotation_hm450_promoter.tsv', sep='\t', index_col='probeID', low_memory=False)
pann = pann[~pann.index.isin(MASK)]
pur_all = pd.read_csv(f'{REF}/ABSOLUTE_pancan_purity.txt', sep='\t')
pur_all = pur_all.assign(s15=pur_all['sample'].astype(str).str[:15]).drop_duplicates('s15').set_index('s15')['purity']
seq = pd.read_excel(f'{REF}/DataS1_samples.xlsx', 'Sequencing Statistics', header=37)
seq['case'] = seq['Tissue_Barcode'].astype(str).str[:15]
seq['log10_dedup_reads'] = np.log10(seq['Final_DeDup_reads'].astype(float))
seq['dup_fraction'] = seq['Picard_Duplicate_reads'].astype(float) / seq['Aligned_reads'].astype(float)
TECH = seq.groupby('case')[['log10_dedup_reads', 'TSS_Enrichment_Score', 'FRIP', 'dup_fraction']].mean()

def promoter_features(T, probe_filter=None):
    """Returns METH (genes x tumours), HCP flag, ACC (genes x tumours), EXP (expressed genes x tumours), cohort."""
    atac = pd.read_csv(f'{P}/atac_{T}.tsv.gz', sep='\t')
    rna = pd.read_csv(f'{P}/rna_{T}.tsv.gz', sep='\t', index_col=0)
    parts = {os.path.basename(f).split('__')[1][:15]: f for f in glob.glob(f'{P}/meth_parts/{T}__*.tsv.gz')}
    cohort = sorted(set(atac.columns[5:]) & set(parts) & set(rna.columns))
    beta = pd.concat({s: pd.read_csv(parts[s], sep='\t', header=None, index_col=0)[1] for s in cohort}, axis=1)
    beta = beta.apply(pd.to_numeric, errors='coerce')
    ann = pann if probe_filter is None else pann[pann.index.isin(probe_filter)]
    beta = beta.reindex(ann.index)
    gg = genes[genes.chrom.isin(set(atac.seqnames))]
    acc = {}
    for c, g in gg.groupby('chrom'):
        sub = atac[atac.seqnames == c]; s, e = sub.start.values, sub.end.values; V = sub[cohort].values
        for gn, row in g.iterrows():
            hit = (e >= row.w0) & (s <= row.w1)
            if hit.any(): acc[gn] = V[hit].max(axis=0)
    ACC = pd.DataFrame(acc, index=cohort).T
    pc = ann['CpG_chrm'].astype(str).values; pb = ann['CpG_beg'].astype(float).values
    isl = (ann['CGIposition'].astype(str) == 'Island').values; Mall = beta[cohort].values
    mf, hcp = {}, {}
    for c, g in gg.groupby('chrom'):
        sel = np.where(pc == c)[0]
        if not len(sel): continue
        pos = pb[sel]
        for gn, row in g.iterrows():
            k = sel[(pos >= row.w0) & (pos <= row.w1)]
            if len(k) < 2: continue
            sub = Mall[k]; nn = np.sum(~np.isnan(sub), axis=0)
            mf[gn] = np.where(nn >= 2, np.nanmean(sub, axis=0), np.nan); hcp[gn] = bool(isl[k].any())
    METH = pd.DataFrame(mf, index=cohort).T
    expressed = rna.index[((rna[cohort] >= 1).mean(axis=1) >= 0.2)]
    EXP = rna.loc[expressed, cohort]
    return METH, pd.Series(hcp, dtype=bool), ACC, EXP, cohort, beta

def h2_matrices(METH, ACC, EXP):
    g2 = METH.index.intersection(ACC.index).intersection(EXP.index)
    M2, A2, E2 = METH.loc[g2], ACC.loc[g2], EXP.loc[g2]
    ok = ~(M2.isna().any(axis=1) | A2.isna().any(axis=1) | E2.isna().any(axis=1))
    M2, A2, E2 = M2[ok], A2[ok], E2[ok]
    var = (M2.std(axis=1) >= 0.05) & (A2.std(axis=1) > A2.std(axis=1).median()) & (E2.std(axis=1) >= 0.5)
    gV = M2.index[var.values]
    return gV, M2.loc[gV].values, A2.loc[gV].values, E2.loc[gV].values, int(ok.sum())

def analyse_type(T, full=True):
    METH, HCP, ACC, EXP, cohort, beta = promoter_features(T)
    res = {'n_tumours': len(cohort), 'genes': {'with_promoter_peak': int(len(ACC)), 'with_>=2_probes': int(len(METH)),
                                               'expressed': int(len(EXP))}}
    # S1: coupling
    g1 = METH.index.intersection(ACC.index); g1 = g1[(METH.loc[g1].std(axis=1) >= 0.05).values]
    Mv, Av = METH.loc[g1].values, ACC.loc[g1].values; ok = ~np.isnan(Mv).any(1)
    rho = spearman_rows(Mv[ok], Av[ok])
    res['S1'] = {'n_genes': int(ok.sum()), 'median_rho': float(np.median(rho)), 'frac_rho_le_-0.3': float(np.mean(rho <= -0.3))}
    # H2 per type
    gV, M, A, E, ncc = h2_matrices(METH, ACC, EXP)
    hcp = HCP.reindex(gV).fillna(False).values.astype(bool)
    res['genes'].update({'complete_cases': ncc, 'variable': int(len(gV)), 'variable_HCP': int(hcp.sum())})
    if len(gV) < 20:
        return res, None
    m, a, ma = models(M, A, E)
    res['H2'] = jackknife(M, A, E)
    pd.DataFrame({'gene': gV, 'HCP': hcp, 'cvR2_M': m, 'cvR2_A': a, 'cvR2_MA': ma}).to_csv(f'{R}/per_gene_{T}.csv', index=False)
    if not full:
        return res, (m, a, ma)
    res['S3_HCP'] = jackknife(M, A, E, sel=hcp) if hcp.sum() >= 20 else None
    res['S3_LCP'] = jackknife(M, A, E, sel=~hcp) if (~hcp).sum() >= 20 else None
    pur = pur_all.reindex(cohort).values.astype(float); okp = ~np.isnan(pur)
    res['S4_purity'] = jackknife(M[:, okp], A[:, okp], E[:, okp], [pur[okp]]) if okp.sum() >= 10 else None
    if res['S4_purity']: res['S4_purity']['n_tumours'] = int(okp.sum())
    res['mean_purity'] = float(np.nanmean(pur))
    # NC1 technical-only
    tech = TECH.reindex(cohort).copy(); tech['meth_sample_mean_beta'] = beta[cohort].mean(axis=0).values
    tech = tech.loc[:, tech.notna().all()]
    cols = [((tech[c] - tech[c].mean()) / (tech[c].std() + 1e-12)).values for c in tech.columns]
    G, n = E.shape; t = cv_r2(design(cols, G, n), E)
    res['NC1_technical_only'] = {'covariates': list(tech.columns), 'median_cvR2': float(np.median(t)),
                                 'frac_genes_MA_gt_technical': float(np.mean(ma > t)), 'MA_exceeds': bool(np.median(ma) > np.median(t))}
    # NC2 mismatched promoters within expression-mean decile
    dec = pd.qcut(E.mean(1), 10, labels=False, duplicates='drop'); null = []
    for d in range(100):
        perm = np.arange(len(gV))
        for q in np.unique(dec):
            ii = np.where(dec == q)[0]; perm[ii] = rng.permutation(ii)
        null.append(np.median(models(M[perm], A[perm], E)[2]))
    res['NC2_mismatched'] = {'null_mean_median_cvR2_MA': float(np.mean(null)),
                             'null_95pct': [float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
                             'observed_median_cvR2_MA': float(np.median(ma))}
    # S7 probe-set sensitivity: breast probe set (HM27 ∩ HM450)
    METHs, _, _, _, _, _ = promoter_features(T, probe_filter=SHARED)
    gVs, Ms, As, Es, _ = h2_matrices(METHs, ACC, EXP)
    res['S7_shared_probes'] = jackknife(Ms, As, Es) if len(gVs) >= 20 else None
    if res['S7_shared_probes']: res['S7_shared_probes']['n_genes'] = int(len(gVs))
    return res, (m, a, ma)

out = {'plan_sha256': 'c440f092afb11f7e4afb9a7204ef51550ff12d0d88ad6f5c765527fbd542fa3b', 'eligible_types': ELIGIBLE, 'per_type': {}}
types = sorted({os.path.basename(f)[4:-7] for f in glob.glob(f'{P}/rna_*.tsv.gz')})
for T in types:
    r, _ = analyse_type(T, full=T in ELIGIBLE)
    out['per_type'][T] = r
    print(T, r['n_tumours'], r['genes'].get('variable'), flush=True)

# ---------------- R1-R3 ----------------
inc = [T for T in ELIGIBLE if out['per_type'][T]['genes'].get('variable', 0) >= 100]
out['meta_types_included'] = inc
out['meta_types_excluded_lt100_genes'] = [T for T in ELIGIBLE if T not in inc]
D = [out['per_type'][T]['H2']['D_median_A_minus_M'] for T in inc]
Gm = [out['per_type'][T]['H2']['G_median_MA_minus_A'] for T in inc]
out['R1_pooled_D'] = dersimonian_laird([d['estimate'] for d in D], [d['SE'] for d in D])
out['R1_supported'] = bool(out['R1_pooled_D']['estimate'] > 0 and out['R1_pooled_D']['CI95'][0] > 0)
npos = int(sum(d['estimate'] > 0 for d in D))
out['R2_sign'] = {'types_D_positive': npos, 'k': len(D), 'one_sided_p': float(stats.binomtest(npos, len(D), 0.5, alternative='greater').pvalue),
                  'criterion_met_8_of_9': bool(npos >= 8 and len(D) == 9) if len(D) == 9 else 'k != 9; see deviations'}
out['R3_pooled_G'] = dersimonian_laird([g['estimate'] for g in Gm], [g['SE'] for g in Gm])
out['R3_supported'] = bool(out['R3_pooled_G']['CI95'][1] < 0.02)
for lab in ['S3_HCP', 'S3_LCP', 'S4_purity', 'S7_shared_probes']:
    ok = [T for T in inc if out['per_type'][T].get(lab)]
    if len(ok) >= 3:
        out[lab + '_pooled'] = {k: dersimonian_laird([out['per_type'][T][lab][k]['estimate'] for T in ok],
                                                     [out['per_type'][T][lab][k]['SE'] for T in ok]) for k in ['D_median_A_minus_M', 'G_median_MA_minus_A']}
        out[lab + '_pooled']['types'] = ok
# S1 pooled median rho (types weighted by tumour count; descriptive)
s1 = [(out['per_type'][T]['S1']['median_rho'], out['per_type'][T]['n_tumours']) for T in inc]
out['S1_median_rho_by_type'] = {T: out['per_type'][T]['S1']['median_rho'] for T in inc}
out['S1_weighted_mean_of_median_rho'] = float(np.average([x for x, _ in s1], weights=[w for _, w in s1]))
out['S6_heterogeneity'] = {'D_vs_n': [[T, out['per_type'][T]['n_tumours'], out['per_type'][T]['H2']['D_median_A_minus_M']['estimate']] for T in inc],
                           'D_vs_mean_purity': [[T, out['per_type'][T]['mean_purity'], out['per_type'][T]['H2']['D_median_A_minus_M']['estimate']] for T in inc]}

# ---------------- S2: all non-breast tumours pooled, within-type ----------------
feats = {}
for T in types:
    METH, HCP, ACC, EXP, cohort, _ = promoter_features(T)
    if len(cohort) >= 5: feats[T] = (METH, ACC, EXP, cohort)
common = None
for T, (METH, ACC, EXP, cohort) in feats.items():
    g = METH.index.intersection(ACC.index).intersection(EXP.index)
    g = g[~(METH.loc[g].isna().any(axis=1) | ACC.loc[g].isna().any(axis=1) | EXP.loc[g].isna().any(axis=1)).values]
    common = g if common is None else common.intersection(g)
cent = lambda X: X - X.mean(axis=1).values[:, None]
Mc = pd.concat([cent(feats[T][0].loc[common]) for T in feats], axis=1)
Ac = pd.concat([cent(feats[T][1].loc[common]) for T in feats], axis=1)
Ec = pd.concat([cent(feats[T][2].loc[common]) for T in feats], axis=1)
var = (Mc.std(axis=1) >= 0.05) & (Ac.std(axis=1) > Ac.std(axis=1).median()) & (Ec.std(axis=1) >= 0.5)
gS = common[var.values]
labels = np.concatenate([[T] * len(feats[T][3]) for T in feats])
def ztype(frames):  # z-score each gene within each cancer type, then join the types
    return np.concatenate([zs(f.loc[gS].values) for f in frames], axis=1)
Mz = ztype([feats[T][0][feats[T][3]] for T in feats]); Az = ztype([feats[T][1][feats[T][3]] for T in feats])
Ev = np.concatenate([feats[T][2].loc[gS, feats[T][3]].values for T in feats], axis=1)
dummies = [(labels == T).astype(float) for T in list(feats)[1:]]
def s2_stats(Mz, Az, Ev, dummies):
    G, n = Ev.shape
    base = press(design(dummies, G, n), Ev)
    rel = lambda cols: 1 - press(design(cols + dummies, G, n), Ev) / base
    m, a, ma = rel([Mz]), rel([Az]), rel([Mz, Az])
    return summ(m, a, ma)
full = s2_stats(Mz, Az, Ev, dummies); n = Ev.shape[1]; reps = []
for i in range(n):
    k = np.arange(n) != i
    reps.append(s2_stats(Mz[:, k], Az[:, k], Ev[:, k], [d[k] for d in dummies]))
s2 = {'types': list(feats), 'n_tumours': int(n), 'n_genes': int(len(gS))}
for key in full:
    v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
    s2[key] = {'estimate': float(full[key]), 'SE': float(se), 'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
out['S2_pooled_within_type'] = s2

json.dump(out, open(f'{R}/pancan_results.json', 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ['meta_types_included', 'R1_pooled_D', 'R1_supported', 'R2_sign', 'R3_pooled_G', 'R3_supported']}, indent=1))
