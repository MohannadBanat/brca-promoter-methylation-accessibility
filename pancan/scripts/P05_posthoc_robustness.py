"""POST HOC, EXPLORATORY analyses requested by an internal referee-style review (pancan_deviations_log.md, entry 3).
Written after the primary results were viewed; none of these analyses was in the plan. Nine eligible types only.
 (a) gene selection repeated inside every jackknife replicate (D, G and median CV-R2 of M, A, MA)
 (b) strata by number of unmasked HM450 probes in the promoter window (2-3, 4-7, >=8)
 (c) genes where methylation is informative (CV-R2_M > 0.3) and their recurrence across types (uses per_gene_*.csv and breast Table S2)
 (d) two-state methylation model: tumours split at the Otsu threshold of each gene's beta values (both groups >= 3 tumours)
 (e) ATAC-seq quality (mean TSS enrichment, FRiP) against D and median CV-R2_A by type (descriptive)
Run: python3 -I P05_posthoc_robustness.py <processed_pancan dir> <ref dir> <results dir> <breast per-gene csv>
"""
import glob, json, os, sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore', message='Mean of empty slice')
warnings.filterwarnings('ignore', message='Degrees of freedom')
P, REF, R, BREAST_PG = sys.argv[1:5]
ELIGIBLE = ['KIRP', 'PRAD', 'COAD', 'LUAD', 'STAD', 'ESCA', 'LIHC', 'LUSC', 'KIRC']

def cv_r2(Xs, y):
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return 1 - (e ** 2).sum(1) / ((y - y.mean(1, keepdims=True)) ** 2).sum(1)
def zs(X): return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)
def design(cols, G, n): return np.stack([np.ones((G, n))] + [np.broadcast_to(c, (G, n)) for c in cols], axis=2)
def models(M, A, E):
    G, n = E.shape
    return cv_r2(design([zs(M)], G, n), E), cv_r2(design([zs(A)], G, n), E), cv_r2(design([zs(M), zs(A)], G, n), E)
def summ(m, a, ma):
    return {'D_median_A_minus_M': np.median(a - m), 'G_median_MA_minus_A': np.median(ma - a), 'median_cvR2_M': np.median(m),
            'median_cvR2_A': np.median(a), 'median_cvR2_MA': np.median(ma)}
def jk_out(full, reps, n):
    out = {}
    for key in full:
        v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
        out[key] = {'estimate': float(full[key]), 'SE': float(se), 'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
    return out
def dl(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float); w = 1 / se ** 2
    fe = np.sum(w * est) / np.sum(w); Q = np.sum(w * (est - fe) ** 2); k = len(est)
    tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    wr = 1 / (se ** 2 + tau2); re = np.sum(wr * est) / np.sum(wr); sre = np.sqrt(1 / np.sum(wr))
    return {'estimate': float(re), 'SE': float(sre), 'CI95': [float(re - 1.96 * sre), float(re + 1.96 * sre)], 'k': int(k),
            'I2': float(max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)}
def variable(M, A, E):
    sa = A.std(axis=1, ddof=1)  # sample SD, as pandas in P02
    return (M.std(axis=1, ddof=1) >= 0.05) & (sa > np.median(sa)) & (E.std(axis=1, ddof=1) >= 0.5)

genes = pd.read_csv(f'{REF}/genes_tss_gencode_v36.tsv', sep='\t')
genes['w0'] = np.where(genes.strand == '+', genes.tss - 1000, genes.tss - 500)
genes['w1'] = np.where(genes.strand == '+', genes.tss + 500, genes.tss + 1000)
genes = genes.drop_duplicates('gene_name').set_index('gene_name')
mk = pd.read_csv(f'{REF}/HM450.hg38.mask.tsv.gz', sep='\t')
MASK = set(mk.loc[mk.M_general.astype(str).str.upper() == 'TRUE', 'Probe_ID'])
pann = pd.read_csv(f'{P}/probe_annotation_hm450_promoter.tsv', sep='\t', index_col='probeID', low_memory=False)
pann = pann[~pann.index.isin(MASK)]
seq = pd.read_excel(f'{REF}/DataS1_samples.xlsx', 'Sequencing Statistics', header=37)
seq['case'] = seq['Tissue_Barcode'].astype(str).str[:15]
TECH = seq.groupby('case')[['TSS_Enrichment_Score', 'FRIP']].mean()

def features(T):
    """Same definitions as P02; also returns the number of non-missing probes per gene (median across tumours)."""
    atac = pd.read_csv(f'{P}/atac_{T}.tsv.gz', sep='\t')
    rna = pd.read_csv(f'{P}/rna_{T}.tsv.gz', sep='\t', index_col=0)
    parts = {os.path.basename(f).split('__')[1][:15]: f for f in glob.glob(f'{P}/meth_parts/{T}__*.tsv.gz')}
    cohort = sorted(set(atac.columns[5:]) & set(parts) & set(rna.columns))
    beta = pd.concat({s: pd.read_csv(parts[s], sep='\t', header=None, index_col=0)[1] for s in cohort}, axis=1)
    beta = beta.apply(pd.to_numeric, errors='coerce').reindex(pann.index)
    gg = genes[genes.chrom.isin(set(atac.seqnames))]
    acc = {}
    for c, g in gg.groupby('chrom'):
        sub = atac[atac.seqnames == c]; s, e = sub.start.values, sub.end.values; V = sub[cohort].values
        for gn, row in g.iterrows():
            hit = (e >= row.w0) & (s <= row.w1)
            if hit.any(): acc[gn] = V[hit].max(axis=0)
    ACC = pd.DataFrame(acc, index=cohort).T
    pc = pann['CpG_chrm'].astype(str).values; pb = pann['CpG_beg'].astype(float).values; Mall = beta[cohort].values
    mf, npr = {}, {}
    for c, g in gg.groupby('chrom'):
        sel = np.where(pc == c)[0]
        if not len(sel): continue
        pos = pb[sel]
        for gn, row in g.iterrows():
            k = sel[(pos >= row.w0) & (pos <= row.w1)]
            if len(k) < 2: continue
            sub = Mall[k]; nn = np.sum(~np.isnan(sub), axis=0)
            mf[gn] = np.where(nn >= 2, np.nanmean(sub, axis=0), np.nan); npr[gn] = float(np.median(nn))
    METH = pd.DataFrame(mf, index=cohort).T
    EXP = rna.loc[rna.index[((rna[cohort] >= 1).mean(axis=1) >= 0.2)], cohort]
    g2 = METH.index.intersection(ACC.index).intersection(EXP.index)
    M2, A2, E2 = METH.loc[g2], ACC.loc[g2], EXP.loc[g2]
    ok = ~(M2.isna().any(axis=1) | A2.isna().any(axis=1) | E2.isna().any(axis=1))
    gc = g2[ok.values]
    return gc, METH.loc[gc].values, ACC.loc[gc].values, EXP.loc[gc].values, pd.Series(npr).reindex(gc).values, cohort

def otsu_indicator(M, min_group=3):
    """Per gene: threshold maximising between-group variance of beta; returns indicator matrix and validity mask."""
    G, n = M.shape; S = np.sort(M, axis=1); I = np.zeros_like(M); valid = np.zeros(G, bool)
    for g in range(G):
        best, thr = -1, None
        for j in range(min_group, n - min_group + 1):
            lo, hi = S[g, :j], S[g, j:]
            if S[g, j - 1] == S[g, j]: continue
            bv = j * (n - j) * (hi.mean() - lo.mean()) ** 2
            if bv > best: best, thr = bv, (S[g, j - 1] + S[g, j]) / 2
        if thr is not None:
            I[g] = (M[g] > thr).astype(float); valid[g] = True
    return I, valid

out = {'note': 'POST HOC, EXPLORATORY (deviation 3). Not part of the committed plan.', 'per_type': {}}
for T in ELIGIBLE:
    gc, M, A, E, npr, cohort = features(T); n = E.shape[1]; r = {'n_tumours': n}
    # (a) selection inside each jackknife replicate
    v = variable(M, A, E); full = summ(*models(M[v], A[v], E[v])); reps = []; ng = []
    for i in range(n):
        k = np.arange(n) != i; Mi, Ai, Ei = M[:, k], A[:, k], E[:, k]; vi = variable(Mi, Ai, Ei); ng.append(int(vi.sum()))
        reps.append(summ(*models(Mi[vi], Ai[vi], Ei[vi])))
    r['a_reselect_jackknife'] = jk_out(full, reps, n); r['a_reselect_jackknife']['n_genes_range'] = [min(ng), max(ng)]
    r['a_reselect_jackknife']['n_genes_full'] = int(v.sum())
    # (b) probe-count strata (variable genes, full-data selection as in P02)
    Mv, Av, Ev, nv = M[v], A[v], E[v], npr[v]; m, a, ma = models(Mv, Av, Ev)
    r['b_probe_strata'] = {}
    for lab, lo, hi in [('2-3', 2, 3.5), ('4-7', 3.5, 7.5), ('>=8', 7.5, 1e9)]:
        s = (nv >= lo) & (nv < hi)
        if s.sum() >= 20:
            r['b_probe_strata'][lab] = {'n_genes': int(s.sum()), **{k2: float(x) for k2, x in summ(m[s], a[s], ma[s]).items()}}
    # (d) two-state methylation
    I, valid = otsu_indicator(Mv)
    Gv = int(valid.sum()); Iv = I[valid]; Ev2 = Ev[valid]; Av2 = Av[valid]
    b = cv_r2(design([Iv], Gv, n), Ev2); ab = cv_r2(design([Iv, zs(Av2)], Gv, n), Ev2)
    r['d_two_state'] = {'n_genes': Gv, 'median_cvR2_two_state': float(np.median(b)), 'median_cvR2_linear_M': float(np.median(m[valid])),
                        'median_cvR2_A': float(np.median(a[valid])), 'median_gain_two_state_added_to_A': float(np.median(ab - a[valid])),
                        'frac_two_state_gt_linear': float(np.mean(b > m[valid]))}
    # (e) ATAC quality
    t = TECH.reindex(cohort)
    r['e_atac_quality'] = {'mean_TSS_enrichment': float(t.TSS_Enrichment_Score.mean()), 'mean_FRiP': float(t.FRIP.mean())}
    out['per_type'][T] = r
    print(T, n, int(v.sum()), {k2: round(x['estimate'], 3) for k2, x in r['a_reselect_jackknife'].items() if isinstance(x, dict)}, flush=True)

pt = out['per_type']
out['a_pooled'] = {k: dl([pt[T]['a_reselect_jackknife'][k]['estimate'] for T in ELIGIBLE], [pt[T]['a_reselect_jackknife'][k]['SE'] for T in ELIGIBLE])
                   for k in ['D_median_A_minus_M', 'G_median_MA_minus_A']}
out['a_types_D_positive'] = int(sum(pt[T]['a_reselect_jackknife']['D_median_A_minus_M']['estimate'] > 0 for T in ELIGIBLE))
out['a_types_D_CI_excludes_0'] = int(sum(pt[T]['a_reselect_jackknife']['D_median_A_minus_M']['CI95'][0] > 0 for T in ELIGIBLE))
# (c) methylation-informative genes and recurrence (P02 per-gene CV-R2, plus breast)
pg = {T: pd.read_csv(f'{R}/per_gene_{T}.csv') for T in ELIGIBLE}
br = pd.read_csv(BREAST_PG).rename(columns=str.strip); pg['BRCA'] = br
inf = {}
for T, d in pg.items():
    s = d[d.cvR2_M > 0.3]
    inf[T] = {'n_genes': int(len(d)), 'n_M_gt_0.3': int(len(s)), 'frac_M_gt_0.3': float(len(s) / len(d)),
              'n_M_beats_A_by_0.1': int(((d.cvR2_M - d.cvR2_A) > 0.1).sum()), 'genes': sorted(s.gene.astype(str))}
cnt = pd.Series([g for T in inf for g in inf[T]['genes']]).value_counts()
tested = pd.Series([g for T in pg for g in pg[T].gene.astype(str)]).value_counts()
rec = pd.DataFrame({'types_M_gt_0.3': cnt, 'types_tested': tested.reindex(cnt.index)}).sort_values(['types_M_gt_0.3', 'types_tested'], ascending=False)
rec.index.name = 'gene'; rec.to_csv(f'{R}/posthoc_methylation_informative_recurrence.csv')
out['c_methylation_informative'] = {T: {k: v for k, v in x.items() if k != 'genes'} for T, x in inf.items()}
rec['chrom'] = genes['chrom'].reindex(rec.index).values; rec.to_csv(f'{R}/posthoc_methylation_informative_recurrence.csv')
out['c_recurrent_in_>=3_types'] = rec[rec['types_M_gt_0.3'] >= 3].reset_index().to_dict('records')
# recurrence null: within each cohort, the same number of 'informative' genes drawn at random from its tested genes (1000 draws)
rng = np.random.default_rng(20261126); tested_by = {T: pg[T].gene.astype(str).values for T in pg}; nulls = {5: [], 3: []}
for _ in range(1000):
    c = pd.Series(np.concatenate([rng.choice(tested_by[T], inf[T]['n_M_gt_0.3'], replace=False) for T in pg])).value_counts()
    for th in nulls: nulls[th].append(int((c >= th).sum()))
out['c_recurrence_vs_null'] = {f'>={th}_cohorts': {'observed': int((cnt >= th).sum()), 'null_mean': float(np.mean(v)),
                               'null_max': int(np.max(v)), 'perm_p': float((1 + np.sum(np.array(v) >= (cnt >= th).sum())) / 1001)}
                               for th, v in nulls.items()}
print(out['c_recurrence_vs_null'])
# (e) correlation across the 9 types (descriptive, k = 9)
P2 = json.load(open(f'{R}/pancan_results.json'))['per_type']
tss = [pt[T]['e_atac_quality']['mean_TSS_enrichment'] for T in ELIGIBLE]
Dv = [P2[T]['H2']['D_median_A_minus_M']['estimate'] for T in ELIGIBLE]; Av_ = [P2[T]['H2']['median_cvR2_A']['estimate'] for T in ELIGIBLE]
out['e_spearman_TSS_enrichment_vs_D'] = float(stats.spearmanr(tss, Dv)[0])
out['e_spearman_TSS_enrichment_vs_median_cvR2_A'] = float(stats.spearmanr(tss, Av_)[0])
json.dump(out, open(f'{R}/posthoc_robustness.json', 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ['a_pooled', 'a_types_D_positive', 'a_types_D_CI_excludes_0', 'e_spearman_TSS_enrichment_vs_D',
                                      'e_spearman_TSS_enrichment_vs_median_cvR2_A']}, indent=1))
