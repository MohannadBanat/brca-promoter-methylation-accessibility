"""POST HOC (deviations log entry 5): corrections after an internal pre-submission review.
 (g) CV-R2 relative to the leave-one-out covariate-only model, 1 - PRESS/PRESS0 (intercept only for non-breast types;
     intercept + array platform for breast); per type, jackknife CIs, pooled D* and G*.
 (h) S1 as specified in the plan: DerSimonian-Laird pooled median Spearman rho with delete-one-tumor jackknife SEs.
 (i) Hartung-Knapp-Sidik-Jonkman CIs and 95% prediction intervals for pooled D and G (R1, R3).
 (j) Recurrence comparator: the random-gene null applied to accessibility-informative genes (CV-R2_A > 0.3).
Run: python3 -I P08_review_corrections.py <processed_pancan> <ref> <results> <breast processed dir> <breast per-gene csv>"""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
P, REF, R, BP, BPG = sys.argv[1:6]
src = open(__file__.replace('P08_review_corrections.py', 'P05_posthoc_robustness.py')).read().split("out = {'note'")[0]
sys.argv = ['P05', P, REF, R, '']; ns = {}; exec(compile(src, 'P05_posthoc_robustness.py', 'exec'), ns)
ELIGIBLE, features, variable, zs, design, dl = ns['ELIGIBLE'], ns['features'], ns['variable'], ns['zs'], ns['design'], ns['dl']

def press(Xs, y):
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return (e ** 2).sum(1)
def rel_models(M, A, E, cov=()):
    G, n = E.shape; cov = list(cov); p0 = press(design(cov, G, n), E)
    f = lambda cols: 1 - press(design(cols + cov, G, n), E) / p0
    return f([zs(M)]), f([zs(A)]), f([zs(M), zs(A)])
def summ(m, a, ma):
    return {'median_rel_M': np.median(m), 'median_rel_A': np.median(a), 'median_rel_MA': np.median(ma),
            'D_rel': np.median(a - m), 'G_rel': np.median(ma - a)}
def jk(M, A, E, cov=()):
    n = E.shape[1]; full = summ(*rel_models(M, A, E, cov)); reps = []
    for i in range(n):
        k = np.arange(n) != i; reps.append(summ(*rel_models(M[:, k], A[:, k], E[:, k], [c[k] for c in cov])))
    out = {}
    for key in full:
        v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
        out[key] = {'estimate': float(full[key]), 'SE': float(se), 'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
    return out
def spearman_rows(X, Y):
    rx = stats.rankdata(X, axis=1); ry = stats.rankdata(Y, axis=1)
    rx -= rx.mean(1, keepdims=True); ry -= ry.mean(1, keepdims=True)
    return (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
def hksj(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float); k = len(est); d = dl(est, se)
    tau2 = max(0.0, (np.sum(est ** 2 / se ** 2) - np.sum(est / se ** 2) ** 2 / np.sum(1 / se ** 2) - (k - 1)) /
               (np.sum(1 / se ** 2) - np.sum(1 / se ** 4) / np.sum(1 / se ** 2)))
    w = 1 / (se ** 2 + tau2); mu = np.sum(w * est) / np.sum(w)
    q = np.sum(w * (est - mu) ** 2) / (k - 1); sh = np.sqrt(q / np.sum(w)); t1 = stats.t.ppf(0.975, k - 1)
    t2 = stats.t.ppf(0.975, k - 2); pse = np.sqrt(tau2 + d['SE'] ** 2)
    return {'estimate': float(mu), 'tau2': float(tau2), 'HKSJ_CI95': [float(mu - t1 * sh), float(mu + t1 * sh)],
            'prediction_interval95': [float(mu - t2 * pse), float(mu + t2 * pse)]}

out = {'note': 'POST HOC (deviations log entry 5).', 'per_type': {}}
p = json.load(open(f'{R}/pancan_results.json'))
for T in ELIGIBLE:
    gc, M, A, E, npr, cohort = features(T); v = variable(M, A, E); n = E.shape[1]
    r = {'n_tumours': n, 'n_genes': int(v.sum()), 'loo_mean_baseline_on_plan_scale': float(1 - (n / (n - 1)) ** 2)}
    r['g_relative'] = jk(M[v], A[v], E[v])
    # (h) S1 jackknife SE of the median rho (genes with methylation SD >= 0.05 and a peak; complete cases)
    ok = (np.nanstd(M, axis=1, ddof=1) >= 0.05)
    Ms, As = M[ok], A[ok]; full = float(np.median(spearman_rows(Ms, As))); reps = []
    for i in range(n):
        k = np.arange(n) != i; reps.append(np.median(spearman_rows(Ms[:, k], As[:, k])))
    reps = np.array(reps); se = float(np.sqrt((n - 1) / n * np.sum((reps - reps.mean()) ** 2)))
    r['h_S1'] = {'n_genes_complete_cases': int(ok.sum()), 'median_rho': full, 'SE': se,
                 'P02_median_rho_all_available_genes': p['per_type'][T]['S1']['median_rho']}
    out['per_type'][T] = r
    print(T, n, {k: round(x['estimate'], 3) for k, x in r['g_relative'].items()}, 'rho', round(full, 3), round(se, 3), flush=True)

# breast: same definitions as P03 (masked probes, platform covariate), baseline = intercept + platform
genes = pd.read_csv(f'{BP}/genes_tss_gencode_v36.tsv', sep='\t')
atac = pd.read_csv(f'{BP}/atac_brca_promoter_peaks_by_tumour.tsv.gz', sep='\t')
meth = pd.read_csv(f'{BP}/meth_shared_promoter_probes.tsv.gz', sep='\t', index_col=0)
pann0 = pd.read_csv(f'{BP}/probe_annotation_shared.tsv', sep='\t', index_col='probeID', low_memory=False)
plat = pd.read_csv(f'{BP}/meth_platform.tsv', sep='\t', index_col=0)['platform']
rna = pd.read_csv(f'{BP}/rna_rsem_log2_cohort.tsv.gz', sep='\t', index_col=0)
mk = pd.read_csv(f'{REF}/HM450.hg38.mask.tsv.gz', sep='\t'); MASK = set(mk.loc[mk.M_general.astype(str).str.upper() == 'TRUE', 'Probe_ID'])
cohort = sorted(set(atac.columns[5:]) & set(meth.columns) & set(rna.columns)); is450 = plat.loc[cohort].str.contains('450').astype(float).values
genes = genes[genes.chrom.isin(set(atac.seqnames))].copy()
genes['w0'] = np.where(genes.strand == '+', genes.tss - 1000, genes.tss - 500); genes['w1'] = np.where(genes.strand == '+', genes.tss + 500, genes.tss + 1000)
genes = genes.drop_duplicates('gene_name').set_index('gene_name'); acc = {}
for c, g in genes.groupby('chrom'):
    sub = atac[atac.seqnames == c]; s, e = sub.start.values, sub.end.values; V = sub[cohort].values
    for gn, row in g.iterrows():
        hit = (e >= row.w0) & (s <= row.w1)
        if hit.any(): acc[gn] = V[hit].max(axis=0)
ACC = pd.DataFrame(acc, index=cohort).T; EXP = rna.loc[rna.index[((rna[cohort] >= 1).mean(axis=1) >= 0.2)], cohort]
pann = pann0.loc[meth.index]; pann = pann[~pann.index.isin(MASK)]
pc = pann['CpG_chrm'].astype(str).values; pb = pann['CpG_beg'].astype(float).values; Mall = meth.loc[pann.index, cohort].values; mf = {}
for c, g in genes.groupby('chrom'):
    sel = np.where(pc == c)[0]
    if not len(sel): continue
    pos = pb[sel]
    for gn, row in g.iterrows():
        k = sel[(pos >= row.w0) & (pos <= row.w1)]
        if len(k) < 2: continue
        sub = Mall[k]; nn = np.sum(~np.isnan(sub), axis=0); mf[gn] = np.where(nn >= 2, np.nanmean(sub, axis=0), np.nan)
METH = pd.DataFrame(mf, index=cohort).T
g2 = METH.index.intersection(ACC.index).intersection(EXP.index); M2, A2, E2 = METH.loc[g2], ACC.loc[g2], EXP.loc[g2]
ok = ~(M2.isna().any(axis=1) | A2.isna().any(axis=1) | E2.isna().any(axis=1)); M2, A2, E2 = M2[ok], A2[ok], E2[ok]
var = (M2.std(axis=1) >= 0.05) & (A2.std(axis=1) > A2.std(axis=1).median()) & (E2.std(axis=1) >= 0.5); gV = M2.index[var.values]
assert len(gV) == 1193, len(gV)
out['BRCA'] = {'n_tumours': len(cohort), 'n_genes': int(len(gV)), 'g_relative': jk(M2.loc[gV].values, A2.loc[gV].values, E2.loc[gV].values, [is450])}
print('BRCA', {k: round(x['estimate'], 3) for k, x in out['BRCA']['g_relative'].items()})

pt = out['per_type']
out['g_pooled'] = {k: dl([pt[T]['g_relative'][k]['estimate'] for T in ELIGIBLE], [pt[T]['g_relative'][k]['SE'] for T in ELIGIBLE]) for k in ['D_rel', 'G_rel']}
out['g_types_median_rel_M_gt_0'] = [T for T in ELIGIBLE if pt[T]['g_relative']['median_rel_M']['estimate'] > 0]
out['h_S1_pooled_DL'] = dl([pt[T]['h_S1']['median_rho'] for T in ELIGIBLE], [pt[T]['h_S1']['SE'] for T in ELIGIBLE])
out['i_R1_D'] = hksj([p['per_type'][T]['H2']['D_median_A_minus_M']['estimate'] for T in ELIGIBLE], [p['per_type'][T]['H2']['D_median_A_minus_M']['SE'] for T in ELIGIBLE])
out['i_R3_G'] = hksj([p['per_type'][T]['H2']['G_median_MA_minus_A']['estimate'] for T in ELIGIBLE], [p['per_type'][T]['H2']['G_median_MA_minus_A']['SE'] for T in ELIGIBLE])
# (j) accessibility comparator for the recurrence analysis
pg = {T: pd.read_csv(f'{R}/per_gene_{T}.csv') for T in ELIGIBLE}; pg['BRCA'] = pd.read_csv(BPG)
rng = np.random.default_rng(20261126); res = {}
for lab, col in [('methylation', 'cvR2_M'), ('accessibility', 'cvR2_A')]:
    inf = {T: set(d.gene[d[col] > 0.3].astype(str)) for T, d in pg.items()}; tested = {T: d.gene.astype(str).values for T, d in pg.items()}
    cnt = pd.Series([g for T in inf for g in inf[T]]).value_counts(); nul = []
    for _ in range(1000):
        c = pd.Series(np.concatenate([rng.choice(tested[T], len(inf[T]), replace=False) for T in pg])).value_counts(); nul.append(int((c >= 5).sum()))
    res[lab] = {'genes_informative_in_>=5_of_10': int((cnt >= 5).sum()), 'null_mean': float(np.mean(nul)), 'null_max': int(np.max(nul)),
                'ratio_observed_to_null_mean': float((cnt >= 5).sum() / max(np.mean(nul), 1e-9))}
out['j_recurrence_comparator'] = res
json.dump(out, open(f'{R}/posthoc_review_corrections.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ['g_pooled', 'g_types_median_rel_M_gt_0', 'h_S1_pooled_DL', 'i_R1_D', 'i_R3_G', 'j_recurrence_comparator']}, indent=1))
