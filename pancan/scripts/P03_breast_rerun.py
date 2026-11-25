"""S5: breast (hypothesis-generating cohort) re-run through the pan-cancer statistics.
Same breast inputs, probe set (HM27 ∩ HM450) and platform covariate as the published analysis; reported with and without the
explicit M_general probe mask (pan-cancer deviation 1), so any difference from the published values is visible.
Run: python3 -I P03_breast_rerun.py <breast processed dir> <ref dir> <results dir>"""
import json, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore', message='Mean of empty slice')
BP, REF, R = sys.argv[1:4]

def cv_r2(Xs, y):
    XtX = Xs.transpose(0, 2, 1) @ Xs + 1e-9 * np.eye(Xs.shape[2]); Ai = np.linalg.inv(XtX)
    b = Ai @ (Xs.transpose(0, 2, 1) @ y[..., None]); res = y - (Xs @ b)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ai, Xs); e = res / np.clip(1 - h, 1e-6, None)
    return 1 - (e ** 2).sum(1) / ((y - y.mean(1, keepdims=True)) ** 2).sum(1)
def zs(X): return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)
def design(cols, G, n): return np.stack([np.ones((G, n))] + [np.broadcast_to(c, (G, n)) for c in cols], axis=2)
def models(M, A, E, cov):
    G, n = E.shape
    return (cv_r2(design([zs(M)] + cov, G, n), E), cv_r2(design([zs(A)] + cov, G, n), E), cv_r2(design([zs(M), zs(A)] + cov, G, n), E))
def summ(m, a, ma):
    return {'D_median_A_minus_M': np.median(a - m), 'G_median_MA_minus_A': np.median(ma - a), 'median_cvR2_M': np.median(m),
            'median_cvR2_A': np.median(a), 'median_cvR2_MA': np.median(ma), 'median_MA_minus_M': np.median(ma - m), 'frac_A_gt_M': np.mean(a > m)}
def jackknife(M, A, E, cov):
    n = E.shape[1]; full = summ(*models(M, A, E, cov)); reps = []
    for i in range(n):
        k = np.arange(n) != i; reps.append(summ(*models(M[:, k], A[:, k], E[:, k], [c[k] for c in cov])))
    out = {}
    for key in full:
        v = np.array([r[key] for r in reps]); se = np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))
        out[key] = {'estimate': float(full[key]), 'SE': float(se), 'CI95': [float(full[key] - 1.96 * se), float(full[key] + 1.96 * se)]}
    return out

genes = pd.read_csv(f'{BP}/genes_tss_gencode_v36.tsv', sep='\t')
atac = pd.read_csv(f'{BP}/atac_brca_promoter_peaks_by_tumour.tsv.gz', sep='\t')
meth = pd.read_csv(f'{BP}/meth_shared_promoter_probes.tsv.gz', sep='\t', index_col=0)
pann0 = pd.read_csv(f'{BP}/probe_annotation_shared.tsv', sep='\t', index_col='probeID', low_memory=False)
plat = pd.read_csv(f'{BP}/meth_platform.tsv', sep='\t', index_col=0)['platform']
rna = pd.read_csv(f'{BP}/rna_rsem_log2_cohort.tsv.gz', sep='\t', index_col=0)
mk = pd.read_csv(f'{REF}/HM450.hg38.mask.tsv.gz', sep='\t')
MASK = set(mk.loc[mk.M_general.astype(str).str.upper() == 'TRUE', 'Probe_ID'])
cohort = sorted(set(atac.columns[5:]) & set(meth.columns) & set(rna.columns))
is450 = plat.loc[cohort].str.contains('450').astype(float).values
genes = genes[genes.chrom.isin(set(atac.seqnames))].copy()
genes['w0'] = np.where(genes.strand == '+', genes.tss - 1000, genes.tss - 500)
genes['w1'] = np.where(genes.strand == '+', genes.tss + 500, genes.tss + 1000)
genes = genes.drop_duplicates('gene_name').set_index('gene_name')
acc = {}
for c, g in genes.groupby('chrom'):
    sub = atac[atac.seqnames == c]; s, e = sub.start.values, sub.end.values; V = sub[cohort].values
    for gn, row in g.iterrows():
        hit = (e >= row.w0) & (s <= row.w1)
        if hit.any(): acc[gn] = V[hit].max(axis=0)
ACC = pd.DataFrame(acc, index=cohort).T
EXP = rna.loc[rna.index[((rna[cohort] >= 1).mean(axis=1) >= 0.2)], cohort]

def run(mask):
    pann = pann0.loc[meth.index]
    if mask: pann = pann[~pann.index.isin(MASK)]
    pc = pann['CpG_chrm'].astype(str).values; pb = pann['CpG_beg'].astype(float).values; Mall = meth.loc[pann.index, cohort].values
    mf = {}
    for c, g in genes.groupby('chrom'):
        sel = np.where(pc == c)[0]
        if not len(sel): continue
        pos = pb[sel]
        for gn, row in g.iterrows():
            k = sel[(pos >= row.w0) & (pos <= row.w1)]
            if len(k) < 2: continue
            sub = Mall[k]; nn = np.sum(~np.isnan(sub), axis=0)
            mf[gn] = np.where(nn >= 2, np.nanmean(sub, axis=0), np.nan)
    METH = pd.DataFrame(mf, index=cohort).T
    g2 = METH.index.intersection(ACC.index).intersection(EXP.index)
    M2, A2, E2 = METH.loc[g2], ACC.loc[g2], EXP.loc[g2]
    ok = ~(M2.isna().any(axis=1) | A2.isna().any(axis=1) | E2.isna().any(axis=1)); M2, A2, E2 = M2[ok], A2[ok], E2[ok]
    var = (M2.std(axis=1) >= 0.05) & (A2.std(axis=1) > A2.std(axis=1).median()) & (E2.std(axis=1) >= 0.5)
    gV = M2.index[var.values]
    r = jackknife(M2.loc[gV].values, A2.loc[gV].values, E2.loc[gV].values, [is450]); r['n_genes'] = int(len(gV)); r['n_tumours'] = len(cohort)
    return r

out = {'unmasked_as_published': run(False), 'masked_M_general': run(True)}
json.dump(out, open(f'{R}/breast_rerun.json', 'w'), indent=1)
for k, v in out.items():
    print(k, v['n_genes'], {kk: round(vv['estimate'], 4) for kk, vv in v.items() if isinstance(vv, dict)})
