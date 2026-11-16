"""Main analysis for the IJMS paper, following plan/analysis_plan_2026-11-16.md (SHA-256 b2b7079a...).
Run:  python3 -I scripts/04_analysis.py processed results 1000
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore', message='Mean of empty slice')  # genes with <2 probes in a tumour are set to NaN below
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

P, R = sys.argv[1], sys.argv[2]
NBOOT = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
os.makedirs(R, exist_ok=True)
rng = np.random.default_rng(20261116)
out = {}

# ---------------- load ----------------
genes = pd.read_csv(f'{P}/genes_tss_gencode_v36.tsv', sep='\t')
atac = pd.read_csv(f'{P}/atac_brca_promoter_peaks_by_tumour.tsv.gz', sep='\t')
meth = pd.read_csv(f'{P}/meth_shared_promoter_probes.tsv.gz', sep='\t', index_col=0)
pann = pd.read_csv(f'{P}/probe_annotation_shared.tsv', sep='\t', index_col='probeID', low_memory=False)
plat = pd.read_csv(f'{P}/meth_platform.tsv', sep='\t', index_col=0)['platform']
rna = pd.read_csv(f'{P}/rna_rsem_log2_cohort.tsv.gz', sep='\t', index_col=0)
clin = pd.read_csv(f'{P}/clinical_cohort.tsv', sep='\t', index_col=0, low_memory=False)
cohort = sorted(set(atac.columns[5:]) & set(meth.columns) & set(rna.columns))
out['n_tumours'] = len(cohort)
out['n_by_platform'] = plat.loc[cohort].value_counts().to_dict()
pam = clin['PAM50Call_RNAseq'].reindex(cohort)
out['pam50'] = pam.fillna('NA').value_counts().to_dict()
is450 = (plat.loc[cohort].str.contains('450')).astype(float).values

# ---------------- promoter features ----------------
genes = genes[genes.chrom.isin(set(atac.seqnames))].copy()
genes['w0'] = np.where(genes.strand == '+', genes.tss - 1000, genes.tss - 500)
genes['w1'] = np.where(genes.strand == '+', genes.tss + 500, genes.tss + 1000)
genes = genes.drop_duplicates('gene_name').set_index('gene_name')

A_mat = atac[cohort].values
acc = {}
for c, g in genes.groupby('chrom'):
    sub = atac[atac.seqnames == c]
    s, e = sub.start.values, sub.end.values
    idx = sub.index.values
    for gn, row in g.iterrows():
        hit = idx[(e >= row.w0) & (s <= row.w1)]
        if len(hit): acc[gn] = atac.loc[hit, cohort].values.max(axis=0)
ACC = pd.DataFrame(acc, index=cohort).T

pann = pann.loc[meth.index]
pc = pann['CpG_chrm'].astype(str).values; pb = pann['CpG_beg'].astype(float).values
isl = (pann['CGIposition'].astype(str) == 'Island').values if 'CGIposition' in pann else pann['CGI'].astype(str).str.contains('Island', na=False).values
M_all = meth[cohort].values
mfeat, hcp = {}, {}
for c, g in genes.groupby('chrom'):
    sel = np.where(pc == c)[0]
    if not len(sel): continue
    pos = pb[sel]
    for gn, row in g.iterrows():
        k = sel[(pos >= row.w0) & (pos <= row.w1)]
        if len(k) < 2: continue
        sub = M_all[k]; nn = np.sum(~np.isnan(sub), axis=0)
        v = np.where(nn >= 2, np.nanmean(np.where(np.isnan(sub), np.nan, sub), axis=0), np.nan)
        mfeat[gn] = v; hcp[gn] = bool(isl[k].any())
METH = pd.DataFrame(mfeat, index=cohort).T
HCP = pd.Series(hcp)

expressed = rna.index[((rna[cohort] >= 1).mean(axis=1) >= 0.2)]
EXP = rna.loc[expressed, cohort]
out['genes'] = {'protein_coding_tss': int(len(genes)), 'with_promoter_peak': int(len(ACC)),
                'with_>=2_probes': int(len(METH)), 'expressed': int(len(EXP))}

# ---------------- H1 ----------------
g1 = METH.index.intersection(ACC.index)
g1 = g1[(METH.loc[g1].std(axis=1) >= 0.05).values]
Mv, Av = METH.loc[g1].values, ACC.loc[g1].values
def spearman_rows(X, Y):
    rx = stats.rankdata(X, axis=1); ry = stats.rankdata(Y, axis=1)
    rx -= rx.mean(1, keepdims=True); ry -= ry.mean(1, keepdims=True)
    return (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
ok = ~np.isnan(Mv).any(1)
rho = spearman_rows(Mv[ok], Av[ok]); gH1 = g1[ok]
frac = float(np.mean(rho <= -0.3))
bf = []
for b in range(NBOOT):
    ii = rng.integers(0, len(cohort), len(cohort))
    r_b = spearman_rows(Mv[ok][:, ii], Av[ok][:, ii]); bf.append(np.nanmean(r_b <= -0.3))
h1 = {'n_genes': int(ok.sum()), 'median_rho': float(np.median(rho)),
      'frac_rho_le_-0.3': frac, 'frac_CI95': [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))]}
for lab, m in [('HCP', HCP.reindex(gH1).fillna(False).values), ('LCP', ~HCP.reindex(gH1).fillna(False).values)]:
    h1[lab] = {'n': int(m.sum()), 'median_rho': float(np.median(rho[m])) if m.any() else None,
               'frac_rho_le_-0.3': float(np.mean(rho[m] <= -0.3)) if m.any() else None}
# discordant states over gene-tumour pairs
t_hi = np.percentile(Av[ok], 200 / 3, axis=1, keepdims=True); t_lo = np.percentile(Av[ok], 100 / 3, axis=1, keepdims=True)
h1['pairs'] = int(Mv[ok].size)
h1['methylated_yet_accessible'] = int(((Mv[ok] >= 0.6) & (Av[ok] >= t_hi)).sum())
h1['unmethylated_yet_closed'] = int(((Mv[ok] <= 0.2) & (Av[ok] <= t_lo)).sum())
out['H1'] = h1
pd.DataFrame({'gene': gH1, 'rho': rho, 'HCP': HCP.reindex(gH1).values}).to_csv(f'{R}/H1_per_gene_rho.csv', index=False)

# ---------------- H2 ----------------
g2 = METH.index.intersection(ACC.index).intersection(EXP.index)
M2, A2, E2 = METH.loc[g2], ACC.loc[g2], EXP.loc[g2]
complete = ~(M2.isna().any(axis=1) | A2.isna().any(axis=1) | E2.isna().any(axis=1))
M2, A2, E2 = M2[complete], A2[complete], E2[complete]
sd_a_med = A2.std(axis=1).median()
var = (M2.std(axis=1) >= 0.05) & (A2.std(axis=1) > sd_a_med) & (E2.std(axis=1) >= 0.5)
gV = M2.index[var.values]
out['H2_genes'] = {'complete_cases': int(complete.sum()), 'variable': int(len(gV)),
                   'HCP': int(HCP.reindex(gV).fillna(False).sum())}
Mx, Ax, Ex = M2.loc[gV].values, A2.loc[gV].values, E2.loc[gV].values
_pur = pd.read_csv(f'{P}/purity_absolute.tsv', sep='\t', index_col=0)['purity'].reindex(cohort).values if os.path.exists(f'{P}/purity_absolute.tsv') else np.full(len(cohort), np.nan)
np.savez_compressed(f'{R}/H2_matrices.npz', M=Mx, A=Ax, E=Ex, is450=is450, purity=_pur, genes=np.array(gV), hcp=HCP.reindex(gV).fillna(False).values.astype(bool), cohort=np.array(cohort))

def cv_r2(Xs, y, w=None):
    """Weighted leave-one-unit-out R2 for a batch of genes. Xs: (G,n,k) design incl. intercept; y: (G,n); w: (n,) multiplicities."""
    G, n, k = Xs.shape
    w = np.ones(n) if w is None else w.astype(float)
    keep = w > 0; Xs, y, w = Xs[:, keep], y[:, keep], w[keep]
    XtW = Xs.transpose(0, 2, 1) * w
    A = XtW @ Xs + 1e-9 * np.eye(k)
    Ainv = np.linalg.inv(A)
    beta = Ainv @ (XtW @ y[..., None])
    res = y - (Xs @ beta)[..., 0]
    h = np.einsum('gnk,gkl,gnl->gn', Xs, Ainv, Xs)
    eloo = res / np.clip(1 - w * h, 1e-6, None)
    ybar = (y * w).sum(1, keepdims=True) / w.sum()
    return 1 - (w * eloo ** 2).sum(1) / (w * (y - ybar) ** 2).sum(1)

def design(cols, n, G):
    return np.stack([np.ones((G, n))] + cols, axis=2)

def zs(X): return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)

def h2_models(Mx, Ax, Ex, cov=None, w=None):
    G, n = Ex.shape; extra = [] if cov is None else [np.broadcast_to(c, (G, n)) for c in cov]
    m = cv_r2(design([zs(Mx)] + extra, n, G), Ex, w)
    a = cv_r2(design([zs(Ax)] + extra, n, G), Ex, w)
    ma = cv_r2(design([zs(Mx), zs(Ax)] + extra, n, G), Ex, w)
    base = cv_r2(design(extra, n, G), Ex, w) if extra else None
    return m, a, ma, base

def summarise(m, a, ma):
    return {'median_cvR2_M': float(np.median(m)), 'median_cvR2_A': float(np.median(a)), 'median_cvR2_MA': float(np.median(ma)),
            'median_MA_minus_A': float(np.median(ma - a)), 'median_MA_minus_M': float(np.median(ma - m)),
            'frac_A_gt_M': float(np.mean(a > m))}

def boot(fn, n, nb=NBOOT):
    keys, vals = None, []
    for b in range(nb):
        w = np.bincount(rng.integers(0, n, n), minlength=n)
        s = fn(w); keys = list(s); vals.append([s[k] for k in keys])
    v = np.array(vals)
    return {k: [float(np.percentile(v[:, i], 2.5)), float(np.percentile(v[:, i], 97.5))] for i, k in enumerate(keys)}

covP = [is450]
m, a, ma, base = h2_models(Mx, Ax, Ex, covP)
pooled = summarise(m - 0, a, ma); pooled['median_cvR2_platform_only'] = float(np.median(base))
pooled['CI95'] = boot(lambda w: summarise(*h2_models(Mx, Ax, Ex, covP, w)[:3]), len(cohort))
hcpV = HCP.reindex(gV).fillna(False).values
pooled['HCP'] = summarise(m[hcpV], a[hcpV], ma[hcpV]); pooled['LCP'] = summarise(m[~hcpV], a[~hcpV], ma[~hcpV])
out['H2_pooled'] = pooled
pd.DataFrame({'gene': gV, 'HCP': hcpV, 'cvR2_M': m, 'cvR2_A': a, 'cvR2_MA': ma}).to_csv(f'{R}/H2_per_gene_cvR2.csv', index=False)

# platform split
split = {}
for lab, sel in [('HM450', is450 == 1), ('HM27', is450 == 0)]:
    mm, aa, mama, _ = h2_models(Mx[:, sel], Ax[:, sel], Ex[:, sel])
    s = summarise(mm, aa, mama)
    s['n_tumours'] = int(sel.sum())
    s['CI95'] = boot(lambda w, sel=sel: summarise(*h2_models(Mx[:, sel], Ax[:, sel], Ex[:, sel], None, w)[:3]), int(sel.sum()))
    split[lab] = s
split['direction_agrees'] = bool(np.sign(split['HM450']['median_MA_minus_A']) == np.sign(split['HM27']['median_MA_minus_A'])
                                 and np.sign(split['HM450']['median_MA_minus_M']) == np.sign(split['HM27']['median_MA_minus_M']))
out['H2_platform_split'] = split

# purity sensitivity (ABSOLUTE if supplied)
pur_file = f'{P}/purity_absolute.tsv'
if os.path.exists(pur_file):
    pur = pd.read_csv(pur_file, sep='\t', index_col=0)['purity'].reindex(cohort)
    okp = ~pur.isna().values
    pv = pur.values[okp]
    mm, aa, mama, _ = h2_models(Mx[:, okp], Ax[:, okp], Ex[:, okp], [is450[okp], pv])
    s = summarise(mm, aa, mama); s['n_tumours'] = int(okp.sum())
    s['CI95'] = boot(lambda w: summarise(*h2_models(Mx[:, okp], Ax[:, okp], Ex[:, okp], [is450[okp], pv], w)[:3]), int(okp.sum()))
    out['H2_purity_adjusted'] = s
else:
    out['H2_purity_adjusted'] = 'purity file not available'

# negative control 1: technical-only model
tech_file = f'{P}/technical_covariates.tsv'
if os.path.exists(tech_file):
    tech = pd.read_csv(tech_file, sep='\t', index_col=0).reindex(cohort)
    tech = tech.loc[:, tech.notna().all()]
    cols = [((tech[c] - tech[c].mean()) / (tech[c].std() + 1e-12)).values for c in tech.columns] + [is450]
    G, n = Ex.shape
    t = cv_r2(design([np.broadcast_to(c, (G, n)) for c in cols], n, G), Ex)
    out['NC1_technical_only'] = {'covariates': list(tech.columns) + ['platform'], 'median_cvR2': float(np.median(t)),
                                 'frac_genes_MA_gt_technical': float(np.mean(ma > t))}
else:
    out['NC1_technical_only'] = 'technical covariate file not available'

# negative control 2: mismatched promoters within expression-mean decile
dec = pd.qcut(Ex.mean(1), 10, labels=False, duplicates='drop')
null_med = []
for d in range(100):
    perm = np.arange(len(gV))
    for q in np.unique(dec):
        ii = np.where(dec == q)[0]; perm[ii] = rng.permutation(ii)
    _, _, ma_n, _ = h2_models(Mx[perm], Ax[perm], Ex, covP)
    null_med.append(np.median(ma_n))
out['NC2_mismatched_promoter'] = {'median_cvR2_MA_null_mean': float(np.mean(null_med)),
                                  'null_95pct': [float(np.percentile(null_med, 2.5)), float(np.percentile(null_med, 97.5))],
                                  'observed_median_cvR2_MA': float(np.median(ma))}

# ---------------- H3 ----------------
def h3(sel_mask, label, nperm=1000):
    basal = (pam == 'Basal').values & sel_mask; lum = pam.isin(['LumA', 'LumB']).values & sel_mask
    E = EXP.values
    t, p = stats.ttest_ind(E[:, basal], E[:, lum], axis=1, equal_var=False)
    lfc = E[:, basal].mean(1) - E[:, lum].mean(1)
    ok_ = ~np.isnan(p); fdr = np.full(len(p), np.nan); fdr[ok_] = multipletests(p[ok_], method='fdr_bh')[1]
    res = pd.DataFrame({'gene': EXP.index, 'log2FC_basal_vs_lum': lfc, 'p': p, 'FDR': fdr})
    def classify(genes_, direction, b, l):
        rows = []
        for gname in genes_:
            md = (METH.loc[gname].values[b].mean() - METH.loc[gname].values[l].mean()) if gname in METH.index else np.nan
            ad = (ACC.loc[gname].values[b].mean() - ACC.loc[gname].values[l].mean()) if gname in ACC.index else np.nan
            mflag = (md >= 0.2) if direction == 'basal' else (md <= -0.2)
            aflag = (ad <= -1) if direction == 'basal' else (ad >= 1)
            mflag = bool(mflag) if not np.isnan(md) else False; aflag = bool(aflag) if not np.isnan(ad) else False
            rows.append((gname, md, ad, 'both' if mflag and aflag else 'methylation' if mflag else 'accessibility' if aflag else 'neither',
                         gname in METH.index, gname in ACC.index))
        return pd.DataFrame(rows, columns=['gene', 'meth_diff', 'acc_diff', 'class', 'has_meth', 'has_peak'])
    bs = res[(res.log2FC_basal_vs_lum <= -1) & (res.FDR < 0.05)].gene
    ls = res[(res.log2FC_basal_vs_lum >= 1) & (res.FDR < 0.05)].gene
    cb, cl = classify(bs, 'basal', basal, lum), classify(ls, 'lum', basal, lum)
    def props(df):
        d = {'n': int(len(df))}
        for k in ['methylation', 'accessibility', 'both', 'neither']:
            x = int((df['class'] == k).sum()); lo, hi = proportion_confint(x, max(len(df), 1), method='wilson')
            d[k] = {'n': x, 'prop': x / max(len(df), 1), 'CI95': [float(lo), float(hi)]}
        d['with_meth_and_peak'] = int((df.has_meth & df.has_peak).sum())
        return d
    stat = lambda cb_, cl_: (cb_['class'].isin(['methylation', 'both']).mean() if len(cb_) else np.nan) - (cl_['class'].isin(['methylation', 'both']).mean() if len(cl_) else np.nan)
    obs = stat(cb, cl)
    o = {'n_basal': int(basal.sum()), 'n_luminal': int(lum.sum()), 'basal_silenced': props(cb), 'luminal_silenced': props(cl),
         'diff_prop_methylation_assoc_basal_minus_luminal': float(obs)}
    if nperm:
        idx = np.where(basal | lum)[0]; nb = int(basal.sum()); null = []
        for i in range(nperm):
            pi = rng.permutation(idx); b2 = np.zeros(len(cohort), bool); l2 = np.zeros(len(cohort), bool)
            b2[pi[:nb]] = True; l2[pi[nb:]] = True
            t_, p_ = stats.ttest_ind(E[:, b2], E[:, l2], axis=1, equal_var=False); lf = E[:, b2].mean(1) - E[:, l2].mean(1)
            okk = ~np.isnan(p_); f_ = np.full(len(p_), np.nan); f_[okk] = multipletests(p_[okk], method='fdr_bh')[1]
            bsg = EXP.index[(lf <= -1) & (f_ < 0.05)]; lsg = EXP.index[(lf >= 1) & (f_ < 0.05)]
            if len(bsg) and len(lsg): null.append(stat(classify(bsg, 'basal', b2, l2), classify(lsg, 'lum', b2, l2)))
        null = np.array(null)
        o['perm_n_valid'] = int(len(null))
        o['perm_p_two_sided'] = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)) if len(null) else None
    res.to_csv(f'{R}/H3_{label}_DE.csv', index=False)
    pd.concat([cb.assign(direction='basal_silenced'), cl.assign(direction='luminal_silenced')]).to_csv(f'{R}/H3_{label}_classes.csv', index=False)
    return o

NPERM = int(os.environ.get('NPERM', '1000'))
out['H3_all'] = h3(np.ones(len(cohort), bool), 'all', NPERM)
out['H3_HM450_only'] = h3(is450 == 1, 'HM450', NPERM)

json.dump(out, open(f'{R}/results.json', 'w'), indent=2, default=float)
print(json.dumps({k: out[k] for k in ['n_tumours', 'n_by_platform', 'genes', 'H2_genes']}, indent=1, default=float))
