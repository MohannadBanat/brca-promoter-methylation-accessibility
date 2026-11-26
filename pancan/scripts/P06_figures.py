"""Figures for the pan-cancer manuscript (Figure 1 uses posthoc_review_corrections.json from P08), sized for the MDPI text column (13.8 cm = 5.43 in). 600 dpi PNG + TIFF.
Run: python3 -I P06_figures.py <results dir> <figure dir>"""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R, F = sys.argv[1:3]; os.makedirs(F, exist_ok=True)
p = json.load(open(f'{R}/pancan_results.json')); b = json.load(open(f'{R}/breast_rerun.json'))['masked_M_general']
s5 = json.load(open(f'{R}/S5_meta_with_breast.json')); ph = json.load(open(f'{R}/posthoc_robustness.json'))
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'legend.fontsize': 6.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False})
C = {'M': '#4C72B0', 'A': '#DD8452', 'MA': '#55A868', 'grey': '#9E9E9E'}
NAME = {'KIRP': 'Kidney papillary', 'PRAD': 'Prostate', 'COAD': 'Colon', 'LUAD': 'Lung adeno.', 'STAD': 'Stomach', 'ESCA': 'Esophagus',
        'LIHC': 'Liver', 'LUSC': 'Lung squamous', 'KIRC': 'Kidney clear cell', 'BRCA': 'Breast'}
T9 = p['meta_types_included']
def save(fig, name):
    fig.savefig(f'{F}/{name}.png', dpi=600, bbox_inches='tight')
    fig.savefig(f'{F}/{name}.tiff', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(fig)
def tag(ax, t, x=-0.04): ax.text(x, 1.04, t, transform=ax.transAxes, fontweight='bold', fontsize=9, ha='right')

# ---------- Figure 1: forest plots of D and G ----------
rows = [(f"{NAME[T]} ({T}, n = {p['per_type'][T]['n_tumours']})", p['per_type'][T]['H2'], 'type') for T in T9]
fig, axs = plt.subplots(1, 2, figsize=(5.43, 3.5), sharey=True, gridspec_kw={'wspace': 0.08})
for ax, key, pooled, lab, title in [(axs[0], 'D_median_A_minus_M', p['R1_pooled_D'], 'D = median(CV-R²$_A$ − CV-R²$_M$)', 'Accessibility vs methylation'),
                                    (axs[1], 'G_median_MA_minus_A', p['R3_pooled_G'], 'G = median(CV-R²$_{MA}$ − CV-R²$_A$)', 'Methylation added to accessibility')]:
    y = list(range(len(rows) + 4, 4, -1))
    for yy, (nm, h, _) in zip(y, rows):
        e = h[key]; ax.errorbar(e['estimate'], yy, xerr=[[e['estimate'] - e['CI95'][0]], [e['CI95'][1] - e['estimate']]], fmt='s', ms=3.5, c='k', lw=.8, capsize=0)
    e = pooled; yy = 3.6
    ax.fill([e['CI95'][0], e['estimate'], e['CI95'][1], e['estimate']], [yy, yy + .32, yy, yy - .32], c='k')
    e = s5[key]; yy = 2.6
    ax.fill([e['CI95'][0], e['estimate'], e['CI95'][1], e['estimate']], [yy, yy + .32, yy, yy - .32], c=C['grey'])
    e = b[key]; yy = 1.4
    ax.errorbar(e['estimate'], yy, xerr=[[e['estimate'] - e['CI95'][0]], [e['CI95'][1] - e['estimate']]], fmt='o', ms=3.5, mfc='white', c='k', lw=.8)
    e = p['S2_pooled_within_type'][key]; yy = 0.4
    ax.errorbar(e['estimate'], yy, xerr=[[e['estimate'] - e['CI95'][0]], [e['CI95'][1] - e['estimate']]], fmt='D', ms=3, c=C['A'], lw=.8)
    ax.axvline(0, c='k', lw=.6)
    if key.startswith('G'): ax.axvline(0.02, c='k', lw=.6, ls='--'); ax.text(0.022, len(rows) + 4.6, 'margin\n0.02', fontsize=6, va='top')
    ax.set_xlabel(lab); ax.set_title(title, fontsize=7.5)
    ax.axhline(4.3, c=C['grey'], lw=.4); ax.axhline(2.0, c=C['grey'], lw=.4, ls=':')
ylab = [r[0] for r in rows] + ['Pooled, 9 test types (k = 9)', 'Pooled with breast (k = 10)', 'Breast, discovery (n = 73)', f"All types, within-type (n = {p['S2_pooled_within_type']['n_tumours']})"]
axs[0].set_yticks(list(range(len(rows) + 4, 4, -1)) + [3.6, 2.6, 1.4, 0.4], ylab, fontsize=6.5); axs[0].set_ylim(-0.3, len(rows) + 4.8)
tag(axs[0], 'A', -0.02); tag(axs[1], 'B', -0.02)
save(fig, 'Figure2')

# ---------- Figure 1: median CV-R2 of M, A, MA by type, relative to the type mean (post hoc re-basing, deviations log entry 5) ----------
rc = json.load(open(f'{R}/posthoc_review_corrections.json'))
ren = lambda g: {'median_cvR2_M': g['median_rel_M'], 'median_cvR2_A': g['median_rel_A'], 'median_cvR2_MA': g['median_rel_MA']}
cols = [(T, ren(rc['per_type'][T]['g_relative'])) for T in T9] + [('BRCA', ren(rc['BRCA']['g_relative']))]
s2 = p['S2_pooled_within_type']
fig, ax = plt.subplots(figsize=(5.43, 2.6))
x = np.arange(len(cols) + 1)
for j, (k, lab, c) in enumerate([('median_cvR2_M', 'Methylation (M)', C['M']), ('median_cvR2_A', 'Accessibility (A)', C['A']), ('median_cvR2_MA', 'Both (MA)', C['MA'])]):
    est = [h[k]['estimate'] for _, h in cols] + [s2[k]['estimate']]
    lo = [h[k]['CI95'][0] for _, h in cols] + [s2[k]['CI95'][0]]; hi = [h[k]['CI95'][1] for _, h in cols] + [s2[k]['CI95'][1]]
    ax.errorbar(x + (j - 1) * 0.24, est, yerr=[np.subtract(est, lo), np.subtract(hi, est)], fmt='o', ms=3.2, c=c, lw=.8, label=lab)
ax.axhline(0, c='k', lw=.6); ax.axvline(len(T9) - 0.5, c=C['grey'], lw=.5, ls=':'); ax.axvline(len(T9) + 0.5, c=C['grey'], lw=.5, ls=':')
ax.set_xticks(x, [T for T, _ in cols] + ['All types\n(within-type)'], fontsize=6.5)
ax.set_ylabel('Median per-gene CV-R²\nrelative to type mean (95% CI)'); ax.legend(frameon=False, ncol=3, loc='upper left', bbox_to_anchor=(0, 1.12))
save(fig, 'Figure1')

# ---------- Figure 3: genes where methylation is informative ----------
rec = pd.read_csv(f'{R}/posthoc_methylation_informative_recurrence.csv')
pg = {T: pd.read_csv(f'{R}/per_gene_{T}.csv') for T in T9}
pg['BRCA'] = pd.read_csv(sys.argv[3]) if len(sys.argv) > 3 else None
coh = T9 + ['BRCA']
inf = {T: set(pg[T].gene[pg[T].cvR2_M > 0.3].astype(str)) for T in coh}; tested = {T: pg[T].gene.astype(str).values for T in coh}
cnt = pd.Series([g for T in coh for g in inf[T]]).value_counts()
rng = np.random.default_rng(20261126); ks = np.arange(1, 10); null = []
for _ in range(1000):
    c = pd.Series(np.concatenate([rng.choice(tested[T], len(inf[T]), replace=False) for T in coh])).value_counts()
    null.append([(c >= k).sum() for k in ks])
null = np.array(null); obs = [(cnt >= k).sum() for k in ks]
fig = plt.figure(figsize=(5.43, 3.9)); gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.35], wspace=0.45)
ax = fig.add_subplot(gs[0, 0])
ax.plot(ks, obs, 'o-', c=C['M'], ms=3.5, lw=1, label='Observed')
ax.plot(ks, null.mean(0), 's--', c=C['grey'], ms=3, lw=1, label='Random genes (mean)')
ax.fill_between(ks, np.percentile(null, 2.5, axis=0), np.percentile(null, 97.5, axis=0), color=C['grey'], alpha=.25, lw=0)
ax.set_yscale('symlog', linthresh=1); ax.set_xlabel('Methylation informative in ≥ k\nof 10 cohorts (CV-R²$_M$ > 0.3)'); ax.set_ylabel('Genes')
ax.set_xticks(ks)
ax.text(6.2, 45, 'Observed', color=C['M'], fontsize=6.5)
ax.text(0.85, 0.25, 'Random genes\n(mean and\n95% range)', color='#616161', fontsize=6.5); tag(ax, 'A', -0.12)
ax = fig.add_subplot(gs[0, 1]); top = rec[rec['types_M_gt_0.3'] >= 6].gene.astype(str).tolist()
for i, g in enumerate(top):
    for j, T in enumerate(coh):
        if g in inf[T]: ax.plot(j, i, 's', c=C['M'], ms=4.2)
        elif g in set(tested[T]): ax.plot(j, i, 's', mfc='white', mec=C['grey'], ms=4.2, mew=.6)
ax.set_yticks(range(len(top)), [f'$\\it{{{g}}}$' for g in top], fontsize=6); ax.invert_yaxis()
ax.set_xticks(range(len(coh)), coh, rotation=90, fontsize=6.5); ax.set_xlim(-.6, len(coh) - .4)
ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False); ax.tick_params(length=0)
ax.plot([], [], 's', c=C['M'], ms=4, label='CV-R²$_M$ > 0.3'); ax.plot([], [], 's', mfc='white', mec=C['grey'], ms=4, label='Tested, ≤ 0.3')
ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2); tag(ax, 'B', -0.2)
save(fig, 'Figure3')
print('obs', obs, 'null mean', null.mean(0).round(1).tolist(), 'top genes', len(top))
