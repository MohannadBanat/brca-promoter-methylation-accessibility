"""Figures. Reads results/ from 04_analysis.py and 05_jackknife_h2.py; writes 600 dpi PNG and TIFF.
Run: python3 -I scripts/06_figures.py results figures"""
import json, sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

R, F = sys.argv[1], sys.argv[2]; os.makedirs(F, exist_ok=True)
res = json.load(open(f'{R}/results.json')); jk = json.load(open(f'{R}/H2_jackknife.json'))
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.spines.top': False, 'axes.spines.right': False})
C = {'M': '#4C72B0', 'A': '#DD8452', 'MA': '#55A868', 'grey': '#BDBDBD', 'HCP': '#4C72B0', 'LCP': '#C44E52'}
def save(fig, name):
    fig.savefig(f'{F}/{name}.png', dpi=600, bbox_inches='tight')
    fig.savefig(f'{F}/{name}.tiff', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(fig)
def tag(ax, t): ax.text(-0.18, 1.06, t, transform=ax.transAxes, fontweight='bold', fontsize=9)

# Figure 1: H1
h1 = pd.read_csv(f'{R}/H1_per_gene_rho.csv')
fig, ax = plt.subplots(1, 2, figsize=(6.7, 2.7), gridspec_kw={'width_ratios': [1.5, 1]})
bins = np.linspace(-1, 1, 41)
for lab, m, col in [('CpG-island promoters', h1.HCP == True, C['HCP']), ('Non-island promoters', h1.HCP != True, C['LCP'])]:
    ax[0].hist(h1.rho[m], bins=bins, alpha=.55, label=f'{lab} (n = {int(m.sum())})', color=col)
ax[0].axvline(-0.3, ls='--', c='k', lw=.7)
ax[0].set_xlabel('Spearman ρ across tumours\n(promoter methylation vs promoter accessibility)'); ax[0].set_ylabel('Genes')
ax[0].legend(frameon=False, fontsize=6.5, loc='upper right', bbox_to_anchor=(1.02, 1.0)); ax[0].set_xlim(-1, 1); ax[0].set_ylim(0, 150); tag(ax[0], 'A')
H = res['H1']; vals = [H['methylated_yet_accessible'], H['unmethylated_yet_closed']]
ax[1].bar([0, 1], [v / H['pairs'] * 100 for v in vals], color=[C['M'], C['A']], width=.6)
ax[1].set_xticks([0, 1], ['Methylated\n(β ≥ 0.6), open', 'Unmethylated\n(β ≤ 0.2), closed'], fontsize=6.5)
ax[1].set_ylabel('% of gene–tumour pairs'); tag(ax[1], 'B')
for i, v in enumerate(vals): ax[1].text(i, v / H['pairs'] * 100 + 0.2, f'{v:,}', ha='center', va='bottom', fontsize=6.5)
fig.tight_layout(); save(fig, 'Figure1')

# Figure 2: H2
h2 = pd.read_csv(f'{R}/H2_per_gene_cvR2.csv')
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.8), gridspec_kw={'width_ratios': [1.1, 1, 1.25]})
bp = ax[0].boxplot([h2.cvR2_M, h2.cvR2_A, h2.cvR2_MA], showfliers=False, patch_artist=True, widths=.6, medianprops={'color': 'k'})
for p, c in zip(bp['boxes'], [C['M'], C['A'], C['MA']]): p.set_facecolor(c); p.set_alpha(.75)
ax[0].set_xticks([1, 2, 3], ['Methyl-\nation', 'Accessi-\nbility', 'Both'])
nc2 = res['NC2_mismatched_promoter']; nc1 = res['NC1_technical_only']
ax[0].axhline(nc2['median_cvR2_MA_null_mean'], ls=':', c='k', lw=.8, label='Mismatched-promoter null')
ax[0].axhline(nc1['median_cvR2'], ls='--', c='#7F7F7F', lw=.8, label='Technical-only model')
ax[0].axhline(0, c='k', lw=.4); ax[0].set_ylabel('Cross-validated R² per gene'); ax[0].set_ylim(-0.45, 1.3)
ax[0].legend(frameon=False, fontsize=5.8, loc='upper left', title='Medians of controls', title_fontsize=5.8); tag(ax[0], 'A')
lim = (-0.5, 1.0)
xm, ya = h2.cvR2_M.clip(*lim), h2.cvR2_A.clip(*lim)
ax[1].scatter(xm, ya, s=2.5, alpha=.35, c=np.where(h2.HCP, C['HCP'], C['LCP']), rasterized=True, lw=0)
ax[1].plot(lim, lim, c='k', lw=.5); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel('CV-R², methylation only'); ax[1].set_ylabel('CV-R², accessibility only')
ax[1].text(0.97, 0.05, f"Accessibility > methylation\nin {100 * np.mean(h2.cvR2_A > h2.cvR2_M):.0f}% of {len(h2):,} genes", transform=ax[1].transAxes, ha='right', fontsize=6)
ax[1].scatter([], [], c=C['HCP'], s=6, label='CpG island'); ax[1].scatter([], [], c=C['LCP'], s=6, label='Non-island')
ax[1].legend(frameon=False, fontsize=5.8, loc='upper left', handletextpad=.2); tag(ax[1], 'B')
labs = ['Pooled', 'Purity-\nadjusted', 'HM450\nonly', 'HM27\nonly']; S = [jk['pooled'], jk['purity_adjusted'], jk['HM450'], jk['HM27']]
x = np.arange(4); wdt = .38
for j, (k, col, nm) in enumerate([('median_MA_minus_M', C['A'], 'Adding accessibility\nto methylation'), ('median_MA_minus_A', C['M'], 'Adding methylation\nto accessibility')]):
    v = np.array([s[k]['estimate'] for s in S]); lo = np.array([s[k]['CI95'][0] for s in S]); hi = np.array([s[k]['CI95'][1] for s in S])
    ax[2].bar(x + (j - .5) * wdt, v, wdt, color=col, alpha=.85, label=nm, yerr=[v - lo, hi - v], capsize=2, error_kw={'lw': .6})
ax[2].axhline(0, c='k', lw=.4); ax[2].set_xticks(x, labs, fontsize=6.3); ax[2].set_ylabel('Median gain in CV-R² (95% CI)')
ax[2].set_ylim(-0.06, 0.40); ax[2].legend(frameon=False, fontsize=5.8, loc='upper right'); tag(ax[2], 'C')
fig.tight_layout(w_pad=1.2); save(fig, 'Figure2')

# Figure 3: H3 among silenced genes with both a methylation value and a promoter peak
t = json.load(open(f'{R}/H3_both_measures_table.json'))
fig, ax = plt.subplots(1, 2, figsize=(6.0, 2.7), sharey=True)
cats = ['methylation', 'both', 'accessibility', 'neither']; cols = [C['M'], C['MA'], C['A'], C['grey']]
names = {'methylation': 'Methylation only', 'both': 'Both', 'accessibility': 'Accessibility only', 'neither': 'Neither'}
for a, key, ttl, hk in [(ax[0], 'all', 'All tumours', 'H3_all'), (ax[1], 'HM450', 'HM450 tumours only', 'H3_HM450_only')]:
    for i, grp in enumerate(['basal_silenced', 'luminal_silenced']):
        g = t[key][grp]; bottom = 0
        for c, col in zip(cats, cols):
            p = g[c][1] * 100; a.bar(i, p, bottom=bottom, color=col, width=.6, label=names[c] if i == 0 else None); bottom += p
        a.text(i, 101.5, f"n = {g['n']}", ha='center', fontsize=6.3)
    a.set_xticks([0, 1], ['Lower in Basal', 'Lower in Luminal']); a.set_ylim(0, 108)
    a.set_title(f"{ttl}\n(Basal {res[hk]['n_basal']}, Luminal {res[hk]['n_luminal']})", fontsize=7)
ax[0].set_ylabel('% of silenced genes'); tag(ax[0], 'A'); ax[1].text(-0.08, 1.06, 'B', transform=ax[1].transAxes, fontweight='bold', fontsize=9)
ax[1].legend(frameon=False, fontsize=6, loc='center left', bbox_to_anchor=(1.0, 0.5))
fig.tight_layout(); save(fig, 'Figure3')
print('figures written')
