"""Supplementary Tables S1-S6 (one workbook) for the pan-cancer manuscript.
Run: python3 -I build_supplementary.py <pancan results> <pancan plan dir> <breast ijms dir> <out.xlsx>"""
import json, re, sys
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side
from openpyxl.utils import get_column_letter
R, PLAN, BIJ, OUT = sys.argv[1:5]
p = json.load(open(f'{R}/pancan_results.json')); b = json.load(open(f'{R}/breast_rerun.json'))
s5 = json.load(open(f'{R}/S5_meta_with_breast.json')); ph = json.load(open(f'{R}/posthoc_robustness.json'))
nl = json.load(open(f'{R}/posthoc_single_predictor_null.json'))['per_type']
rc = json.load(open(f'{R}/posthoc_review_corrections.json'))
T9 = p['meta_types_included']; pt = p['per_type']
TITLE = ('Promoter DNA Methylation Adds Little to Promoter Chromatin Accessibility in Explaining Gene Expression Differences '
         'Between Tumors: A Pre-Specified Analysis of Ten Cancer Types')
AUTH = 'Hamza Abu Owida, Ayoub Alsarhan, Mohannad Banat, Dimah Al-Fraihat and Abdulmohsen Algarni'
F = 'Arial'; fT = Font(name=F, size=12, bold=True); fB = Font(name=F, size=10, bold=True); f = Font(name=F, size=10); fN = Font(name=F, size=9, italic=True)
thin = Side(style='thin', color='000000'); wrap = Alignment(wrap_text=True, vertical='top')
wb = Workbook()

def sheet(name, caption, note, headers, rows, widths, numfmt='0.000', wrapcols=()):
    ws = wb.create_sheet(name); ws['A1'] = caption; ws['A1'].font = fT
    if note: ws['A2'] = note; ws['A2'].font = fN
    for c, (h, w) in enumerate(zip(headers, widths), 1):
        x = ws.cell(row=4, column=c, value=h); x.font = fB; x.alignment = Alignment(wrap_text=True, vertical='center')
        x.border = Border(top=thin, bottom=thin); ws.column_dimensions[get_column_letter(c)].width = w
    for i, r in enumerate(rows):
        for c, v in enumerate(r, 1):
            x = ws.cell(row=5 + i, column=c, value=v); x.font = f
            if isinstance(v, float): x.number_format = numfmt
            if c in wrapcols: x.alignment = wrap
    last = 4 + len(rows)
    for c in range(1, len(headers) + 1): ws.cell(row=last, column=c).border = Border(bottom=thin)
    ws.freeze_panes = 'A5'; ws.sheet_view.showGridLines = False
    if len(rows) > 30: ws.auto_filter.ref = f'A4:{get_column_letter(len(headers))}{last}'
    return ws

# ---------- Contents ----------
ws = wb.active; ws.title = 'Contents'; ws.sheet_view.showGridLines = False
ws['A1'] = 'Supplementary Materials'; ws['A1'].font = Font(name=F, size=14, bold=True)
ws['A3'] = TITLE; ws['A3'].font = fB; ws['A3'].alignment = wrap; ws.merge_cells('A3:C3'); ws.row_dimensions[3].height = 30
ws['A4'] = AUTH; ws['A4'].font = f
contents = [('Table S1', 'Deviations from the two committed analysis plans and their timing, with the full SHA-256 hashes.'),
            ('Table S2', 'Per-gene leave-one-out cross-validated R² in every cohort (9 test types and breast).'),
            ('Table S3', 'Breast discovery cohort: detailed results, including the Basal-versus-luminal comparison (exploratory in this paper).'),
            ('Table S4', 'Post hoc analyses (deviations log entries 3-5): gene selection within jackknife, probe-number strata, two-state methylation model, ATAC-seq quality, single-predictor null, CV-R² relative to the type mean, pooled coupling, Hartung-Knapp-Sidik-Jonkman and prediction intervals.'),
            ('Table S5', 'Full per-type results: gene counts, variance analysis, coupling, negative controls, secondary analyses S1-S7; exploratory results for the 8 non-breast types with 5-14 tumors.'),
            ('Table S6', 'Methylation-informative genes (CV-R² of methylation alone > 0.3): counts per cohort, recurrence across cohorts, comparison with random genes and with accessibility-informative genes.')]
for i, (a, d) in enumerate([('Sheet', 'Description')] + contents):
    for c, v in enumerate((a, d), 1):
        x = ws.cell(row=6 + i, column=c, value=v); x.font = fB if i == 0 else f; x.alignment = wrap
ws.column_dimensions['A'].width = 12; ws.column_dimensions['B'].width = 120
ws['A14'] = ('Code, both hashed analysis plans, deviations logs and all outputs: https://github.com/MohannadBanat/brca-promoter-methylation-accessibility '
             '(folder pancan/ for this paper).'); ws['A14'].font = fN

# ---------- S1 deviations ----------
rows = []
md = open(f'{BIJ}/manuscript/TableS1_deviations.md', encoding='utf-8').read()
for num, head, text in re.findall(r'^(\d+)\. \*\*(.+?)\*\*\s*(.*)$', md, flags=re.M):
    rows.append(['Breast plan (16 Nov 2026; SHA-256 b2b7079a…8d2d)', int(num), head.rstrip('.'), 'Before results viewed' if int(num) <= 6 else 'After results viewed',
                 text.replace('R2', 'R²').replace('+/-', '±')])
lg = open(f'{PLAN}/pancan_deviations_log.md', encoding='utf-8').read()
for num, head, text in re.findall(r'^(\d+)\. \*\*(.+?)\*\*\s*(.*)$', lg, flags=re.M):
    timing = 'Before any pan-cancer result was viewed' if int(num) == 1 else ('After primary results; before the analyses it adds were run' if int(num) in (3, 4, 5) else 'After primary results viewed')
    rows.append(['Pan-cancer plan (25 Nov 2026; SHA-256 c440f092…fa3b)', int(num), head.rstrip('.'), timing, text.replace('`', '')])
sheet('Table S1', 'Table S1. Deviations from the two committed analysis plans and their timing.',
      'Breast plan SHA-256: b2b7079a364df1e80b5e9c5c0c46947204045ea7e401b5e1b1391ce3e8b48d2d. Pan-cancer plan SHA-256: c440f092afb11f7e4afb9a7204ef51550ff12d0d88ad6f5c765527fbd542fa3b (committed 2026-11-25T00:03:36Z). Entries are copied verbatim from the deviations logs archived with the code; each states when it was written relative to the results it concerns.',
      ['Plan', 'No.', 'Topic', 'Written', 'Description'], rows, [34, 5, 40, 26, 120], wrapcols=(1, 3, 4, 5))

# ---------- S2 per type ----------
def E(d, k): return None if not d or k not in d else float(d[k]['estimate'])
def L(d, k): return None if not d or k not in d else float(d[k]['CI95'][0])
def H(d, k): return None if not d or k not in d else float(d[k]['CI95'][1])
types = [T for T in pt if pt[T]['n_tumours'] >= 5 and pt[T].get('H2')]
types = T9 + sorted(t for t in types if t not in T9)
hdr = ['Cancer type', 'Role', 'Tumors', 'Mean purity', 'Genes with promoter peak', 'Genes with ≥2 probes', 'Expressed genes', 'Complete cases', 'Variable genes', 'Variable, CpG-island']
for k, lab in [('median_cvR2_M', 'M'), ('median_cvR2_A', 'A'), ('median_cvR2_MA', 'MA'), ('D_median_A_minus_M', 'D'), ('G_median_MA_minus_A', 'G'), ('frac_A_gt_M', 'Fraction A > M')]:
    hdr += [f'{lab}', f'{lab} 95% CI low', f'{lab} 95% CI high']
hdr += ['S1 genes', 'S1 median ρ', 'S1 fraction ρ ≤ −0.3', 'NC1 technical-only median CV-R²', 'NC1 fraction genes MA > technical',
        'NC2 null median CV-R²(MA), mean', 'NC2 null 2.5%', 'NC2 null 97.5%', 'S3 CpG-island D', 'S3 CpG-island G', 'S3 non-island D', 'S3 non-island G',
        'S4 purity tumors', 'S4 purity-adjusted D', 'S4 purity-adjusted G', 'S7 breast-probe-set genes', 'S7 D', 'S7 G']
rows = []
for T in types:
    v = pt[T]; h = v['H2']; row = [T, 'Test (pre-specified)' if T in T9 else 'Exploratory (5-14 tumors)', v['n_tumours'], v.get('mean_purity'),
                                   v['genes']['with_promoter_peak'], v['genes']['with_>=2_probes'], v['genes']['expressed'], v['genes']['complete_cases'], v['genes']['variable'], v['genes']['variable_HCP']]
    for k in ['median_cvR2_M', 'median_cvR2_A', 'median_cvR2_MA', 'D_median_A_minus_M', 'G_median_MA_minus_A', 'frac_A_gt_M']: row += [E(h, k), L(h, k), H(h, k)]
    nc1, nc2 = v.get('NC1_technical_only') or {}, v.get('NC2_mismatched') or {}
    row += [v['S1']['n_genes'], v['S1']['median_rho'], v['S1']['frac_rho_le_-0.3'], nc1.get('median_cvR2'), nc1.get('frac_genes_MA_gt_technical'),
            nc2.get('null_mean_median_cvR2_MA'), (nc2.get('null_95pct') or [None, None])[0], (nc2.get('null_95pct') or [None, None])[1],
            E(v.get('S3_HCP'), 'D_median_A_minus_M'), E(v.get('S3_HCP'), 'G_median_MA_minus_A'), E(v.get('S3_LCP'), 'D_median_A_minus_M'), E(v.get('S3_LCP'), 'G_median_MA_minus_A'),
            (v.get('S4_purity') or {}).get('n_tumours'), E(v.get('S4_purity'), 'D_median_A_minus_M'), E(v.get('S4_purity'), 'G_median_MA_minus_A'),
            (v.get('S7_shared_probes') or {}).get('n_genes'), E(v.get('S7_shared_probes'), 'D_median_A_minus_M'), E(v.get('S7_shared_probes'), 'G_median_MA_minus_A')]
    rows.append(row)
ws = sheet('Table S5', 'Table S5. Full per-type results.',
           'M, A, MA: median per-gene leave-one-out CV-R² of models on promoter methylation, accessibility or both; D = median(A − M), G = median(MA − A); 95% delete-one-tumor jackknife CIs. '
           'Exploratory types (5-14 tumors) are not part of any pre-specified test and their small samples make CV-R² strongly negative; they contribute to the pooled within-type analysis (S2). '
           'Types with fewer than 5 tumors (CESC, CHOL, UCEC, SKCM) are omitted (Table S1, pan-cancer entry 2). CV-R² on the scale of the plan (1 − PRESS/SST); see Table S4 for values relative to the type mean.',
           hdr, rows, [10, 22] + [11] * (len(hdr) - 2))
r0 = 6 + len(rows)
ws.cell(row=r0, column=1, value='Heterogeneity (S6), pooled coupling (S1) and the pooled within-type analysis (S2)').font = fB
het = [('Spearman ρ, D vs tumor count', float(pd.Series([pt[T]['n_tumours'] for T in T9]).corr(pd.Series([pt[T]['H2']['D_median_A_minus_M']['estimate'] for T in T9]), method='spearman'))),
       ('Spearman ρ, D vs mean purity', float(pd.Series([pt[T]['mean_purity'] for T in T9]).corr(pd.Series([pt[T]['H2']['D_median_A_minus_M']['estimate'] for T in T9]), method='spearman'))),
       ('I² of pooled D (R1)', p['R1_pooled_D']['I2']), ('τ² of pooled D', p['R1_pooled_D']['tau2']), ('I² of pooled G (R3)', p['R3_pooled_G']['I2']), ('τ² of pooled G', p['R3_pooled_G']['tau2'])]
s2p = p['S2_pooled_within_type']; sp = rc['h_S1_pooled_DL']
het += [('S1: pooled median ρ, DerSimonian–Laird (entry 5)', sp['estimate']), ('S1: pooled median ρ, 95% CI low', sp['CI95'][0]), ('S1: pooled median ρ, 95% CI high', sp['CI95'][1]),
        ('S1: tumor-weighted mean of median ρ (as computed by the analysis code)', p['S1_weighted_mean_of_median_rho'])]
for k, lab in [('median_cvR2_M', 'M'), ('median_cvR2_A', 'A'), ('median_cvR2_MA', 'MA'), ('D_median_A_minus_M', 'D'), ('G_median_MA_minus_A', 'G'), ('median_MA_minus_M', 'MA − M'), ('frac_A_gt_M', 'fraction A > M')]:
    het += [(f"S2, all types pooled within type ({s2p['n_tumours']} tumors, {s2p['n_genes']} genes): {lab}", s2p[k]['estimate']), (f'S2 {lab}: 95% CI low', s2p[k]['CI95'][0]), (f'S2 {lab}: 95% CI high', s2p[k]['CI95'][1])]
for i, (a, v) in enumerate(het):
    ws.cell(row=r0 + 1 + i, column=1, value=a).font = f; x = ws.cell(row=r0 + 1 + i, column=3, value=float(v)); x.font = f; x.number_format = '0.000'

# ---------- S3 per gene ----------
rows = []
for T in T9:
    d = pd.read_csv(f'{R}/per_gene_{T}.csv')
    rows += [[T, r.gene, 'Yes' if r.HCP else 'No', float(r.cvR2_M), float(r.cvR2_A), float(r.cvR2_MA)] for r in d.itertuples()]
d = pd.read_csv(f'{BIJ}/manuscript/TableS2_per_gene_cvR2.csv')
rows += [['BRCA', r.gene, 'Yes' if bool(r.HCP) else 'No', float(r.cvR2_M), float(r.cvR2_A), float(r.cvR2_MA)] for r in d.itertuples()]
sheet('Table S2', 'Table S2. Per-gene leave-one-out cross-validated R² (CV-R²) in every cohort.',
      'Genes variable in all three layers in that cohort. BRCA models include array platform. CV-R² = 1 − PRESS/SST as in the plan, where SST is about the mean of all tumors; predicting each held-out tumor by the mean of the others scores 1 − (n/(n − 1))² on this scale.',
      ['Cohort', 'Gene', 'CpG-island promoter', 'CV-R² (M)', 'CV-R² (A)', 'CV-R² (MA)'], rows, [9, 14, 12, 11, 11, 11], numfmt='0.0000')

# ---------- S4 post hoc ----------
hdr = ['Cancer type', 'Tumors', '(a) D, selection within jackknife', '(a) 95% CI low', '(a) 95% CI high', '(a) G', '(a) G 95% CI low', '(a) G 95% CI high', '(a) genes per replicate, range']
for s in ['2-3', '4-7', '>=8']: hdr += [f'(b) {s} probes: genes', f'(b) {s}: median CV-R² (M)', f'(b) {s}: D', f'(b) {s}: G']
hdr += ['(d) two-state: median CV-R²', '(d) linear M: median CV-R²', '(d) fraction of genes two-state > linear', '(d) median gain, two-state added to A',
        '(f) observed median CV-R² (M)', '(f) unrelated-gene null, M (mean)', '(f) null M 2.5%', '(f) null M 97.5%',
        '(f) observed median CV-R² (A)', '(f) unrelated-gene null, A (mean)', '(f) null A 2.5%', '(f) null A 97.5%', '(e) mean TSS enrichment', '(e) mean FRiP',
        '(g) baseline on plan scale, 1 − (n/(n − 1))²']
for k, lab in [('median_rel_M', 'M'), ('median_rel_A', 'A'), ('median_rel_MA', 'MA'), ('D_rel', 'D'), ('G_rel', 'G')]:
    hdr += [f'(g) {lab} relative to type mean', f'(g) {lab} 95% CI low', f'(g) {lab} 95% CI high']
hdr += ['(h) median ρ, complete cases', '(h) jackknife SE', '(k) fraction of genes with CV-R²(M) − CV-R²(A) > 0.1']
rows = []
for T in T9:
    q = ph['per_type'][T]; a = q['a_reselect_jackknife']
    row = [T, q['n_tumours'], E(a, 'D_median_A_minus_M'), L(a, 'D_median_A_minus_M'), H(a, 'D_median_A_minus_M'), E(a, 'G_median_MA_minus_A'), L(a, 'G_median_MA_minus_A'),
           H(a, 'G_median_MA_minus_A'), f"{a['n_genes_range'][0]}–{a['n_genes_range'][1]}"]
    for s in ['2-3', '4-7', '>=8']:
        st = q['b_probe_strata'].get(s, {}); row += [st.get('n_genes'), st.get('median_cvR2_M'), st.get('D_median_A_minus_M'), st.get('G_median_MA_minus_A')]
    t = q['d_two_state']; n = nl[T]
    row += [t['median_cvR2_two_state'], t['median_cvR2_linear_M'], t['frac_two_state_gt_linear'], t['median_gain_two_state_added_to_A'],
            n['observed_median_M'], n['null_median_M_mean'], n['null_median_M_95'][0], n['null_median_M_95'][1],
            n['observed_median_A'], n['null_median_A_mean'], n['null_median_A_95'][0], n['null_median_A_95'][1],
            q['e_atac_quality']['mean_TSS_enrichment'], q['e_atac_quality']['mean_FRiP']]
    g = rc['per_type'][T]; row += [g['loo_mean_baseline_on_plan_scale']]
    for k in ['median_rel_M', 'median_rel_A', 'median_rel_MA', 'D_rel', 'G_rel']: row += [E(g['g_relative'], k), L(g['g_relative'], k), H(g['g_relative'], k)]
    dpg = pd.read_csv(f'{R}/per_gene_{T}.csv')
    row += [g['h_S1']['median_rho'], g['h_S1']['SE'], float(((dpg.cvR2_M - dpg.cvR2_A) > 0.1).mean())]
    rows.append(row)
pa = ph['a_pooled']
rows.append(['Pooled (random effects)', 191, E(pa, 'D_median_A_minus_M'), L(pa, 'D_median_A_minus_M'), H(pa, 'D_median_A_minus_M'), E(pa, 'G_median_MA_minus_A'),
             L(pa, 'G_median_MA_minus_A'), H(pa, 'G_median_MA_minus_A')] + [None] * (len(hdr) - 8))
gb = rc['BRCA']['g_relative']; brow = ['BRCA (discovery)', 73] + [None] * (len(hdr) - 2); i0 = hdr.index('(g) M relative to type mean')
for j, k in enumerate(['median_rel_M', 'median_rel_A', 'median_rel_MA', 'D_rel', 'G_rel']): brow[i0 + 3 * j:i0 + 3 * j + 3] = [E(gb, k), L(gb, k), H(gb, k)]
rows.append(brow)
gp = rc['g_pooled']; prow = ['Pooled, relative scale (random effects)', 191] + [None] * (len(hdr) - 2)
prow[hdr.index('(g) D relative to type mean'):hdr.index('(g) D relative to type mean') + 3] = [E(gp, 'D_rel'), L(gp, 'D_rel'), H(gp, 'D_rel')]
prow[hdr.index('(g) G relative to type mean'):hdr.index('(g) G relative to type mean') + 3] = [E(gp, 'G_rel'), L(gp, 'G_rel'), H(gp, 'G_rel')]
rows.append(prow)
ws = sheet('Table S4', 'Table S4. Post hoc analyses (exploratory; added after the primary results were seen; Table S1, pan-cancer entries 3-5). Letters follow the deviations log.',
           '(a) Gene selection repeated within each delete-one-tumor jackknife replicate. (b) Strata by the median number of non-missing HM450 probes in the promoter window. '
           '(d) Two-state model: tumors split per gene at the threshold maximizing between-group variance of beta (Otsu), ≥3 tumors per group. (e) ATAC-seq quality from Corces et al. Data S1. '
           '(f) Each gene given the methylation or accessibility of a random gene from the same expression-mean decile (100 draws). (g) CV-R² relative to the leave-one-out type-mean model (breast: platform-only model), 1 − PRESS/PRESS₀. '
           '(h) Median Spearman ρ among genes with complete data and its delete-one-tumor jackknife SE, for the pre-specified random-effects pooling of S1 (Table S5). (k) Exploratory count. (c, j) Recurrence: Table S6. '
           f"Spearman ρ across the nine types: TSS enrichment vs D {ph['e_spearman_TSS_enrichment_vs_D']:.2f}; TSS enrichment vs median CV-R² (A) {ph['e_spearman_TSS_enrichment_vs_median_cvR2_A']:.2f}.",
           hdr, rows, [12, 8] + [12] * (len(hdr) - 2))
r0 = 6 + len(rows); ws.cell(row=r0, column=1, value='(i) Sensitivity of the pre-specified meta-analyses to the method (k = 9)').font = fB
for i, (lab, d) in enumerate([('R1, pooled D', rc['i_R1_D']), ('R3, pooled G', rc['i_R3_G'])]):
    vals = [lab, None, d['estimate'], d['HKSJ_CI95'][0], d['HKSJ_CI95'][1], d['prediction_interval95'][0], d['prediction_interval95'][1]]
    for j, v in enumerate(vals, 1):
        x = ws.cell(row=r0 + 2 + i, column=j, value=v); x.font = f
        if isinstance(v, float): x.number_format = '0.000'
for j, h in enumerate(['Criterion', '', 'Estimate', 'HKSJ 95% CI low', 'HKSJ 95% CI high', '95% prediction interval low', '95% prediction interval high'], 1):
    ws.cell(row=r0 + 1, column=j, value=h).font = fB

# ---------- S5 methylation-informative genes ----------
rec = pd.read_csv(f'{R}/posthoc_methylation_informative_recurrence.csv')
coh = T9 + ['BRCA']; pg = {T: pd.read_csv(f'{R}/per_gene_{T}.csv') for T in T9}
pg['BRCA'] = pd.read_csv(f'{BIJ}/manuscript/TableS2_per_gene_cvR2.csv')
info = {T: dict(zip(pg[T].gene.astype(str), pg[T].cvR2_M)) for T in coh}
rr = rec[rec['types_M_gt_0.3'] >= 3]
rows = []
for r in rr.itertuples():
    rows.append([r.gene, r.chrom, int(r._2), int(r.types_tested)] + [(float(info[T][r.gene]) if r.gene in info[T] else None) for T in coh])
ws = sheet('Table S6', 'Table S6. Genes whose promoter methylation is informative (methylation-only CV-R² > 0.3) in at least three of the ten cohorts (exploratory).',
           'Cohort columns give the methylation-only CV-R² where the gene was tested (blank: not variable in all three layers in that cohort). '
           'Summary by cohort, the comparison with random genes and the comparator for accessibility-informative genes are given below the gene list. The random-gene null treats genes as interchangeable and ignores differences between genes.',
           ['Gene', 'Chromosome', 'Cohorts with CV-R²(M) > 0.3', 'Cohorts tested'] + [f'CV-R²(M), {T}' for T in coh], rows, [12, 10, 12, 10] + [11] * len(coh))
r0 = 6 + len(rows); c = ph['c_methylation_informative']
for j, h in enumerate(['Cohort', 'Genes tested', 'Genes with CV-R²(M) > 0.3', 'Fraction', 'Genes where M beats A by > 0.1'], 1):
    x = ws.cell(row=r0, column=j, value=h); x.font = fB
for i, T in enumerate(coh):
    for j, v in enumerate([T, c[T]['n_genes'], c[T]['n_M_gt_0.3'], float(c[T]['frac_M_gt_0.3']), c[T]['n_M_beats_A_by_0.1']], 1):
        x = ws.cell(row=r0 + 1 + i, column=j, value=v); x.font = f
        if isinstance(v, float): x.number_format = '0.000'
r1 = r0 + 2 + len(coh)
for j, h in enumerate(['Recurrence threshold', 'Observed genes', 'Random genes, mean of 1000 draws', 'Random genes, maximum', 'Permutation p'], 1):
    ws.cell(row=r1, column=j, value=h).font = fB
jc = rc['j_recurrence_comparator']
extra = {'≥ 5 cohorts, accessibility-informative genes (comparator)': {'observed': jc['accessibility']['genes_informative_in_>=5_of_10'], 'null_mean': jc['accessibility']['null_mean'], 'null_max': jc['accessibility']['null_max'], 'perm_p': float('nan')}}
for i, (k, v) in enumerate(list(ph['c_recurrence_vs_null'].items()) + list(extra.items())):
    for j, val in enumerate([k.replace('>=', '≥ ').replace('_', ' '), v['observed'], float(v['null_mean']), v['null_max'], float(v['perm_p'])], 1):
        x = ws.cell(row=r1 + 1 + i, column=j, value=val); x.font = f
        if isinstance(val, float): x.number_format = '0.000'

# ---------- S6 breast detail ----------
old = open(f'{BIJ}/manuscript/manuscript_source.md', encoding='utf-8').read()
def md_table(num):
    blk = old.split(f'**Table {num}.**')[1].split('**Table ')[0].split('## Figure legends')[0]; lines = [l for l in blk.split('\n') if l.startswith('|')]
    out = []
    for l in lines:
        cells = [x.strip() for x in l.strip('|').split('|')]
        if set(''.join(cells)) <= set('-: '): continue
        out.append(cells)
    note = [l for l in blk.split('\n')[1:] if l and not l.startswith('|') and not l.startswith('**')][0]
    return out, note
ws = wb.create_sheet('Table S3'); ws.sheet_view.showGridLines = False
ws['A1'] = 'Table S3. Breast discovery cohort (73 primary tumors; internal analysis plan committed on 16 November 2026), detailed results.'; ws['A1'].font = fT
ws['A2'] = ('Values as obtained in the breast analysis (archived as version 1.0.1, https://doi.org/10.5281/zenodo.24294127). The Basal-versus-luminal comparison (Part B) '
            'was pre-specified for breast, but its permutation test could not be evaluated; it is descriptive and is not part of the pan-cancer test.'); ws['A2'].font = fN
r = 4
for part, num in [('Part A. Median per-gene CV-R² (95% jackknife CI); models adjusted for array platform.', 2),
                  ('Part B. Promoter changes accompanying lower expression in Basal versus luminal tumors, among genes with both a methylation value and a promoter peak; n (%) [Wilson 95% CI].', 3)]:
    ws.cell(row=r, column=1, value=part).font = fB; r += 1
    rows_, note = md_table(num)
    for i, cells in enumerate(rows_):
        for j, v in enumerate(cells, 1):
            x = ws.cell(row=r, column=j, value=v); x.font = fB if i == 0 else f; x.alignment = wrap
            if i == 0: x.border = Border(top=thin, bottom=thin)
        r += 1
    ws.cell(row=r, column=1, value=note).font = fN; r += 2
ws.cell(row=r, column=1, value='Part C. Coupling between promoter methylation and accessibility (breast).').font = fB; r += 1
for line in ['1781 genes with variable promoter methylation: median Spearman ρ −0.19; 29.8% with ρ ≤ −0.3 (95% bootstrap CI 26.8–39.5%).',
             'CpG-island promoters (n = 1184): median ρ −0.17, 26.2% with ρ ≤ −0.3; non-island promoters (n = 597): −0.22, 37.0%.',
             'Re-run through the pan-cancer code (Table S1, pan-cancer entry 1): identical values for 1193 genes with and without the explicit probe mask '
             f"(D {b['masked_M_general']['D_median_A_minus_M']['estimate']:.4f}; G {b['masked_M_general']['G_median_MA_minus_A']['estimate']:.4f})."]:
    ws.cell(row=r, column=1, value=line).font = f; r += 1
for col, w in zip('ABCDEFG', [30, 22, 22, 22, 22, 22, 22]): ws.column_dimensions[col].width = w
wb._sheets = [wb['Contents']] + sorted([w for w in wb.worksheets if w.title != 'Contents'], key=lambda w: w.title)
wb.save(OUT); print('saved', OUT, [s.title for s in wb.worksheets])
