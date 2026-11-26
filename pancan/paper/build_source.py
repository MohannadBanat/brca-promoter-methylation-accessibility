"""Fill tables and references in the keyed manuscript, then check every quoted number against the result files.
Run: python3 -I build_source.py <pancan results dir> <ref dir> <keyed md> <out md> <version> <doi>"""
import json, re, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from refs import REFS
R, REF, SRC, OUT, VER, DOI = sys.argv[1:7]
p = json.load(open(f'{R}/pancan_results.json')); b = json.load(open(f'{R}/breast_rerun.json'))['masked_M_general']
s5 = json.load(open(f'{R}/S5_meta_with_breast.json')); ph = json.load(open(f'{R}/posthoc_robustness.json'))
nl = json.load(open(f'{R}/posthoc_single_predictor_null.json'))['per_type']
rc = json.load(open(f'{R}/posthoc_review_corrections.json'))
T9 = p['meta_types_included']; pt = p['per_type']
NAME = {'KIRP': 'Kidney papillary', 'PRAD': 'Prostate', 'COAD': 'Colon', 'LUAD': 'Lung adenocarcinoma', 'STAD': 'Stomach',
        'ESCA': 'Esophagus', 'LIHC': 'Liver', 'LUSC': 'Lung squamous', 'KIRC': 'Kidney clear cell'}
def n3(x): return f'{(0.0 if abs(x) < 0.0005 else x):.3f}'.replace('-', '−')
def ci(e): return f"{n3(e['estimate'])} ({n3(e['CI95'][0])}, {n3(e['CI95'][1])})"

# breast purity (descriptive)
coh = pd.read_csv(f'{REF}/cohort_all3_primary.tsv', sep='\t'); pur = pd.read_csv(f'{REF}/ABSOLUTE_pancan_purity.txt', sep='\t')
pur = pur.assign(s15=pur['sample'].astype(str).str[:15]).drop_duplicates('s15').set_index('s15')['purity']
bp = pur.reindex(coh.loc[coh.type == 'BRCA', 'sample'].str[:15]).astype(float)

t1 = []
for T in T9:
    v = pt[T]
    t1.append(f"| {NAME[T]} ({T}) | {v['n_tumours']} | {v['mean_purity']:.2f} | {v['genes']['variable']} ({v['genes']['variable_HCP']}) | "
              f"{n3(v['S1']['median_rho']).replace('0.', '0.')[:5] if False else n3(round(v['S1']['median_rho'], 2))[:-1]} | {100 * v['NC1_technical_only']['frac_genes_MA_gt_technical']:.1f}% |")
t1.append(f"| Breast, discovery (BRCA) | 73 | {np.nanmean(bp):.2f} | 1193 (747) | −0.19 | 95.2% |")
def base(n): return n3(1 - (n / (n - 1)) ** 2)
t2 = [f"| {NAME[T]} ({T}) | {base(pt[T]['n_tumours'])} | {ci(pt[T]['H2']['median_cvR2_M'])} | {ci(pt[T]['H2']['median_cvR2_A'])} | {ci(pt[T]['H2']['median_cvR2_MA'])} | "
      f"{ci(pt[T]['H2']['D_median_A_minus_M'])} | {ci(pt[T]['H2']['G_median_MA_minus_A'])} |" for T in T9]
t2.append(f"| Pooled, 9 test types | — | — | — | — | {ci(p['R1_pooled_D'])} | {ci(p['R3_pooled_G'])} |")
t2.append(f"| Breast, discovery (BRCA) | {base(73)} | {ci(b['median_cvR2_M'])} | {ci(b['median_cvR2_A'])} | {ci(b['median_cvR2_MA'])} | {ci(b['D_median_A_minus_M'])} | {ci(b['G_median_MA_minus_A'])} |")
s2 = p['S2_pooled_within_type']
rt = lambda g: f"{ci(g['median_rel_M'])} | {ci(g['median_rel_A'])} | {ci(g['median_rel_MA'])} | {ci(g['D_rel'])} | {ci(g['G_rel'])}"
t3r = [f"| {NAME[T]} ({T}) | " + rt(rc['per_type'][T]['g_relative']) + " |" for T in T9]
t3r.append(f"| Pooled, 9 test types | — | — | — | {ci(rc['g_pooled']['D_rel'])} | {ci(rc['g_pooled']['G_rel'])} |")
t3r.append("| Breast, discovery (BRCA) | " + rt(rc['BRCA']['g_relative']) + " |")
t3r.append(f"| All types, within type | {ci(s2['median_cvR2_M'])} | {ci(s2['median_cvR2_A'])} | {ci(s2['median_cvR2_MA'])} | {ci(s2['D_median_A_minus_M'])} | {ci(s2['G_median_MA_minus_A'])} |")
npur = sum(pt[T]['S4_purity']['n_tumours'] for T in T9)
pa = ph['a_pooled']
rows3 = [('R1, primary (pre-specified)', '9 types', ci(p['R1_pooled_D']), '—', 'Met'),
         ('R2, co-primary (pre-specified)', '9 types', f"D > 0 in {p['R2_sign']['types_D_positive']} of 9", '—', f"Not met (p = {p['R2_sign']['one_sided_p']:.3f})"),
         ('R3, co-primary (pre-specified)', '9 types', '—', ci(p['R3_pooled_G']), 'Met'),
         ('S2, all types pooled within type', f"17 types, {s2['n_tumours']} tumors", ci(s2['D_median_A_minus_M']), ci(s2['G_median_MA_minus_A']), '—'),
         ('S3, CpG-island promoters', '9 types', ci(p['S3_HCP_pooled']['D_median_A_minus_M']), ci(p['S3_HCP_pooled']['G_median_MA_minus_A']), '—'),
         ('S3, non-island promoters', '9 types', ci(p['S3_LCP_pooled']['D_median_A_minus_M']), ci(p['S3_LCP_pooled']['G_median_MA_minus_A']), '—'),
         ('S4, purity-adjusted', f"9 types, {npur} tumors", ci(p['S4_purity_pooled']['D_median_A_minus_M']), ci(p['S4_purity_pooled']['G_median_MA_minus_A']), '—'),
         ('S5, including breast', '10 cohorts', ci(s5['D_median_A_minus_M']), ci(s5['G_median_MA_minus_A']), f"D > 0 in {s5['D_positive_including_breast']} of 10"),
         ('S7, breast probe set', '9 types', ci(p['S7_shared_probes_pooled']['D_median_A_minus_M']), ci(p['S7_shared_probes_pooled']['G_median_MA_minus_A']), '—'),
         ('Post hoc: gene reselection', '9 types', ci(pa['D_median_A_minus_M']), ci(pa['G_median_MA_minus_A']), f"D > 0 in {ph['a_types_D_positive']} of 9"),
         ('Post hoc: HKSJ interval', '9 types', f"{n3(rc['i_R1_D']['estimate'])} ({n3(rc['i_R1_D']['HKSJ_CI95'][0])}, {n3(rc['i_R1_D']['HKSJ_CI95'][1])})",
          f"{n3(rc['i_R3_G']['estimate'])} ({n3(rc['i_R3_G']['HKSJ_CI95'][0])}, {n3(rc['i_R3_G']['HKSJ_CI95'][1])})", '—'),
         ('Post hoc: relative to type mean', '9 types', ci(rc['g_pooled']['D_rel']), ci(rc['g_pooled']['G_rel']), '—')]
t3 = ['| ' + ' | '.join(r) + ' |' for r in rows3]

src = open(SRC, encoding='utf-8').read()
src = src.replace('{{TABLE1}}', '\n'.join(t1)).replace('{{TABLE2}}', '\n'.join(t2)).replace('{{TABLE3}}', '\n'.join(t3r)).replace('{{TABLE4}}', '\n'.join(t3))
src = src.replace('{{VERSION}}', VER).replace('{{DOI}}', DOI)
# citations: number by first appearance (body before the reference list)
body = src.split('## References')[0]; order = []
for grp in re.findall(r'\[(@[^\]]+)\]', body):
    for k in re.findall(r'@([\w]+)', grp):
        if k not in order: order.append(k)
missing = [k for k in order if k not in REFS]; assert not missing, missing
unused = [k for k in REFS if k not in order]
def rep(m):
    nums = sorted(order.index(k) + 1 for k in re.findall(r'@([\w]+)', m.group(1)))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        out.append(f'{nums[i]}–{nums[j]}' if j - i >= 2 else (f'{nums[i]},{nums[j]}' if j == i + 1 else f'{nums[i]}')); i = j + 1
    return '[' + ','.join(out) + ']'
src = re.sub(r'\[(@[^\]]+)\]', rep, src)
src = src.replace('{{REFERENCES}}', '\n'.join(f'{i + 1}. {REFS[k]}' for i, k in enumerate(order)))
open(OUT, 'w', encoding='utf-8').write(src)

# ---------------- checks: every number quoted in the text against the outputs ----------------
txt = src.split('## Tables')[0]; fails = []
def has(s_, why):
    if s_ not in txt: fails.append(f'MISSING "{s_}"  ({why})')
def num(x, d=3): return f'{x:.{d}f}'.replace('-', '−')
def pc(x, d=1): return f'{100 * x:.{d}f}'.replace('-', '−')
e = lambda d: d['estimate']; lo = lambda d: d['CI95'][0]; hi = lambda d: d['CI95'][1]
rng_ = lambda d: f"{num(e(d))} ({num(lo(d))} to {num(hi(d))})"
# pre-specified
has(f"{num(e(p['R3_pooled_G']))} (95% CI {num(lo(p['R3_pooled_G']))} to {num(hi(p['R3_pooled_G']))}", 'R3')
has(f"{num(e(p['R1_pooled_D']))} (95% CI {num(lo(p['R1_pooled_D']))} to {num(hi(p['R1_pooled_D']))})", 'R1')
has(f"I² {round(100 * p['R1_pooled_D']['I2'])}%", 'R1 I2'); has(f"I² {round(100 * p['R3_pooled_G']['I2'])}%", 'R3 I2')
has(f"p = {p['R2_sign']['one_sided_p']:.3f}", 'R2 p'); assert p['R2_sign']['types_D_positive'] == 7; has('seven of the nine types', 'R2 count')
has(rng_(s2['G_median_MA_minus_A']), 'S2 G'); has(rng_(s2['D_median_A_minus_M']), 'S2 D')
has(f"{100 * e(s2['frac_A_gt_M']):.1f}% of genes", 'S2 frac')
has(f"{pc(e(s2['median_cvR2_M']))}% ({pc(lo(s2['median_cvR2_M']))}% to {pc(hi(s2['median_cvR2_M']))}%)", 'S2 M %')
has(f"{pc(e(s2['median_cvR2_A']))}% ({pc(lo(s2['median_cvR2_A']))}% to {pc(hi(s2['median_cvR2_A']))}%)", 'S2 A %')
has(f"{pc(e(s2['median_cvR2_A']))}% for accessibility", 'abstract A'); has(f"{pc(e(s2['median_cvR2_M']))}% pooled", 'abstract M')
has(f"gain was {round(100 * e(s2['G_median_MA_minus_A']), 1)} percentage points", 'discussion G pp')
# plan-scale M range
Ms = {T: e(pt[T]['H2']['median_cvR2_M']) for T in T9}; assert all(v < 0 for v in Ms.values())
has(f"from {num(min(Ms.values()))} in {min(Ms, key=Ms.get)} to {num(max(Ms.values()))} in {max(Ms, key=Ms.get)}", 'M range plan scale')
has('−0.063 for 33 tumors to −0.138 for 16', 'baseline'); assert num(1 - (33 / 32) ** 2) == '−0.063' and num(1 - (16 / 15) ** 2) == '−0.138'
# relative scale
g = {T: rc['per_type'][T]['g_relative'] for T in T9}
Mr = {T: e(g[T]['median_rel_M']) for T in T9}; Ar = {T: e(g[T]['median_rel_A']) for T in T9}
has(f"{pc(min(Mr.values()))}% ({min(Mr, key=Mr.get)}) to {pc(max(Mr.values()))}% ({max(Mr, key=Mr.get)})", 'M rel range')
has(f"{pc(min(Ar.values()))}% ({min(Ar, key=Ar.get)}) to {pc(max(Ar.values()))}% ({max(Ar, key=Ar.get)})", 'A rel range')
pos = [T for T in T9 if Mr[T] > 0]; assert pos == ['KIRP', 'PRAD', 'STAD', 'ESCA'], pos; has('four types (KIRP, PRAD, STAD and ESCA) and worse in five', 'M>0 types')
has(f"exceeded {num(max(hi(g[T]['median_rel_M']) for T in T9))}", 'M rel upper')
has(f"−{round(-100 * min(Mr.values()))}% to {round(100 * max(Mr.values()))}%", 'abstract/discussion rounded range')
br = rc['BRCA']['g_relative']['median_rel_M']; has(f"{pc(e(br))}% ({pc(lo(br))}% to {pc(hi(br))}%)", 'breast M rel')
has(f"accessibility explained {pc(Ar['PRAD'])}% and methylation {pc(Mr['PRAD'])}% in PRAD, and {pc(Ar['STAD'])}% and {pc(Mr['STAD'])}% in STAD", 'PRAD/STAD rel')
has(f"({pc(Ar['LIHC'])}% and {pc(Ar['LUSC'])}%)", 'LIHC/LUSC rel')
has(f"pooled D was {num(e(rc['g_pooled']['D_rel']))} ({num(lo(rc['g_pooled']['D_rel']))} to {num(hi(rc['g_pooled']['D_rel']))}) and pooled G {num(e(rc['g_pooled']['G_rel']))} ({num(lo(rc['g_pooled']['G_rel']))} to {num(hi(rc['g_pooled']['G_rel']))})", 'rel pooled')
# single-predictor null
nM = {T: nl[T]['null_median_M_mean'] for T in T9}; dM = {T: nl[T]['observed_median_M'] - nl[T]['null_median_M_mean'] for T in T9}
dA = {T: nl[T]['observed_median_A'] - nl[T]['null_median_A_mean'] for T in T9}
has(f"between {num(min(nM.values()))} ({min(nM, key=nM.get)}) and {num(max(nM.values()))} ({max(nM, key=nM.get)})", 'null M range')
has(f"by {num(min(dM.values()))} to {num(max(dM.values()))}", 'M above null'); has(f"by {num(min(dA.values()))} in {min(dA, key=dA.get)} to {num(max(dA.values()))} in {max(dA, key=dA.get)}", 'A above null')
# per-type G
for T in ['PRAD', 'STAD']:
    d = pt[T]['H2']['G_median_MA_minus_A']; has(f"{T} ({num(e(d))}; {num(lo(d))} to {num(hi(d))})", f'G {T}'); assert hi(d) > 0.02
neg = [T for T in T9 if hi(pt[T]['H2']['G_median_MA_minus_A']) < 0]; assert len(neg) == 7, neg; has('In seven types, G was below zero with a CI excluding zero', 'G neg 7')
has(rng_(b['G_median_MA_minus_A']), 'breast G')
for T in ['PRAD', 'STAD']:
    d = pt[T]['H2']['D_median_A_minus_M']; has(f"{T} (D {num(e(d))}; {num(lo(d))} to {num(hi(d))})" if T == 'PRAD' else f"{T} ({num(e(d))}; {num(lo(d))} to {num(hi(d))})", f'D {T}')
pg = {T: pd.read_csv(f'{R}/per_gene_{T}.csv') for T in T9}; beat = {T: 100 * ((pg[T].cvR2_M - pg[T].cvR2_A) > 0.1).mean() for T in T9}
oth = [beat[T] for T in T9 if T not in ('PRAD', 'STAD')]
has(f"{round(beat['PRAD'])}% of PRAD genes and {round(beat['STAD'])}% of STAD genes, against {round(min(oth))}–{round(max(oth))}%", 'M beats A')
has(rng_(s5['D_median_A_minus_M']), 'S5 D'); has(rng_(s5['G_median_MA_minus_A']), 'S5 G'); assert s5['D_positive_including_breast'] == 8
hk1, hk3 = rc['i_R1_D'], rc['i_R3_G']
has(f"prediction interval {num(hk1['prediction_interval95'][0])} to {num(hk1['prediction_interval95'][1])}", 'PI')
has(f"intervals were {num(hk1['HKSJ_CI95'][0])} to {num(hk1['HKSJ_CI95'][1])} for D and {num(hk3['HKSJ_CI95'][0])} to {num(hk3['HKSJ_CI95'][1])} for G", 'HKSJ')
# coupling
rho = {T: pt[T]['S1']['median_rho'] for T in T9}; fr = [pt[T]['S1']['frac_rho_le_-0.3'] for T in T9]
has(f"from {num(max(rho.values()), 2)} ({max(rho, key=rho.get)}) to {num(min(rho.values()), 2)} ({min(rho, key=rho.get)})", 'rho range')
has(f"from {100 * min(fr):.1f}% to {100 * max(fr):.1f}%", 'rho frac')
sp = rc['h_S1_pooled_DL']; has(f"{num(e(sp))} (95% CI {num(lo(sp))} to {num(hi(sp))}; I² {round(100 * sp['I2'])}%)", 'S1 pooled'); has(f"ρ {num(e(sp), 2)})", 'abstract rho')
nc1 = [pt[T]['NC1_technical_only']['frac_genes_MA_gt_technical'] for T in T9]; has(f"{round(100 * min(nc1))}–{round(100 * max(nc1))}% of genes", 'NC1')
assert all(pt[T]['NC2_mismatched']['observed_median_cvR2_MA'] > pt[T]['NC2_mismatched']['null_95pct'][1] for T in T9)
has(f"pooled D was {num(e(pa['D_median_A_minus_M']))} ({num(lo(pa['D_median_A_minus_M']))} to {num(hi(pa['D_median_A_minus_M']))}) and pooled G {num(e(pa['G_median_MA_minus_A']))} ({num(lo(pa['G_median_MA_minus_A']))} to {num(hi(pa['G_median_MA_minus_A']))})", 'post hoc a')
strata = [x for T in T9 for x in ph['per_type'][T]['b_probe_strata'].values()]; assert len(strata) == 27 and all(x['median_cvR2_M'] < 0 for x in strata)
ts = [ph['per_type'][T]['d_two_state'] for T in T9]; assert all(t['median_cvR2_two_state'] < t['median_cvR2_linear_M'] for t in ts) and all(t['median_gain_two_state_added_to_A'] < 0 for t in ts)
has(f"for {round(100 * min(t['frac_two_state_gt_linear'] for t in ts))}–{round(100 * max(t['frac_two_state_gt_linear'] for t in ts))}% of genes", 'two-state')
dd = [e(pt[T]['H2']['D_median_A_minus_M']) for T in T9]
has(f"ρ {num(stats.spearmanr([pt[T]['n_tumours'] for T in T9], dd)[0], 2)}, {num(stats.spearmanr([pt[T]['mean_purity'] for T in T9], dd)[0], 2)} and {ph['e_spearman_TSS_enrichment_vs_D']:.2f}", 'D associations')
ci_ = ph['c_methylation_informative']; fr2 = {T: 100 * ci_[T]['frac_M_gt_0.3'] for T in ci_}
has(f"{min(fr2.values()):.1f}% ({min(fr2, key=fr2.get)}) to {max(fr2.values()):.1f}% ({max(fr2, key=fr2.get)})", 'informative frac')
rv = ph['c_recurrence_vs_null']['>=5_cohorts']; has(f"against a mean of {rv['null_mean']:.1f} (maximum {rv['null_max']})", 'null'); assert rv['observed'] == 59
jc = rc['j_recurrence_comparator']['accessibility']; has(f"({jc['genes_informative_in_>=5_of_10']} genes against a mean of {jc['null_mean']:.1f})", 'acc comparator')
rec = pd.read_csv(f'{R}/posthoc_methylation_informative_recurrence.csv'); s_ = rec[rec['types_M_gt_0.3'] >= 5]
z = s_[s_.gene.str.match(r'^(ZNF|ZFP|ZIK|ZSCAN)')]; assert len(z) == 17 and (z.chrom == 'chr19').sum() == 12
for gg in ['DDX43', 'LDHC', 'OOEP', 'SPESP1', 'MKRN3', 'FAM50B', 'PTPRCAP', 'HCLS1', 'MYO1G', 'ARHGAP30']: assert gg in set(s_.gene), gg
named = 17 + 4 + 2 + 4; assert named < 59 / 2, named; assert int(rec.set_index('gene').loc['BST2', 'types_M_gt_0.3']) == 9 and rec.set_index('gene').loc['BST2', 'chrom'] == 'chr19'
assert sum(pt[T]['n_tumours'] for T in T9) == 191 and sum(pt[T]['genes']['variable'] for T in T9) == 22490 and s2['n_tumours'] == 266
assert 191 + 73 + (266 - 191) == 339; has('339 tumors were analyzed', 'total'); has('the 75 tumors of eight smaller types', 'n75')
vr = [pt[T]['genes']['variable'] for T in T9]; has(f"between {min(vr)} and {max(vr)} genes", 'gene range')
ss = src.split('## Simple Summary')[1].split('##')[0].split(); ab = src.split('## Abstract')[1].split('**Keywords')[0]
print('simple summary words', len(ss), '| abstract words', len(re.sub(r'\*\*[^*]+\*\*', '', ab).split()))
print('references', len(order), 'unused library keys', unused)
print('CHECKS:', 'ALL PASSED' if not fails else '\n'.join(fails))
