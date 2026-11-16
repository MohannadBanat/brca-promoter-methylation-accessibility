"""ABSOLUTE purity (TCGA PanCancer Atlas) for the cohort, used in the purity-adjusted H2 sensitivity analysis.
Run: python3 -I scripts/02_purity.py data processed"""
import sys, pandas as pd
D, P = sys.argv[1:3]
a = pd.read_csv(f'{D}/ABSOLUTE_pancan_purity.txt', sep='\t')
a = a[a['call status'] == 'called'].drop_duplicates('array').set_index('array')
cohort = [l.strip() for l in open(f'{D}/cohort_74_atac.txt') if l.strip()]
p = a['purity'].reindex(cohort)
pd.DataFrame({'purity': p}).rename_axis('case').to_csv(f'{P}/purity_absolute.tsv', sep='\t')
print('with ABSOLUTE purity:', int(p.notna().sum()), 'of', len(cohort))
