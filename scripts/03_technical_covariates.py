"""Per-tumour technical covariates for negative control 1 (plan Section 6).
Run: python3 -I scripts/03_technical_covariates.py data/DataS1_samples.xlsx data/cohort_74_atac.txt processed/meth_shared_promoter_probes.tsv.gz processed/technical_covariates.tsv"""
import sys, numpy as np, pandas as pd
s1, cohort_f, meth_f, out = sys.argv[1:5]
d = pd.read_excel(s1, 'Sequencing Statistics', header=37)
d = d[d['cohort'] == 'BRCA'].copy()
d['case'] = d['Tissue_Barcode'].astype(str).str[:15]
d['log10_dedup_reads'] = np.log10(d['Final_DeDup_reads'].astype(float))
d['dup_fraction'] = d['Picard_Duplicate_reads'].astype(float) / d['Aligned_reads'].astype(float)
t = d.groupby('case')[['log10_dedup_reads', 'TSS_Enrichment_Score', 'FRIP', 'dup_fraction']].mean()
t['n_atac_reps'] = d.groupby('case').size()
m = pd.read_csv(meth_f, sep='\t', index_col=0)
t['meth_sample_mean_beta'] = m.mean(axis=0).reindex(t.index)   # global level across shared promoter probes (batch proxy)
cohort = [l.strip() for l in open(cohort_f) if l.strip()]
t = t.reindex(cohort); t.index.name = 'case'
t.drop(columns=['n_atac_reps']).to_csv(out, sep='\t')
print(t.describe().T.round(3).to_string()); print('missing:', t.isna().sum().to_dict())
