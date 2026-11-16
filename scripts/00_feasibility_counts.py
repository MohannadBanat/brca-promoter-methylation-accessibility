# Feasibility gate only: counts tumours with all three assays. No outcome data examined.
import gzip, csv, sys, os
D = os.path.join(os.path.dirname(__file__), '..', 'data')
atac = {}
for r in csv.DictReader(open(os.path.join(D,'sample_id_lookup.txt')), delimiter='\t'):
    if r['bam_prefix'].startswith('BRCA-'):
        s = r['Case_ID'][:15]  # TCGA-XX-XXXX-01
        atac.setdefault(s, 0); atac[s] += 1
with gzip.open(os.path.join(D,'BRCA_HiSeqV2.gz'),'rt') as f: rna = set(f.readline().rstrip('\n').split('\t')[1:])
with gzip.open(os.path.join(D,'BRCA_HumanMethylation450.gz'),'rt') as f: meth = set(f.readline().rstrip('\n').split('\t')[1:])
pam = {}
for r in csv.DictReader(open(os.path.join(D,'BRCA_clinicalMatrix.tsv')), delimiter='\t'):
    pam[r['sampleID']] = r.get('PAM50Call_RNAseq','')
prim = {s for s in atac if s.endswith('-01')}
allthree = sorted(s for s in prim if s in rna and s in meth)
from collections import Counter
print('ATAC BRCA samples (15-char):', len(atac), 'primary -01:', len(prim))
print('with RNA:', len([s for s in prim if s in rna]), 'with 450K:', len([s for s in prim if s in meth]))
print('all three assays:', len(allthree))
print('PAM50 among all-three:', dict(Counter(pam.get(s,'NA') or 'NA' for s in allthree)))
open(os.path.join(D,'cohort_allthree.txt'),'w').write('\n'.join(allthree)+'\n')
