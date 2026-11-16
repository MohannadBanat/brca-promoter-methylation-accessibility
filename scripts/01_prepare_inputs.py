# Size reduction (standard library only). Selects promoter-proximal data; computes no associations.
# Outputs go to ../processed/.
import bisect, csv, glob, gzip, io, os, re, zipfile
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
D = os.path.join(ROOT, 'data'); P = os.path.join(ROOT, 'processed'); os.makedirs(P, exist_ok=True)
FLANK = 2000  # generous; the plan's window (-1000/+500) is applied in the cloud analysis

# 1. Protein-coding TSS (most 5' TSS of basic transcripts), GENCODE v36
tss = {}
with gzip.open(os.path.join(D, 'gencode.v36.basic.annotation.gtf.gz'), 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        x = line.rstrip('\n').split('\t')
        if x[2] != 'transcript': continue
        a = x[8]
        if 'gene_type "protein_coding"' not in a or 'tag "basic"' not in a: continue
        gid = re.search(r'gene_id "([^"]+)"', a).group(1); gn = re.search(r'gene_name "([^"]+)"', a).group(1)
        chrom, s, e, strand = x[0], int(x[3]), int(x[4]), x[6]
        t = s if strand == '+' else e
        k = (gid, gn, chrom, strand)
        if k not in tss: tss[k] = t
        else: tss[k] = min(tss[k], t) if strand == '+' else max(tss[k], t)
with open(os.path.join(P, 'genes_tss_gencode_v36.tsv'), 'w') as o:
    o.write('gene_id\tgene_name\tchrom\tstrand\ttss\n')
    for (gid, gn, c, st), t in sorted(tss.items(), key=lambda z: (z[0][2], z[1])): o.write(f'{gid}\t{gn}\t{c}\t{st}\t{t}\n')
win = defaultdict(list)
for (gid, gn, c, st), t in tss.items(): win[c].append((t - FLANK, t + FLANK))
for c in win: win[c].sort()
starts = {c: [w[0] for w in win[c]] for c in win}
def near_tss(c, s, e):
    if c not in win: return False
    i = bisect.bisect_right(starts[c], e); j = bisect.bisect_left(starts[c], s - 2 * FLANK)
    return any(win[c][k][1] >= s and win[c][k][0] <= e for k in range(j, i))
print('protein-coding genes with TSS:', len(tss))

# 2. Cohort and ATAC replicate -> tumour map
cohort = [l.strip() for l in open(os.path.join(D, 'cohort_74_atac.txt')) if l.strip()]
uuid2case = {}
for r in csv.DictReader(open(os.path.join(D, 'sample_id_lookup.txt')), delimiter='\t'):
    if r['bam_prefix'].startswith('BRCA-'): uuid2case[r['stanfordUUID'].replace('-', '_')] = r['Case_ID'][:15]
z = zipfile.ZipFile(os.path.join(D, 'atac_typespecific_normcounts.zip'))
name = [n for n in z.namelist() if 'BRCA' in n][0]; print('ATAC counts file:', name)
with z.open(name) as fb, gzip.open(os.path.join(P, 'atac_brca_promoter_peaks_by_tumour.tsv.gz'), 'wt') as o, \
     open(os.path.join(P, 'atac_replicate_map.tsv'), 'w') as om:
    f = io.TextIOWrapper(fb)
    hdr = f.readline().rstrip('\n').split('\t'); meta = hdr[:5]; reps = hdr[5:]
    rep_case = []; om.write('replicate\tcase\tbatch\n')
    for rp in reps:
        u = '_'.join(rp.split('_')[1:6]); c = uuid2case.get(u, '')
        b = re.search(r'_(X\d+)_', rp); rep_case.append(c); om.write(f'{rp}\t{c}\t{b.group(1) if b else ""}\n')
    cases = [c for c in cohort if c in set(rep_case)]
    idx = {c: [i for i, rc in enumerate(rep_case) if rc == c] for c in cases}
    o.write('\t'.join(meta) + '\t' + '\t'.join(cases) + '\n'); n = 0
    for line in f:
        x = line.rstrip('\n').split('\t')
        if not near_tss(x[0], int(x[1]), int(x[2])): continue
        v = x[5:]
        o.write('\t'.join(x[:5]) + '\t' + '\t'.join(f'{sum(float(v[i]) for i in idx[c]) / len(idx[c]):.4f}' for c in cases) + '\n'); n += 1
print('ATAC replicates:', len(reps), 'tumours mapped:', len(cases), 'promoter-proximal peaks kept:', n)

# 3. Methylation: probes shared by HM27 and HM450, promoter-proximal, all tumours with a file
def manifest(fn):
    with gzip.open(os.path.join(D, fn), 'rt') as f:
        return {row['probeID']: row for row in csv.DictReader(f, delimiter='\t')}
m450 = manifest('HM450.hg38.manifest.gencode.v36.tsv.gz'); m27 = manifest('HM27.hg38.manifest.gencode.v36.tsv.gz')
shared = sorted(set(m450) & set(m27))
keep = [p for p in shared if m450[p].get('CpG_chrm') not in (None, '', 'NA', '*')
        and near_tss(m450[p]['CpG_chrm'], int(m450[p]['CpG_beg']), int(m450[p]['CpG_end']))]
print('shared probes:', len(shared), 'promoter-proximal:', len(keep))
with open(os.path.join(P, 'probe_annotation_shared.tsv'), 'w') as o:
    cols = list(next(iter(m450.values())).keys()); o.write('\t'.join(cols) + '\n')
    for p in keep: o.write('\t'.join(m450[p].get(c, '') or '' for c in cols) + '\n')
keepset = set(keep); beta = {}; plat = {}
for fn in sorted(glob.glob(os.path.join(D, 'gdc_meth_all', '*'))):
    base = os.path.basename(fn); sample, platform = base.split('__')[0], base.split('__')[1]
    vals = {}
    for line in open(fn):
        x = line.rstrip('\n').split('\t')
        if len(x) >= 2 and x[0] in keepset:
            try: vals[x[0]] = float(x[1])
            except ValueError: pass
    beta[sample[:15]] = vals; plat[sample[:15]] = platform
ms = [c for c in cohort if c in beta]
with gzip.open(os.path.join(P, 'meth_shared_promoter_probes.tsv.gz'), 'wt') as o:
    o.write('probeID\t' + '\t'.join(ms) + '\n')
    for p in keep: o.write(p + '\t' + '\t'.join(f"{beta[c][p]:.4f}" if p in beta[c] else 'NA' for c in ms) + '\n')
with open(os.path.join(P, 'meth_platform.tsv'), 'w') as o:
    o.write('case\tplatform\n'); [o.write(f'{c}\t{plat[c]}\n') for c in ms]
print('methylation tumours:', len(ms))

# 4. RNA-seq and clinical subset
with gzip.open(os.path.join(D, 'BRCA_HiSeqV2.gz'), 'rt') as f, gzip.open(os.path.join(P, 'rna_rsem_log2_cohort.tsv.gz'), 'wt') as o:
    hdr = f.readline().rstrip('\n').split('\t'); ii = [hdr.index(c) for c in cohort if c in hdr]
    o.write('gene\t' + '\t'.join(hdr[i] for i in ii) + '\n')
    for line in f:
        x = line.rstrip('\n').split('\t'); o.write(x[0] + '\t' + '\t'.join(x[i] for i in ii) + '\n')
with open(os.path.join(D, 'BRCA_clinicalMatrix.tsv')) as f, open(os.path.join(P, 'clinical_cohort.tsv'), 'w') as o:
    r = csv.reader(f, delimiter='\t'); hdr = next(r); o.write('\t'.join(hdr) + '\n'); cs = set(cohort)
    for row in r:
        if row[0] in cs: o.write('\t'.join(row) + '\n')
print('done')
