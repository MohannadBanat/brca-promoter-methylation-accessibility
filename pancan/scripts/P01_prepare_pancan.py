"""Pan-cancer input preparation (standard library only; runs on the Mac). Selects promoter-proximal data; computes no associations.
Streams every download and keeps only promoter-proximal rows and cohort columns, so no full raw file is stored.
Usage: python3 -I P01_prepare_pancan.py <IJMS_paper dir> <cohort_all3_primary.tsv> <gdc_methylation_pancan_manifest.tsv> <out dir>
"""
import bisect, csv, gzip, hashlib, io, json, os, re, subprocess, sys, time, zipfile, zlib
from collections import defaultdict

ROOT, COHORT, MANIFEST, OUT = sys.argv[1:5]
D = os.path.join(ROOT, 'data'); os.makedirs(OUT, exist_ok=True)
FLANK = 2000  # generous; the plan's window (-1000/+500) is applied in the cloud analysis
LOG = open(os.path.join(OUT, 'prepare_log.txt'), 'a')
def log(*a):
    s = time.strftime('%Y-%m-%dT%H:%M:%S ') + ' '.join(map(str, a)); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()

# 1. Protein-coding TSS windows (same as the breast pipeline)
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
        tss[k] = t if k not in tss else (min(tss[k], t) if strand == '+' else max(tss[k], t))
win = defaultdict(list)
for (gid, gn, c, st), t in tss.items(): win[c].append((t - FLANK, t + FLANK))
for c in win: win[c].sort()
starts = {c: [w[0] for w in win[c]] for c in win}
def near_tss(c, s, e):
    if c not in win: return False
    i = bisect.bisect_right(starts[c], e); j = bisect.bisect_left(starts[c], s - 2 * FLANK)
    return any(win[c][k][1] >= s and win[c][k][0] <= e for k in range(j, i))
log('protein-coding genes with TSS:', len(tss))

# 2. Cohort (non-breast; HM450 primary tumours with all three assays; breast already prepared)
coh = [r for r in csv.DictReader(open(COHORT), delimiter='\t') if r['type'] != 'BRCA' and r['hm450'] == 'True']
types = sorted({r['type'] for r in coh}); bytype = {t: [r['sample'] for r in coh if r['type'] == t] for t in types}
log('cohort tumours:', len(coh), 'types:', len(types))

def curl_stream(url, extra=()):
    return subprocess.Popen(['curl', '-sS', '-L', '--fail', '-m', '900', *extra, url], stdout=subprocess.PIPE)

# 3. RNA-seq: stream each type's HiSeqV2, keep cohort columns
for t in types:
    fn = os.path.join(OUT, f'rna_{t}.tsv.gz')
    if os.path.exists(fn): continue
    for attempt in range(3):
        p = curl_stream(f'https://tcga.xenahubs.net/download/TCGA.{t}.sampleMap/HiSeqV2.gz')
        sha = hashlib.sha256(); d = zlib.decompressobj(16 + zlib.MAX_WBITS); buf = ''; out_rows = []; hdr = None
        for chunk in iter(lambda: p.stdout.read(1 << 20), b''):
            sha.update(chunk); buf += d.decompress(chunk).decode()
            *lines, buf = buf.split('\n')
            for line in lines:
                x = line.split('\t')
                if hdr is None:
                    hdr = x; ii = [hdr.index(s) for s in bytype[t] if s in hdr]
                    out_rows.append('gene\t' + '\t'.join(hdr[i] for i in ii)); continue
                out_rows.append(x[0] + '\t' + '\t'.join(x[i] for i in ii))
        if buf.strip():
            x = buf.split('\t'); out_rows.append(x[0] + '\t' + '\t'.join(x[i] for i in ii))
        if p.wait() == 0 and hdr is not None:
            with gzip.open(fn, 'wt') as o: o.write('\n'.join(out_rows) + '\n')
            log('RNA', t, 'samples', len(ii), 'genes', len(out_rows) - 1, 'download sha256', sha.hexdigest()); break
        log('RNA retry', t, attempt)

# 4. ATAC: type-specific normalized counts, promoter-proximal peaks, replicates averaged per tumour
uuid2case = {}
for r in csv.DictReader(open(os.path.join(D, 'sample_id_lookup.txt')), delimiter='\t'):
    uuid2case[r['stanfordUUID'].replace('-', '_')] = r['Case_ID'][:15]
z = zipfile.ZipFile(os.path.join(D, 'atac_typespecific_normcounts.zip'))
with open(os.path.join(OUT, 'atac_replicate_map.tsv'), 'w') as om:
    om.write('type\treplicate\tcase\tbatch\n')
    for t in types:
        fn = os.path.join(OUT, f'atac_{t}.tsv.gz')
        name = [n for n in z.namelist() if n.startswith(t + '_')][0]
        with z.open(name) as fb:
            f = io.TextIOWrapper(fb)
            hdr = f.readline().rstrip('\n').split('\t'); meta = hdr[:5]; reps = hdr[5:]
            rep_case = []
            for rp in reps:
                u = '_'.join(rp.split('_')[1:6]); c = uuid2case.get(u, '')
                b = re.search(r'_(X\d+)_', rp); rep_case.append(c); om.write(f'{t}\t{rp}\t{c}\t{b.group(1) if b else ""}\n')
            cases = [c for c in bytype[t] if c in set(rep_case)]
            idx = {c: [i for i, rc in enumerate(rep_case) if rc == c] for c in cases}
            if os.path.exists(fn): continue
            n = 0
            with gzip.open(fn, 'wt') as o:
                o.write('\t'.join(meta) + '\t' + '\t'.join(cases) + '\n')
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    if not near_tss(x[0], int(x[1]), int(x[2])): continue
                    v = x[5:]
                    o.write('\t'.join(x[:5]) + '\t' + '\t'.join(f'{sum(float(v[i]) for i in idx[c]) / len(idx[c]):.4f}' for c in cases) + '\n'); n += 1
        log('ATAC', t, 'replicates', len(reps), 'tumours', len(cases), 'promoter peaks', n)

# 5. Methylation: all promoter-proximal HM450 probes (masks applied later, in the analysis)
with gzip.open(os.path.join(D, 'HM450.hg38.manifest.gencode.v36.tsv.gz'), 'rt') as f:
    m450 = {row['probeID']: row for row in csv.DictReader(f, delimiter='\t')}
keep = [p for p in sorted(m450) if m450[p].get('CpG_chrm') not in (None, '', 'NA', '*')
        and near_tss(m450[p]['CpG_chrm'], int(m450[p]['CpG_beg']), int(m450[p]['CpG_end']))]
keepset = set(keep); log('HM450 promoter-proximal probes:', len(keep))
with open(os.path.join(OUT, 'probe_annotation_hm450_promoter.tsv'), 'w') as o:
    cols = list(next(iter(m450.values())).keys()); o.write('\t'.join(cols) + '\n')
    for p in keep: o.write('\t'.join(m450[p].get(c, '') or '' for c in cols) + '\n')
man = defaultdict(list)
for r in csv.DictReader(open(MANIFEST), delimiter='\t'):
    if '450' in r['platform'] and r['sample_type'] == 'Primary Tumor': man[r['sample'][:15]].append(r)
md5 = {}
ids = sorted({sorted(man[s], key=lambda r: r['file_id'])[0]['file_id'] for t in types for s in bytype[t] if s in man})
body = json.dumps({"filters": {"op": "in", "content": {"field": "file_id", "value": ids}}, "fields": "file_id,md5sum,file_size", "size": 2000})
res = json.loads(subprocess.run(['curl', '-sS', '-X', 'POST', '-H', 'Content-Type: application/json', '--data-binary', body,
                                 'https://api.gdc.cancer.gov/files'], capture_output=True).stdout)
for h in res['data']['hits']: md5[h['file_id']] = h['md5sum']
done = os.path.join(OUT, 'meth_parts'); os.makedirs(done, exist_ok=True)
for t in types:
    for s in bytype[t]:
        if s not in man: log('no HM450 file', t, s); continue
        r = sorted(man[s], key=lambda r: r['file_id'])[0]; part = os.path.join(done, f'{t}__{s}.tsv.gz')
        if os.path.exists(part): continue
        for attempt in range(4):
            p = curl_stream(f'https://api.gdc.cancer.gov/data/{r["file_id"]}')
            h = hashlib.md5(); buf = b''; vals = []
            for chunk in iter(lambda: p.stdout.read(1 << 20), b''):
                h.update(chunk); buf += chunk
                *lines, buf = buf.split(b'\n')
                for line in lines:
                    x = line.decode().split('\t')
                    if len(x) >= 2 and x[0] in keepset: vals.append(f'{x[0]}\t{x[1]}')
            if buf:
                x = buf.decode().split('\t')
                if len(x) >= 2 and x[0] in keepset: vals.append(f'{x[0]}\t{x[1]}')
            if p.wait() == 0 and h.hexdigest() == md5.get(r['file_id']):
                with gzip.open(part, 'wt') as o: o.write('\n'.join(vals) + '\n')
                log('METH', t, s, r['file_id'], 'probes', len(vals), 'md5 ok'); break
            log('METH retry', t, s, attempt, 'md5', h.hexdigest(), 'expected', md5.get(r['file_id']))
            time.sleep(5)
log('done')
