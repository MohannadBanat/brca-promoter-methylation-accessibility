"""Pan-cancer feasibility gate: sample counts only. No methylation, accessibility or expression VALUES are read.
For each TCGA ATAC-seq tumor sample (15-character barcode, as in the breast analysis) it checks:
  RNA-seq  : sample present in the header of the Xena TCGA hub HiSeqV2 file of its cancer type (header line only)
  Methylation: a GDC harmonised 'Methylation Beta Value' file (HM450 or HM27) for the same 15-character sample
  Purity   : ABSOLUTE estimate present
Run on the Mac: python3 -I P00_feasibility_counts.py atac_samples.tsv ABSOLUTE_pancan_purity.txt out_dir"""
import csv, json, subprocess, sys, os, zlib
atac_tsv, purity_txt, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
S = list(csv.DictReader(open(atac_tsv), delimiter='\t'))
for r in S: r['s15'] = r['sample'][:15]
types = sorted({r['type'] for r in S})
def curl(args, data=None):
    cmd = ['curl', '-sS', '-m', '120', '--retry', '3'] + args
    return subprocess.run(cmd, input=data, capture_output=True).stdout
# RNA-seq: header line of each cohort file (stream, decompress only the start)
XENA = 'https://tcga.xenahubs.net/download/TCGA.{t}.sampleMap/HiSeqV2.gz'
rna = {}
for t in types:
    raw = curl(['-L', '-r', '0-400000', XENA.format(t=t)])
    try:
        txt = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw).decode()
        rna[t] = set(txt.split('\n')[0].split('\t')[1:])
    except Exception as e:
        rna[t] = set(); print('RNA header failed', t, e)
# Methylation: GDC files endpoint, by patient barcode
patients = sorted({r['patient'] for r in S})
filt = {"op": "and", "content": [
    {"op": "in", "content": {"field": "cases.submitter_id", "value": patients}},
    {"op": "=", "content": {"field": "data_type", "value": "Methylation Beta Value"}}]}
fields = 'file_id,file_name,platform,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.project.project_id'
body = json.dumps({"filters": filt, "fields": fields, "format": "JSON", "size": 5000}).encode()
res = json.loads(curl(['-X', 'POST', '-H', 'Content-Type: application/json', '--data-binary', '@-',
                       'https://api.gdc.cancer.gov/files'], body))
hits = res['data']['hits']
meth = {}
with open(os.path.join(out, 'gdc_methylation_pancan_manifest.tsv'), 'w') as fo:
    fo.write('case\tsample\tsample_type\tplatform\tfile_id\tfile_name\tproject\n')
    for h in hits:
        c = h['cases'][0]; smp = c['samples'][0]
        fo.write('\t'.join([c['submitter_id'], smp['submitter_id'], smp.get('sample_type', ''), h.get('platform', ''),
                            h['file_id'], h['file_name'], c['project']['project_id']]) + '\n')
        meth.setdefault(smp['submitter_id'][:15], set()).add(h.get('platform', ''))
pur = set()
for r in csv.DictReader(open(purity_txt), delimiter='\t'):
    pur.add(r['sample'][:15])
rows = []
for r in S:
    p = meth.get(r['s15'], set())
    rows.append(dict(type=r['type'], sample=r['s15'], primary=r['stype'] == '01', rna=r['s15'] in rna[r['type']],
                     hm450=any('450' in x for x in p), hm27=any('27' in x for x in p), purity=r['s15'] in pur))
with open(os.path.join(out, 'pancan_feasibility_samples.tsv'), 'w') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
print('GDC methylation files returned:', len(hits), '| total reported:', res['data']['pagination']['total'])
print('type\tATAC\tprimary\tRNA\tmeth_any\tHM450\tall3_primary\tall3+purity')
tot = [0] * 7
for t in types:
    g = [r for r in rows if r['type'] == t]
    a3 = [r for r in g if r['primary'] and r['rna'] and (r['hm450'] or r['hm27'])]
    v = [len(g), sum(r['primary'] for r in g), sum(r['rna'] for r in g), sum(r['hm450'] or r['hm27'] for r in g),
         sum(r['hm450'] for r in g), len(a3), sum(r['purity'] for r in a3)]
    tot = [x + y for x, y in zip(tot, v)]
    print(t + '\t' + '\t'.join(map(str, v)))
print('TOTAL\t' + '\t'.join(map(str, tot)))
