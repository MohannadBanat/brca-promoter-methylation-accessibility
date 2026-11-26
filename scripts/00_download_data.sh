#!/usr/bin/env bash
# Download all public inputs into ./data (run from the repository root). About 1.6 GB.
# GDC sometimes resets long transfers; the large ATAC archive is fetched in 20 MB ranges.
set -euo pipefail
mkdir -p data/gdc_meth_all && cd data
get(){ for i in 1 2 3 4 5; do curl -fsSL --retry 3 -o "$1" "$2" && return 0; sleep 5; done; echo "FAILED $1" >&2; return 1; }
chunked(){ out=$1; url=$2; total=$3; step=20000000; mkdir -p ".parts_$out"; s=0; i=0
  while [ $s -lt $total ]; do e=$((s+step-1)); [ $e -ge $total ] && e=$((total-1)); p=".parts_$out/$(printf %04d $i)"
    for t in $(seq 1 10); do [ -f "$p" ] && [ "$(wc -c < "$p")" -eq $((e-s+1)) ] && break; curl -sS -m 300 -r $s-$e -o "$p" "$url" || true; sleep 2; done
    s=$((e+1)); i=$((i+1)); done
  cat ".parts_$out"/* > "$out" && rm -r ".parts_$out"; }
GDC=https://api.gdc.cancer.gov/data
# The Xena S3 address used for v1.0.x (tcga-xena-hub.s3.us-east-1.amazonaws.com) returned HTTP 403 from November 2026.
# The hub address below serves byte-identical files (SHA-256 of HiSeqV2.gz, BRCA_clinicalMatrix and HumanMethylation450.gz checked against data/SHA256SUMS_downloads.txt on 2026-11-26).
XENA=https://tcga.xenahubs.net/download
# Corces et al. 2018 (GDC ATACseq-AWG publication page)
chunked atac_typespecific_normcounts.zip $GDC/38b8f311-f3a4-4746-9829-b8e3edb9c157 632769396
get atac_typespecific_peaks.zip   $GDC/71ccfc55-b428-4a04-bb5a-227f7f3bf91c
get sample_id_lookup.txt          $GDC/7a3d7067-09d6-4acf-82c8-a1a81febf72c
get DataS1_samples.xlsx           $GDC/148b8d1e-a132-48f6-8e40-9ca96aa2f9b2
# TCGA PanCancer Atlas ABSOLUTE purity
get ABSOLUTE_pancan_purity.txt    $GDC/4f277128-f793-4354-a13d-30cc7fe9f6b5
# UCSC Xena TCGA BRCA hub: RNA-seq (RSEM log2) and clinical matrix (PAM50)
get BRCA_HiSeqV2.gz               "$XENA/TCGA.BRCA.sampleMap/HiSeqV2.gz"
get BRCA_clinicalMatrix.tsv       "$XENA/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
# Legacy Xena HM450 matrix: used only by 00_feasibility_counts.py (sample identifiers), not by the analysis
get BRCA_HumanMethylation450.gz   "$XENA/TCGA.BRCA.sampleMap/HumanMethylation450.gz"
# Probe and gene annotation
get HM450.hg38.manifest.gencode.v36.tsv.gz https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz
get HM27.hg38.manifest.gencode.v36.tsv.gz  https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/HM27.hg38.manifest.gencode.v36.tsv.gz
get gencode.v36.basic.annotation.gtf.gz    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.basic.annotation.gtf.gz
# GDC harmonised SeSAMe methylation betas, one primary-tumour file per case (manifest in data/)
tail -n +2 gdc_methylation_all74_manifest.tsv | while IFS=$'\t' read -r c s st pl id fn; do
  get "gdc_meth_all/${s}__${pl// /_}__${fn}" "$GDC/$id"; done
echo "Done. Compare against data/SHA256SUMS_*.txt"
