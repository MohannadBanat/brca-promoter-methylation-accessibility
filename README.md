# Promoter DNA Methylation and Promoter Chromatin Accessibility as Explanations of Gene Expression in Primary Tumors (TCGA)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24293230.svg)](https://doi.org/10.5281/zenodo.24293230)

Code, committed analysis plans and results for two linked analyses:

1. **Breast discovery analysis** (repository root): 73 primary breast tumors; internal analysis plan committed on 16 November 2026.
2. **Pan-cancer test** (`pancan/`, added in version 2.0.0): nine further cancer types (191 tumors) under a second plan committed on 25 November 2026, before any non-breast data were read. See [`pancan/README.md`](pancan/README.md).

> Abu Owida, H.; Alsarhan, A.; Banat, M.; Al-Fraihat, D.; Algarni, A. Promoter DNA Methylation Adds Little to Promoter Chromatin Accessibility in Explaining Gene Expression Differences Between Tumors: A Pre-Specified Analysis of Ten Cancer Types. Submitted to *Cancers* (2026).

**Study type:** secondary, observational, integrative multi-omic analysis of public TCGA data. No new samples were collected.

## Archived versions

| Version | Content | DOI |
|---|---|---|
| 2.0.0 | Adds the pan-cancer test (`pancan/`); fixes the RNA-seq download address in `scripts/00_download_data.sh` (the Xena S3 address now returns HTTP 403; the replacement serves byte-identical files) | https://doi.org/10.5281/zenodo.24627804 |
| 1.0.1 | Breast analysis; American English spelling, figure sizes for the IJMS template | https://doi.org/10.5281/zenodo.24294127 |
| 1.0.0 | Breast analysis, first release | https://doi.org/10.5281/zenodo.24293231 |

The breast code, data manifests and numerical results are unchanged since 1.0.0. The DOI 10.5281/zenodo.24293230 always resolves to the latest version. The analysis plans and deviations logs keep their original spelling, because editing them would change their recorded SHA-256 hashes.

# Breast discovery analysis

## What the analysis does

For 73 primary breast tumors with matched ATAC-seq (Corces et al. 2018), DNA methylation (HM450 or HM27) and RNA-seq, it asks:

1. **H1:** how closely do promoter methylation and promoter accessibility track each other across tumors?
2. **H2:** how much cross-tumor expression variance does each explain, alone and together? Measured as leave-one-out cross-validated R² per gene.
3. **H3:** do genes expressed at lower levels in Basal than in luminal tumors show promoter hypermethylation, reduced accessibility, both or neither?

## Analysis plan and deviations

| File | Content |
|---|---|
| `plan/analysis_plan_2026-11-16.md` | Internal analysis plan committed on 16 November 2026, before any association was computed |
| `plan/analysis_plan_2026-11-16.sha256` | SHA-256 `b2b7079a364df1e80b5e9c5c0c46947204045ea7e401b5e1b1391ce3e8b48d2d` and commit time (13:48:48 UTC) |
| `plan/deviations_log.md` | Every change after the commit, with its timing. Entries 1–6 were written before results were viewed; entries 7–8 after |
| `plan/deviations_log_pre_results.sha256`, `plan/deviations_log_v2.sha256` | Hashes of the deviations log at those two points |

The plan was not registered on an external registry. Its hash was recorded locally at the time of writing.

## Reproducing the results

Requires Python ≥ 3.11 and about 2 GB of disk space.

```bash
pip install -r requirements.txt
bash scripts/00_download_data.sh                    # all public inputs into data/ (~1.6 GB)
python3 -I scripts/00_feasibility_counts.py         # feasibility gate: sample counts only
python3 -I scripts/01_prepare_inputs.py             # promoter-proximal subsets into processed/
python3 -I scripts/02_purity.py data processed      # ABSOLUTE purity for the cohort
python3 -I scripts/03_technical_covariates.py data/DataS1_samples.xlsx data/cohort_74_atac.txt \
    processed/meth_shared_promoter_probes.tsv.gz processed/technical_covariates.tsv
NPERM=1000 python3 -I scripts/04_analysis.py processed results 1000   # H1, H2, H3, negative controls
python3 -I scripts/05_jackknife_h2.py results       # jackknife CIs for H2 (deviation 7)
python3 -I scripts/05b_h3_both_measures.py results  # Table 3 proportions (deviation 8)
python3 -I scripts/06_figures.py results figures    # Figures 1–3
```

Random seeds are fixed (20261116). Starting from the processed inputs, a clean rerun reproduces every file in `results/` exactly. Checksums of all downloaded inputs are in `data/SHA256SUMS_*.txt`.

Tested with Python 3.13.18, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, statsmodels 0.15.0, Matplotlib 3.11.2 and openpyxl 3.1.5.

## Data sources

All inputs are public, and the download script fetches them directly from these sources:

- **ATAC-seq:** breast-cancer peak calls, normalized counts and sequencing statistics, from the [GDC ATACseq-AWG publication page](https://gdc.cancer.gov/about-data/publications/ATACseq-AWG) (Corces et al., *Science* 2018).
- **DNA methylation:** GDC harmonized SeSAMe level-3 beta values (HM450 and HM27), one primary-tumor file per case; the file IDs are in `data/gdc_methylation_all74_manifest.tsv`.
- **RNA-seq and PAM50:** the UCSC Xena TCGA BRCA hub.
- **Tumor purity:** ABSOLUTE estimates from the TCGA PanCancer Atlas.
- **Annotation:** Zhou lab HM450/HM27 GENCODE v36 manifests, and GENCODE v36 basic annotation.

## Main results

Pooled across 73 tumors, adjusted for array platform; 95% CIs from a delete-one-tumor jackknife:

| | Median per-gene CV-R² (1193 genes) |
|---|---|
| Methylation only | 0.005 (−0.025 to 0.036) |
| Accessibility only | 0.196 (0.134 to 0.258) |
| Both | 0.219 (0.160 to 0.277) |
| Gain from adding methylation to accessibility | −0.004 (−0.011 to 0.002) |

The full results are in `results/results.json` and `results/H2_jackknife.json`.

## Repository layout

```
pancan/    pan-cancer test: plan, deviations log, scripts, results, figures, manuscript build (see pancan/README.md)
scripts/   breast analysis code, in run order
plan/      analysis plan, deviations log and their hashes
results/   all numerical outputs
figures/   Figures 1–3 (600 dpi PNG)
data/      download manifest, cohort list and checksums (raw data are not redistributed)
```

## Use of generative AI

Claude (Anthropic) was used to help design both analysis plans, write and run the code, prepare the figures and draft the manuscripts. The authors reviewed the analyses and take full responsibility for them.

## Licence

MIT (see `LICENSE`).
