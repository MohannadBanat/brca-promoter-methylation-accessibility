# Promoter Chromatin Accessibility Outperforms Promoter DNA Methylation in Explaining Gene Expression in Primary Breast Tumours

Code, analysis plan and results for:

> Abu Owida, H.; Alsarhan, A.; Banat, M.; Al-Fraihat, D.; Algarni, A. Promoter Chromatin Accessibility Outperforms Promoter DNA Methylation in Explaining Gene Expression in Primary Breast Tumours: An Integrative Multi-Omic Analysis of TCGA Data. Submitted to *International Journal of Molecular Sciences* (2026).

**Study type:** secondary, observational, integrative multi-omic analysis of public TCGA data. No new samples were collected.

## What the analysis does

For 73 primary breast tumours with matched ATAC-seq (Corces et al. 2018), DNA methylation (HM450 or HM27) and RNA-seq, it asks:

1. **H1:** how closely do promoter methylation and promoter accessibility track each other across tumours?
2. **H2:** how much cross-tumour expression variance does each explain, alone and together? Measured as leave-one-out cross-validated R² per gene.
3. **H3:** do genes expressed at lower levels in Basal than in luminal tumours show promoter hypermethylation, reduced accessibility, both or neither?

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

- **ATAC-seq:** breast-cancer peak calls, normalised counts and sequencing statistics, from the [GDC ATACseq-AWG publication page](https://gdc.cancer.gov/about-data/publications/ATACseq-AWG) (Corces et al., *Science* 2018).
- **DNA methylation:** GDC harmonised SeSAMe level-3 beta values (HM450 and HM27), one primary-tumour file per case; the file IDs are in `data/gdc_methylation_all74_manifest.tsv`.
- **RNA-seq and PAM50:** the UCSC Xena TCGA BRCA hub.
- **Tumour purity:** ABSOLUTE estimates from the TCGA PanCancer Atlas.
- **Annotation:** Zhou lab HM450/HM27 GENCODE v36 manifests, and GENCODE v36 basic annotation.

## Main results

Pooled across 73 tumours, adjusted for array platform; 95% CIs from a delete-one-tumour jackknife:

| | Median per-gene CV-R² (1193 genes) |
|---|---|
| Methylation only | 0.005 (−0.025 to 0.036) |
| Accessibility only | 0.196 (0.134 to 0.258) |
| Both | 0.219 (0.160 to 0.277) |
| Gain from adding methylation to accessibility | −0.004 (−0.011 to 0.002) |

The full results are in `results/results.json` and `results/H2_jackknife.json`.

## Repository layout

```
scripts/   analysis code, in run order
plan/      analysis plan, deviations log and their hashes
results/   all numerical outputs
figures/   Figures 1–3 (600 dpi PNG)
data/      download manifest, cohort list and checksums (raw data are not redistributed)
```

## Use of generative AI

Claude (Anthropic) was used to help design the analysis plan, write and run the code, and prepare the figures. The authors reviewed the analyses and take full responsibility for them.

## Licence

MIT (see `LICENSE`).
