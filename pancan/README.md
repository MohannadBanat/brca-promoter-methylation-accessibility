# Pan-cancer extension (version 2.0.0)

Code, plan and results for:

> Abu Owida, H.; Alsarhan, A.; Banat, M.; Al-Fraihat, D.; Algarni, A. Promoter DNA Methylation Adds Little to Promoter Chromatin Accessibility in Explaining Gene Expression Differences Between Tumors: A Pre-Specified Test in Nine Cancer Types. Submitted to *Cancers* (2026).

**Study type:** secondary, observational analysis of public TCGA data. The breast analysis in the repository root is the discovery cohort; the nine other cancer types with at least 15 primary tumors with ATAC-seq, HM450 methylation and RNA-seq are the test.

## Plan and deviations

| File | Content |
|---|---|
| `plan/pancan_analysis_plan_2026-11-25.md` | Internal analysis plan committed on 25 November 2026 (00:03:36 UTC), before any non-breast methylation, accessibility or expression value was read |
| `plan/pancan_analysis_plan_2026-11-25.sha256` | SHA-256 `c440f092afb11f7e4afb9a7204ef51550ff12d0d88ad6f5c765527fbd542fa3b` |
| `plan/pancan_deviations_log.md` | Five entries. 1: written before any result; 2: after the primary results (reporting only); 3–5: post hoc analyses and corrections, each written after the results it responds to and before the analyses it adds were run |
| `plan/pancan_deviations_log.sha256` | Hash of the log as archived in this version |

## Pre-specified results (nine test types, 191 tumors)

| Criterion | Result | Met |
|---|---|---|
| R1: pooled D = median(CV-R²<sub>A</sub> − CV-R²<sub>M</sub>) > 0, CI excluding 0 | 0.130 (95% CI 0.053 to 0.207), I² 81% | Yes |
| R2: D > 0 in at least 8 of 9 types | 7 of 9 (PRAD and STAD negative), one-sided p = 0.090 | **No** |
| R3: upper 95% limit of pooled G = median(CV-R²<sub>MA</sub> − CV-R²<sub>A</sub>) < 0.02 | −0.046 (−0.060 to −0.032) | Yes |

All 266 non-breast tumors of 17 types pooled within type: D 0.123 (0.100 to 0.147), G 0.003 (−0.001 to 0.008), median CV-R² of methylation alone 0.021 (0.012 to 0.031).

Note on the CV-R² scale: as specified in the plan, CV-R² = 1 − PRESS/SST with SST about the mean of all tumors, so a model that predicts each held-out tumor by the mean of the others scores 1 − (n/(n − 1))² (−0.063 to −0.138 here), not zero. `P08_review_corrections.py` re-expresses per-type values relative to that leave-one-out mean (post hoc; the pooled within-type analysis is on that scale by design).

## Scripts (run order)

| Script | Where it ran | What it does |
|---|---|---|
| `P00_feasibility_counts.py` | local | sample counts only (feasibility gate) |
| `P01_prepare_pancan.py` | local, needs network | streams RNA-seq (Xena), ATAC-seq (GDC) and methylation (GDC, MD5-checked) and keeps promoter-proximal rows only |
| `P02_pancan_analysis.py` | any | the committed plan: per-type models, R1–R3, S1–S7, negative controls |
| `P03_breast_rerun.py` | any | breast through the pan-cancer code (S5), with and without the explicit probe mask |
| `P04_s5_meta_with_breast.py` | any | meta-analysis including breast (S5; deviation 2) |
| `P05_posthoc_robustness.py` | any | post hoc (deviation 3): reselection within jackknife, probe strata, two-state model, ATAC quality, methylation-informative genes |
| `P07_single_predictor_null.py` | any | post hoc (deviation 4): single-predictor mismatched-promoter null |
| `P08_review_corrections.py` | any | post hoc (deviation 5): CV-R² relative to the type mean, pre-specified random-effects pooling of S1, HKSJ and prediction intervals, accessibility comparator for recurrence |
| `P06_figures.py` | any | Figures 1–3 of the manuscript |

```bash
python3 -I pancan/scripts/P02_pancan_analysis.py <processed> <ref> <breast processed> pancan/results
python3 -I pancan/scripts/P03_breast_rerun.py <breast processed> <ref> pancan/results
python3 -I pancan/scripts/P04_s5_meta_with_breast.py pancan/results
python3 -I pancan/scripts/P05_posthoc_robustness.py <processed> <ref> pancan/results <breast per-gene csv>
python3 -I pancan/scripts/P07_single_predictor_null.py <processed> <ref> pancan/results
python3 -I pancan/scripts/P08_review_corrections.py <processed> <ref> pancan/results <breast processed> <breast per-gene csv>
python3 -I pancan/scripts/P06_figures.py pancan/results pancan/figures <breast per-gene csv>
```

`<ref>` holds `genes_tss_gencode_v36.tsv`, `HM450.hg38.mask.tsv.gz` (Zhou lab), `probe_annotation_shared.tsv`, `DataS1_samples.xlsx` (Corces et al.), `ABSOLUTE_pancan_purity.txt` and `cohort_all3_primary.tsv` (from P00). Random seeds are fixed. Re-running P02 from the processed inputs on 26 November 2026 reproduced `results/pancan_results.json` and all 21 `per_gene_*.csv` files byte for byte.

## Manuscript build

`paper/` holds the code that turns the results into the submitted files: `build_source.py` fills the tables from `results/`, numbers the citations and checks every number quoted in the text against the result files; `build_supplementary.py` writes Supplementary Tables S1–S6; `build_cancers_docx.py` builds the Word manuscript in the journal template; `manuscript_keyed.md` is the text with citation keys.
