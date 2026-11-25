# Pan-cancer deviations log

Plan: pancan_analysis_plan_2026-11-25.md (SHA-256 c440f092afb11f7e4afb9a7204ef51550ff12d0d88ad6f5c765527fbd542fa3b, committed 2026-11-25T00:03:36Z).



## Entries written before any pan-cancer result was viewed

1. **Probe masks (implementation detail, 2026-11-25).** The plan says masks come from the Zhou lab GENCODE v36 manifests, but those files have no mask columns. Masks are taken instead from the Zhou lab HM450 hg38 mask file (`HM450.hg38.mask.tsv.gz`, column `M_general`; SHA-256 5de40616d0f6b400204663583bdadfc2cdbe8c57b9c619195f385f258e48c416). Checked on the breast data only: the GDC SeSAMe beta values are already missing for 99.3% of the 388 masked promoter probes (383 of them missing in every tumor), so the breast analysis effectively excluded masked probes. Masked probes are removed explicitly in the pan-cancer analysis; in the breast re-run (S5) they are removed too, and the effect on breast values is reported.
