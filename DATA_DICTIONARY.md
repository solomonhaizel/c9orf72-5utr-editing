# Data dictionary

Conditions are machine-readable as `0.3 uM` / `1 uM` / `scramble`; the `µM`
symbol is used only for display in figures.

## data/editing_by_rep.csv  (source of truth for editing values)
| column | meaning |
|---|---|
| construct | A74, A88, or ACTB |
| site | edited adenosine (A74, A88, ACTB) |
| dose | ASO dose: `0.3 uM`, `1 uM`, or `scramble` |
| rep | biological replicate number |
| run | sequencing run label (see data/runs.csv for the Azenta report ID) |
| **guide_raw_pct** | **RAW** editing 100·G/(A+G) at the target in the guide sample |
| guide_cov | reads at the target in the guide sample |
| noguide_raw_pct | RAW editing at the target in the matched no-guide sample (the baseline) |
| noguide_cov | reads at the target in the no-guide sample |
| **net_pct** | **BACKGROUND-SUBTRACTED** editing = guide_raw_pct − noguide_raw_pct |

## data/bystander_by_rep.csv  (amplicon-wide adenosine scan)

> Note: this scan shows a representative subset of replicates (a per-position QC view) rather than every sequenced sample, so its replicate count can be smaller than that of data/editing_by_rep.csv and the SRA deposit; quantitative editing values and replicate counts come from editing_by_rep.csv.

| column | meaning |
|---|---|
| condition | A74/A88 guide at each dose, No guide, or Scramble |
| rep | biological replicate |
| run | sequencing run label (see data/runs.csv) |
| position | biological 5′UTR position (1–98) of a reference adenosine |
| edit_pct | RAW editing 100·G/(A+G) at that position (no baseline subtraction) |
| site_type | `target` (A74/A88) or `other` |

## data/runs.csv  (sequencing-run manifest)
| column | meaning |
|---|---|
| run_id | run label used in the tables |
| azenta_report | Azenta 30-xxxx report number (stable sequencing-run ID) |
| construct | construct(s) in that run |
| replicates | which biological replicate(s) came from that run |
| note | provenance notes |

## results/stats_table.csv
| column | meaning |
|---|---|
| endpoint | construct + dose tested |
| comparison | guide vs matched no-guide (paired difference vs 0) |
| n | biological replicates |
| paired | yes (net values are paired differences) |
| test | test used, or "none (n<3, descriptive only)" |
| mean_net, sd | mean and SD of net editing |
| statistic, df | t statistic and degrees of freedom |
| p_raw, p_adj | raw and Benjamini–Hochberg-adjusted P |
| cohens_dz | paired effect size (mean/SD of differences) |
| ci95_low, ci95_high | 95% CI on mean net editing |
