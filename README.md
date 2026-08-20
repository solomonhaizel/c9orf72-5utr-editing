# C9ORF72 5′UTR RNA-editing: figure source

Data and code to regenerate the A-to-I editing figures (A74, A88, and the ACTB
positive control) for the C9ORF72 5′UTR manuscript. RNA-editing (amplicon-seq)
only; luminescence data are provided separately.

## How the values are computed (needed to reproduce)

- Editing = 100 × G / (A + G) at the target adenosine (A-to-I read as A→G).
- Net = guide − matched same-run no-guide (UTP), per replicate.
- Coordinates: A74 = ref position 107 (biological 74); A88 = ref 121
  (biological 88); ACTB target = position 113. The bystander scan uses the
  biological 1–98 axis.
- Negative net values are retained (not truncated).
- Statistics (`run_stats.R` → `results/stats_table.csv`): net values are paired
  differences, so the one-sample t vs 0 equals a paired t vs the matched
  no-guide; P values are Benjamini–Hochberg adjusted across the guide-vs-control
  tests; groups with n < 3 are reported descriptively (no test).

## Regenerate

```bash
python pipeline/make_editing_table.py     # -> data/editing_by_rep.csv
python pipeline/make_bystander_table.py   # -> data/bystander_by_rep.csv
Rscript analysis/run_stats.R              # -> results/stats_table.csv
Rscript analysis/plot_editing.R           # -> {A74,A88,ACTB}_editing.{pdf,png}
Rscript analysis/plot_coverage.R          # -> coverage.{pdf,png}
Rscript analysis/plot_bystander.R         # -> bystander_scan.{pdf,png}
```

`make_editing_table.py` and `make_bystander_table.py` read each run's Genotyper
`.var` output, regenerated from the SRA reads.

## Layout

```
pipeline/  analyze_c9_editing.py, analyze_actb_editing.py,
           make_editing_table.py, make_bystander_table.py
analysis/  run_stats.R, plot_editing.R, plot_coverage.R, plot_bystander.R
data/      editing_by_rep.csv, bystander_by_rep.csv, runs.csv
           (column definitions in DATA_DICTIONARY.md)
results/   stats_table.csv
```

## Data

Raw reads: SRA BioProject **PRJNA1514796**. Figure-input tables
are in `data/`.

## License

Code: MIT (`LICENSE`). Citation: `CITATION.cff`.
