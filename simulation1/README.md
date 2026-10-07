# Simulation 1 reproducibility files

This directory reproduces Simulation 1 of the MCNFA-CR manuscript. The study compares MFA-CR and MCNFA-CR under four censoring levels using 80 independent Monte Carlo replications per design cell.

In the R code and CSV files, the abbreviations `MFAC` and `MCNFAC` correspond to MFA-CR and MCNFA-CR, respectively.

## Files

- `00_run_simulation1.R`: master script for the complete simulation and table reconstruction.
- `01_simulation1_analysis.R`: simulation design, data generation, model fitting, and result export.
- `02_reproduce_table1.R`: validates raw results and reconstructs the manuscript table.
- `MCNFAC_Functions.R`: AECM fitting functions used by the simulation.
- `published_results/Study1_recovery_raw.csv`: archived row-level results used for the published table (1,920 rows).
- `published_results/MCNFAC_Lin_FINAL_recovery_summary.csv`: archived parameter-recovery summary.
- `published_results/MCNFAC_Lin_FINAL_clustering_summary.csv`: archived clustering summary.
- `published_results/Study1_FINAL_manuscript_table.csv`: archived table reported in the manuscript.

Files created by a new full run are written to `generated_results/simulation1/`, which is excluded from version control.

## Requirements

Use a current R installation with these packages installed:

```r
install.packages(c("mvtnorm", "MomTrunc", "mclust", "gtools"))
```

The scripts deliberately do not install packages automatically.

## Fast verification of the archived published results

From this directory, run:

```bash
Rscript 02_reproduce_table1.R published_results/Study1_recovery_raw.csv
```

This checks that:

- the raw file contains exactly 1,920 rows;
- each combination of sample size, censoring level, and model contains replications 1 through 80 exactly once; and
- all reconstructed manuscript entries agree with the archived table after the reported rounding.

The reconstructed table and validation report are written to `validation_output/`.

## Complete reproduction

Run:

```bash
Rscript 00_run_simulation1.R
```

The master script fits all Simulation 1 design cells, saves the raw and summary results, reconstructs the manuscript table, and records `sessionInfo()`. Because the simulation performs many censored-mixture fits, a full run can take substantial time.

The random-number seed is fixed at 123. The code uses a PSOCK cluster and leaves one detected CPU core free. To run serially, set `USE_PARALLEL <- FALSE` in `01_simulation1_analysis.R`.

## Reporting conventions

Summary metrics use the finite fitted values available for each metric. `Conv` is the proportion of the 80 replications that converged. The CPU time stored for a replication is the elapsed time of its retained (winning) initialization, rather than the cumulative time over all starting values. CPU times are therefore hardware- and workload-dependent and need not match when the study is rerun on another computer.

