# MCNFA-CR

This repository contains the data and R code used to reproduce the simulation study and the Lake Michigan application presented in the manuscript:

> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression**  
> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.

## Repository structure

- [`simulation1/`](simulation1/): R scripts and archived results for Simulation Study 1.
- [`real_data/`](real_data/): source data, R scripts, archived fitted-model results, tables, and figures for the Lake Michigan application.

Each directory contains a separate `README.md` with detailed instructions for reproducing the corresponding results.

## Simulation study

The simulation study compares MCNFA-CR with its Gaussian counterpart, MFA-CR, for sample sizes of 300, 500, and 1,000 and censoring levels of 0%, 10%, 20%, and 30%. Each design cell contains 80 independent Monte Carlo replications.

To reproduce Table 1 from the archived Simulation Study 1 results without rerunning the computationally intensive Monte Carlo simulation, open the `simulation1` directory in R and run:

```r
source("02_reproduce_table1.R")
```

Alternatively, from a command line opened in the `simulation1` directory, run:

```text
Rscript 02_reproduce_table1.R published_results/Study1_recovery_raw.csv
```

The complete simulation can be rerun using the scripts provided in the `simulation1` directory. Because the full simulation is computationally intensive, archived results are included to facilitate rapid reproduction and verification of the reported results.

See [`simulation1/README.md`](simulation1/README.md) for the complete workflow and additional instructions.

## Lake Michigan application

To reproduce all reported Lake Michigan tables and figures without refitting the candidate models, open the `real_data` directory in R and run:

```r
source("00_QUICK_VALIDATE.R")
```

This script reconstructs the reported tables and figures from the archived fitted-model results and performs the associated structural validation checks.

The complete analysis considers 16 candidate specifications: MFA-CR and MCNFA-CR models with one to four components and either one or two factors. Refitting all candidates and performing the targeted confirmation analysis are considerably more computationally intensive than the rapid validation workflow.

See [`real_data/README.md`](real_data/README.md) for the complete fitting, validation, and paper-to-code mapping.

The source data consist of Lake Michigan water-chemistry measurements collected in 2017 and obtained from the Biological and Chemical Oceanography Data Management Office (BCO-DMO):

<https://www.bco-dmo.org/dataset/737176>

## Requirements

A current installation of R is required. The analyses use the following R packages:

```r
install.packages(c(
  "mvtnorm",
  "MomTrunc",
  "mclust",
  "gtools",
  "ggplot2",
  "GGally",
  "rlang"
))
```

The scripts do not install packages automatically. Requirements specific to each analysis are documented in the corresponding directory-level `README.md`.

## Citation

If you use the methodology or code provided in this repository, please cite the manuscript:

> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.  
> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression.**

Users of the Lake Michigan data should also cite the original BCO-DMO dataset and its associated documentation.

## Contact

**Corresponding author:** Wan-Lun Wang  
**Repository maintainer:** F. Setoudehtazangi ([fariborz.setoudehtazangi@studenti.unipd.it](mailto:fariborz.setoudehtazangi@studenti.unipd.it))
