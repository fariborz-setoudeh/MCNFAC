# MCNFA-CR

This repository contains the data and R code used to reproduce the simulation study and the Lake Michigan application in the manuscript on mixtures of contaminated-normal factor analyzers for censored responses.

> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression**
> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.

## Repository structure

- `simulation1/`: R scripts and saved results for Simulation Study 1.
- `real_data/`: data, R scripts, fitted-model results, tables, and figures for the Lake Michigan application.

Each directory contains its own `README.md` with detailed instructions.

## Simulation study

To reproduce Simulation Study 1, open the `simulation1` directory and follow the instructions in:

`simulation1/README.md`

The study compares MCNFA-CR with its Gaussian counterpart MFA-CR for sample sizes n ∈ {300, 500, 1000} and censoring levels of 0%, 10%, 20% and 30%, using 80 Monte Carlo replications per scenario. Note that the full simulation is computationally heavy (see computation times in Table 1 of the manuscript), so the saved results are provided for quick checking.

## Lake Michigan application

To validate the published Lake Michigan results without refitting the models, open the `real_data` directory and run:

```r
source("00_QUICK_VALIDATE.R")
```

This reproduces the reported tables and figures from the saved fitted-model results. Refitting all candidate models (g ∈ {1, 2, 3, 4}, q ∈ {1, 2}) is much slower; see `real_data/README.md` for the full workflow.

The data are 2017 Lake Michigan water-chemistry measurements from the BCO-DMO repository: <https://www.bco-dmo.org/dataset/737176>.

## Requirements

- R (version 4.x recommended)
- The `MomTrunc` package, used to compute truncated multivariate normal moments

Any additional packages required by a given script are listed in the README of the corresponding directory.

## Citation

If you use this code or data, please cite the manuscript above.

## Contact

Corresponding author: Tsung-I Lin (tilin@nchu.edu.tw).
Repository maintainer: F. Setoudehtazangi (fariborz.setoudehtazangi@studenti.unipd.it).
