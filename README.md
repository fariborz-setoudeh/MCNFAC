# MCNFA-CR

This repository contains the data and R code used to reproduce the simulation study and the Lake Michigan application presented in the manuscript:

> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression**  
> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.

## Repository structure

- `simulation1/`: R scripts and saved results for Simulation Study 1.
- `real_data/`: data, R scripts, fitted-model results, tables, and figures for the Lake Michigan application.

Each directory contains its own `README.md` with detailed instructions for reproducing the corresponding results.

## Simulation study

To reproduce Simulation Study 1, open the `simulation1` directory and follow the instructions in:

`simulation1/README.md`

The study compares MCNFA-CR with its Gaussian counterpart, MFA-CR, for sample sizes \(n \in \{300, 500, 1000\}\) and censoring levels of 0%, 10%, 20%, and 30%, using 80 Monte Carlo replications per scenario.

The full simulation study is computationally intensive (see the computation times reported in Table 1 of the manuscript). Saved simulation results are therefore provided to facilitate reproduction and verification of the reported results without rerunning the complete simulation.

## Lake Michigan application

To reproduce the reported Lake Michigan tables and figures without refitting all candidate models, open the `real_data` directory and run:

```r
source("00_QUICK_VALIDATE.R")
