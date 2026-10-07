# MCNFA-CR

This repository contains the data and R code used to reproduce the simulation study and the Lake Michigan application presented in the manuscript:

> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression**  
> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.

## Repository structure

- `simulation1/`: R scripts and saved results for Simulation Study 1.
- `real_data/`: data, R scripts, fitted-model results, tables, and figures for the Lake Michigan application.

Each directory contains its own `README.md` with detailed instructions for reproducing the corresponding results.

## Simulation study

The simulation study compares MCNFA-CR with its Gaussian counterpart, MFA-CR, for sample sizes \(n \in \{300, 500, 1000\}\) and censoring levels of 0%, 10%, 20%, and 30%, using 80 Monte Carlo replications per scenario.

To reproduce Table 1 from the saved Simulation Study 1 results without rerunning the computationally intensive Monte Carlo simulation, open the `simulation1` directory and run:

```r
source("02_reproduce_table1.R")
```

The complete simulation can be rerun using the scripts provided in the `simulation1` directory. Because the full simulation is computationally intensive (see the computation times reported in Table 1 of the manuscript), saved results are provided to facilitate rapid reproduction and verification of the reported results.

See `simulation1/README.md` for the complete workflow and additional instructions.

## Lake Michigan application

To reproduce the reported Lake Michigan tables and figures without refitting all candidate models, open the `real_data` directory and run:

```r
source("00_QUICK_VALIDATE.R")
```

This script reproduces the reported tables and figures from the saved fitted-model results and performs the associated structural validation checks.

Refitting the complete set of candidate models, with \(g \in \{1,2,3,4\}\) components and \(q \in \{1,2\}\) factors, is considerably more computationally intensive. See `real_data/README.md` for the complete fitting and validation workflow.

The data consist of 2017 Lake Michigan water-chemistry measurements obtained from the Biological and Chemical Oceanography Data Management Office (BCO-DMO):

<https://www.bco-dmo.org/dataset/737176>

## Requirements

- R (version 4.x recommended)
- The `MomTrunc` package for computing moments of truncated multivariate normal distributions

Any additional R packages required by individual scripts are listed in the `README.md` file of the corresponding directory.

## Citation

If you use the code or data provided in this repository, please cite the manuscript:

> F. Setoudehtazangi, T. Manouchehri, Tsung-I Lin, and Wan-Lun Wang.  
> **Mixtures of Contaminated Normal Factor Analyzers: A Robust Approach to Censored Regression.**

## Contact

**Corresponding author:** Wan-Lun Wang  
**Repository maintainer:** F. Setoudehtazangi (fariborz.setoudehtazangi@studenti.unipd.it)
