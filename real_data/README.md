# Lake Michigan real-data reproducibility files

This directory contains the complete Lake Michigan analysis for the MFA-CR and MCNFA-CR manuscript. It includes the source data, fitting engine, full 16-model analysis, targeted confirmation of the higher-component candidates, archived fitted results, and code for every intended Lake Michigan table and figure in the manuscript and Supplementary Material.

In filenames and CSV model keys, `MFAC` and `MCNFAC` correspond to MFA-CR and MCNFA-CR, respectively.

## Files

- `00_run_MCNFAC_LakeMichigan_final.R`: master script for the complete analysis.
- `01_MCNFAC_LakeMichigan_base_analysis.R`: 16-model, multi-start candidate analysis.
- `02_MCNFAC_LakeMichigan_targeted_confirmation.R`: longer confirmation runs and admissibility diagnostics for higher-component candidates.
- `03_LakeMichigan_pairwise_figure.R`: manuscript Figure 1 and its correlation matrix.
- `04_LakeMichigan_publication_outputs.R`: manuscript Tables 2--5, Figures 2--3, and Supplementary Tables D.1--D.8.
- `MCNFAC_Functions.R`: AECM fitting functions.
- `737176_v3_lake_michigan_chemistry.csv`: source dataset downloaded from BCO-DMO.
- `published_results/`: archived CSV results used for the reported criterion winners and sensitivity calculations.

The source dataset is publicly available from the [Biological and Chemical Oceanography Data Management Office](https://www.bco-dmo.org/dataset/737176).

## Required R packages

```r
install.packages(c("mvtnorm", "MomTrunc", "ggplot2", "GGally", "rlang"))
```

The scripts do not install packages automatically.

## Fast recreation of all tables and figures

The archived fitted results allow all publication outputs to be reconstructed without repeating the expensive model fitting:

```bash
Rscript 00_QUICK_VALIDATE.R
```

From RStudio, the equivalent command is:

```r
source("00_QUICK_VALIDATE.R")
```

This does not refit any model. It verifies the archived results and recreates all Lake Michigan tables and figures. Rendering the 600-dpi pairwise figure is normally the slowest step.

Outputs are written under:

```text
publication_outputs/
├── figures/
├── tables/
├── R_SESSION_INFO.txt
└── VALIDATION_REPORT.txt
```

## Complete analysis from the source data

From a clean R session, run:

```bash
Rscript 00_run_MCNFAC_LakeMichigan_final.R
```

This runs all four stages: the candidate grid, targeted confirmation, final consolidation, and publication-output generation. The fitting stages are computationally expensive and use checkpoints to support resumption.

For a completely clean rerun, start from a fresh repository checkout. Do not combine checkpoints created with different data, fitting functions, or numerical settings.

## Mapping between the paper and code

| Paper output | Reproduction script | Output file |
|---|---|---|
| Main Figure 1 | `03_LakeMichigan_pairwise_figure.R` | `Figure_LakeMichigan_Pairwise.pdf` |
| Main Table 2 | `04_LakeMichigan_publication_outputs.R` | `Table_2_convergence_diagnostics.csv` |
| Main Table 3 | `04_LakeMichigan_publication_outputs.R` | `Table_3_criterion_winners.csv` |
| Main Table 4 | `04_LakeMichigan_publication_outputs.R` | `Table_4_MAP_crossclassification.csv` |
| Main Table 5 | `04_LakeMichigan_publication_outputs.R` | `Table_5_AIC_MAP_component_profiles.csv` |
| Main Figure 2 | `04_LakeMichigan_publication_outputs.R` | `Figure_2_component_mean_differences.pdf` |
| Main Figure 3 | `04_LakeMichigan_publication_outputs.R` | `Figure_3_fitted_correlations.pdf` |
| Supplementary Tables D.1--D.8 | `04_LakeMichigan_publication_outputs.R` | `Table_D1_...csv` through `Table_D8_...csv` |

## Interpretation rule

Information criteria are calculated only for converged, admissible solutions. Nonregular and persistent-boundary endpoints remain visible in the diagnostic files but are not treated as regular maximum-likelihood estimates. AIC, BIC, and EDC select different admissible models; therefore, the analysis does not declare a criterion-independent final model.
