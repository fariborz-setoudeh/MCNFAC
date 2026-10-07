# Fast validation of the archived Lake Michigan results.
# This script DOES NOT refit the 16 candidate models.
# Expected running time is normally seconds to a few minutes, dominated by
# rendering the high-resolution pairwise figure.

rm(list = ls())

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if(length(file_arg) == 1L){
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  setwd(dirname(script_path))
}

required_files <- c(
  "737176_v3_lake_michigan_chemistry.csv",
  "03_LakeMichigan_pairwise_figure.R",
  "04_LakeMichigan_publication_outputs.R",
  file.path("published_results", "03_final_all_candidate_models.csv"),
  file.path("published_results", "04_final_criterion_winners_overall.csv"),
  file.path("published_results", "06_affine_invariance_check.csv"),
  file.path("published_results", "AIC_MCNFAC_g2_q1", "observation_results.csv"),
  file.path("published_results", "AIC_MCNFAC_g2_q1", "component_summary.csv"),
  file.path("published_results", "EDC_MCNFAC_g2_q2", "observation_results.csv"),
  file.path("published_results", "EDC_MCNFAC_g2_q2", "component_summary.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files) > 0L){
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

required_packages <- c("ggplot2", "GGally", "rlang")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if(length(missing_packages) > 0L){
  stop(
    "Install the following package(s), then rerun: ",
    paste(missing_packages, collapse = ", ")
  )
}

cat("\n=====================================================\n")
cat("FAST LAKE MICHIGAN REPRODUCIBILITY CHECK\n")
cat("No model fitting will be performed.\n")
cat("=====================================================\n")

# Use separate environments because the pairwise script intentionally clears
# its own workspace at startup.
pairwise_environment <- new.env(parent = globalenv())
sys.source("03_LakeMichigan_pairwise_figure.R", envir = pairwise_environment)

publication_environment <- new.env(parent = globalenv())
sys.source("04_LakeMichigan_publication_outputs.R",
           envir = publication_environment)

expected_outputs <- c(
  file.path("publication_outputs", "VALIDATION_REPORT.txt"),
  file.path("publication_outputs", "R_SESSION_INFO.txt"),
  file.path("publication_outputs", "figures", "Figure_LakeMichigan_Pairwise.pdf"),
  file.path("publication_outputs", "figures", "Figure_2_component_mean_differences.pdf"),
  file.path("publication_outputs", "figures", "Figure_3_fitted_correlations.pdf"),
  file.path("publication_outputs", "tables", "Table_2_convergence_diagnostics.csv"),
  file.path("publication_outputs", "tables", "Table_3_criterion_winners.csv"),
  file.path("publication_outputs", "tables", "Table_4_MAP_crossclassification.csv"),
  file.path("publication_outputs", "tables", "Table_5_AIC_MAP_component_profiles.csv"),
  file.path("publication_outputs", "tables", "Table_D1_missingness.csv"),
  file.path("publication_outputs", "tables", "Table_D2_descriptive_statistics.csv"),
  file.path("publication_outputs", "tables", "Table_D3_distributional_summaries.csv"),
  file.path("publication_outputs", "tables", "Table_D4_response_standardization.csv"),
  file.path("publication_outputs", "tables", "Table_D5_response_spatial_correlations.csv"),
  file.path("publication_outputs", "tables", "Table_D6_affine_likelihood_verification.csv"),
  file.path("publication_outputs", "tables", "Table_D7_high_component_sizes.csv"),
  file.path("publication_outputs", "tables", "Table_D8_EDC_component_summary.csv")
)
missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if(length(missing_outputs) > 0L){
  stop("Validation did not create:\n", paste(missing_outputs, collapse = "\n"))
}

report <- readLines(file.path("publication_outputs", "VALIDATION_REPORT.txt"))
if(!any(grepl("All structural checks: PASSED", report, fixed = TRUE))){
  stop("The validation report does not contain the required PASS result.")
}

cat("\n=====================================================\n")
cat("FAST VALIDATION COMPLETED SUCCESSFULLY\n")
cat("=====================================================\n")
cat(paste(report, collapse = "\n"), "\n")
cat("\nAll tables and figures are in: publication_outputs/\n")

