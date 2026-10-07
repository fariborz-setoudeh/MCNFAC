# Master script for reproducing Simulation 1.
# Run from a terminal with: Rscript 00_run_simulation1.R

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if(length(file_arg) == 1L){
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  setwd(dirname(script_path))
}

required_files <- c(
  "01_simulation1_analysis.R",
  "02_reproduce_table1.R",
  "MCNFAC_Functions.R"
)

missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files) > 0L){
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

source("01_simulation1_analysis.R", echo = FALSE)

RAW_RESULTS_FILE <- file.path(
  "generated_results", "simulation1", "Study1_recovery_raw.csv"
)
TABLE_OUTPUT_DIR <- file.path("generated_results", "simulation1")
source("02_reproduce_table1.R", echo = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(TABLE_OUTPUT_DIR, "R_SESSION_INFO.txt")
)

cat("\nFull Simulation 1 reproduction completed.\n")

