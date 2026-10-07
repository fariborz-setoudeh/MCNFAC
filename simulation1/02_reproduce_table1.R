# Recreate and validate the manuscript's Simulation 1 table from raw results.
# Fast validation of the archived results:
#   Rscript 02_reproduce_table1.R published_results/Study1_recovery_raw.csv

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if(length(file_arg) == 1L){
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  setwd(dirname(script_path))
}

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^--", args)]

if(!exists("RAW_RESULTS_FILE", inherits = FALSE)){
  RAW_RESULTS_FILE <- if(length(args) >= 1L) {
    args[[1L]]
  } else {
    file.path("published_results", "Study1_recovery_raw.csv")
  }
}

if(!exists("TABLE_OUTPUT_DIR", inherits = FALSE)){
  TABLE_OUTPUT_DIR <- if(length(args) >= 2L) args[[2L]] else "validation_output"
}

if(!file.exists(RAW_RESULTS_FILE)){
  stop("Raw Simulation 1 results not found: ", RAW_RESULTS_FILE)
}

dir.create(TABLE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(RAW_RESULTS_FILE, stringsAsFactors = FALSE)

required_columns <- c(
  "study", "model", "rep", "n", "censor_target", "converged",
  "beta_RMSE", "beta_abs_bias", "Sigma_RMSE", "nu_RMSE", "eta_RMSE",
  "ARI", "misclass", "cpu_time"
)
missing_columns <- setdiff(required_columns, names(raw))
if(length(missing_columns) > 0L){
  stop("Raw-results file is missing: ", paste(missing_columns, collapse = ", "))
}

raw <- raw[raw$study == "Study1_recovery", , drop = FALSE]
expected_n <- c(300, 500, 1000)
expected_censor <- c(0, 0.1, 0.2, 0.3)
expected_models <- c("MFAC", "MCNFAC")

if(nrow(raw) != 1920L) stop("Expected 1,920 Simulation 1 rows; found ", nrow(raw), ".")
if(!isTRUE(all.equal(sort(unique(raw$n)), expected_n))) {
  stop("Unexpected sample-size grid.")
}
if(!isTRUE(all.equal(sort(unique(raw$censor_target)), expected_censor))) {
  stop("Unexpected censoring grid.")
}
if(!setequal(unique(raw$model), expected_models)) stop("Unexpected model labels.")

key <- paste(raw$n, raw$censor_target, raw$model, raw$rep, sep = "|")
if(anyDuplicated(key)) stop("Duplicate simulation keys were found.")

cell <- interaction(raw$n, raw$censor_target, raw$model, drop = TRUE)
cell_sizes <- table(cell)
if(any(cell_sizes != 80L)) stop("Every design cell must contain exactly 80 replications.")

rep_ok <- vapply(
  split(raw$rep, cell),
  function(x) identical(sort(as.integer(x)), 1:80),
  logical(1)
)
if(!all(rep_ok)) stop("At least one design cell does not contain replications 1 through 80.")

finite_mean <- function(x){
  x <- x[is.finite(x)]
  if(length(x) == 0L) NA_real_ else mean(x)
}

metrics <- c(
  "beta_RMSE", "beta_abs_bias", "Sigma_RMSE", "nu_RMSE", "eta_RMSE",
  "ARI", "misclass", "converged", "cpu_time"
)

summary_raw <- aggregate(
  raw[, metrics, drop = FALSE],
  by = list(
    n = raw$n,
    censor_target = raw$censor_target,
    model = raw$model
  ),
  FUN = finite_mean
)

summary_raw$model_order <- match(summary_raw$model, expected_models)
summary_raw <- summary_raw[
  order(summary_raw$n, summary_raw$censor_target, summary_raw$model_order),
]
row.names(summary_raw) <- NULL

dash_if_na <- function(x, digits){
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

table_out <- data.frame(
  n = summary_raw$n,
  `Censoring (%)` = as.integer(round(100 * summary_raw$censor_target)),
  Model = summary_raw$model,
  RMSE_beta = formatC(summary_raw$beta_RMSE, format = "f", digits = 4),
  ABias_beta = formatC(summary_raw$beta_abs_bias, format = "f", digits = 4),
  RMSE_Sigma = formatC(summary_raw$Sigma_RMSE, format = "f", digits = 4),
  RMSE_nu = dash_if_na(summary_raw$nu_RMSE, 4),
  RMSE_eta = dash_if_na(summary_raw$eta_RMSE, 4),
  ARI = formatC(summary_raw$ARI, format = "f", digits = 4),
  MCR = formatC(summary_raw$misclass, format = "f", digits = 4),
  Conv = formatC(summary_raw$converged, format = "f", digits = 4),
  `CPU (s)` = formatC(summary_raw$cpu_time, format = "f", digits = 1),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

output_file <- file.path(TABLE_OUTPUT_DIR, "Study1_FINAL_manuscript_table_reproduced.csv")
write.csv(table_out, output_file, row.names = FALSE, quote = FALSE)

published_file <- file.path("published_results", "Study1_FINAL_manuscript_table.csv")
validation_lines <- c(
  paste("Raw file:", normalizePath(RAW_RESULTS_FILE, mustWork = TRUE)),
  paste("Rows:", nrow(raw)),
  "Design: n = 300, 500, 1000; censoring = 0%, 10%, 20%, 30%; 80 replications; two models.",
  "Structural validation: PASSED"
)

if(file.exists(published_file)){
  published <- read.csv(published_file, stringsAsFactors = FALSE, check.names = FALSE)
  if(!identical(names(published), names(table_out))) stop("Published-table columns differ.")

  normalize_cell <- function(x) sub("^-$", "--", trimws(as.character(x)))
  columns_equal <- function(a, b){
    a_chr <- normalize_cell(a)
    b_chr <- normalize_cell(b)
    a_num <- suppressWarnings(as.numeric(a_chr))
    b_num <- suppressWarnings(as.numeric(b_chr))
    numeric_positions <- is.finite(a_num) & is.finite(b_num)
    numeric_ok <- all(abs(a_num[numeric_positions] - b_num[numeric_positions]) < 1e-12)
    text_ok <- identical(a_chr[!numeric_positions], b_chr[!numeric_positions])
    numeric_ok && text_ok
  }
  time_column <- "CPU (s)"
  non_time_columns <- setdiff(names(table_out), time_column)
  same_non_time <- all(
    mapply(
      columns_equal,
      table_out[non_time_columns], published[non_time_columns]
    )
  )
  if(!same_non_time) stop("Reproduced statistical entries do not match the archived manuscript table.")

  archived_input <- identical(
    normalizePath(RAW_RESULTS_FILE, mustWork = TRUE),
    normalizePath(file.path("published_results", "Study1_recovery_raw.csv"), mustWork = TRUE)
  )
  if(archived_input){
    same_time <- columns_equal(table_out[[time_column]], published[[time_column]])
    if(!same_time) stop("Archived CPU-time entries do not match the manuscript table.")
    validation_lines <- c(validation_lines, "Published-table validation (including CPU time): PASSED")
  } else {
    validation_lines <- c(
      validation_lines,
      "Published-table validation (statistical entries, excluding hardware-dependent CPU time): PASSED"
    )
  }
} else {
  validation_lines <- c(validation_lines, "Published-table comparison: NOT RUN (archived table absent).")
}

writeLines(validation_lines, file.path(TABLE_OUTPUT_DIR, "VALIDATION_REPORT.txt"))
cat(paste(validation_lines, collapse = "\n"), "\n")
cat("Reproduced table written to:", output_file, "\n")
