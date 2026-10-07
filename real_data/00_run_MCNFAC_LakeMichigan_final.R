###############################################################
# COMPLETE REPRODUCIBLE LAKE MICHIGAN ANALYSIS
# MFA-CR and MCNFA-CR with left-censored responses
#
# Run this file from a clean R session.  It executes:
#   Stage 1: the full 16-model multi-start analysis;
#   Stage 2: targeted confirmation of unresolved high-component fits;
#   Stage 3: consolidation of all publication-ready final outputs;
#   Stage 4: recreation of all Lake Michigan tables and figures.
#
# Required files in the same directory:
#   MCNFAC_Functions.R
#   737176_v3_lake_michigan_chemistry.csv
#   01_MCNFAC_LakeMichigan_base_analysis.R
#   02_MCNFAC_LakeMichigan_targeted_confirmation.R
#   03_LakeMichigan_pairwise_figure.R
#   04_LakeMichigan_publication_outputs.R
#
# Main final output directory:
#   MCNFAC_LakeMichigan_FINAL_RELEASE
###############################################################

rm(list = ls())

###############################################################
# 0. LOCATE THE PROJECT AND VERIFY REQUIRED FILES
###############################################################

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)

if(length(file_arg)){
  master_file <- sub("^--file=", "", file_arg[1])
  project_dir <- dirname(normalizePath(master_file, mustWork = TRUE))
}else{
  project_dir <- normalizePath(getwd(), mustWork = TRUE)
}

setwd(project_dir)

required_files <- c(
  "MCNFAC_Functions.R",
  "737176_v3_lake_michigan_chemistry.csv",
  "01_MCNFAC_LakeMichigan_base_analysis.R",
  "02_MCNFAC_LakeMichigan_targeted_confirmation.R",
  "03_LakeMichigan_pairwise_figure.R",
  "04_LakeMichigan_publication_outputs.R"
)

missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files)){
  stop(
    "The following required file(s) are missing:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\nPlace them in: ", project_dir
  )
}

cat("\n=====================================================\n")
cat("COMPLETE LAKE MICHIGAN ANALYSIS\n")
cat("Project directory:", project_dir, "\n")
cat("=====================================================\n")

###############################################################
# 1. BASE 16-MODEL ANALYSIS
###############################################################

cat("\n\n=====================================================\n")
cat("STAGE 1 OF 3: BASE 16-MODEL ANALYSIS\n")
cat("=====================================================\n")

source("01_MCNFAC_LakeMichigan_base_analysis.R", chdir = FALSE)

base_rds <- file.path(
  "MCNFAC_LakeMichigan_2017_Lin_revision_v2",
  "LakeMichigan_2017_complete_analysis.rds"
)
if(!file.exists(base_rds)){
  stop("Stage 1 did not create the required base RDS: ", base_rds)
}

###############################################################
# 2. TARGETED CONFIRMATION ANALYSIS
###############################################################

cat("\n\n=====================================================\n")
cat("STAGE 2 OF 3: TARGETED CONFIRMATION\n")
cat("=====================================================\n")

source("02_MCNFAC_LakeMichigan_targeted_confirmation.R", chdir = FALSE)

confirmation_dir <- "MCNFAC_LakeMichigan_final_confirmation"
confirmation_rds <- file.path(
  confirmation_dir,
  "LakeMichigan_final_confirmation.rds"
)
if(!file.exists(confirmation_rds)){
  stop("Stage 2 did not create the required confirmation RDS: ",
       confirmation_rds)
}

###############################################################
# 3. VALIDATE AND CONSOLIDATE THE FINAL RESULTS
###############################################################

cat("\n\n=====================================================\n")
cat("STAGE 3 OF 3: FINAL VALIDATION AND CONSOLIDATION\n")
cat("=====================================================\n")

release_dir <- "MCNFAC_LakeMichigan_FINAL_RELEASE"
if(!dir.exists(release_dir)){
  dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
}

# Stage 2 deliberately clears the global workspace before loading the base
# analysis, so reconstruct paths needed by the consolidation stage here.
base_dir <- "MCNFAC_LakeMichigan_2017_Lin_revision_v2"
base_rds <- file.path(base_dir, "LakeMichigan_2017_complete_analysis.rds")

required_objects <- c(
  "final_table", "criterion_winners_overall",
  "criterion_winners_by_family", "scale_check_table"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if(length(missing_objects)){
  stop(
    "Stage 2 did not leave the required object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

if(nrow(final_table) != 16L){
  stop("Expected 16 candidate models, but final_table has ",
       nrow(final_table), " rows.")
}

if(anyDuplicated(final_table$key)){
  stop("Duplicate model keys were found in final_table.")
}

formal_rows <- final_table$converged & final_table$admissible
if(!any(formal_rows)){
  stop("No admissible converged models are available for formal comparison.")
}

if(any(!is.finite(final_table$AIC[formal_rows])) ||
   any(!is.finite(final_table$BIC[formal_rows])) ||
   any(!is.finite(final_table$EDC[formal_rows]))){
  stop("At least one admissible candidate has a non-finite information criterion.")
}

if(any(is.finite(final_table$AIC[!formal_rows])) ||
   any(is.finite(final_table$BIC[!formal_rows])) ||
   any(is.finite(final_table$EDC[!formal_rows]))){
  stop("Formal information criteria were retained for a nonadmissible endpoint.")
}

formal_scale_rows <- scale_check_table$formal_check
if(!any(formal_scale_rows)){
  stop("No admissible candidate is represented in the affine check.")
}
if(any(!scale_check_table$passed_1e_6[formal_scale_rows])){
  stop("An admissible candidate failed the 1e-6 affine likelihood check.")
}

# Protect the published analysis from silently changing if a different data
# file, fitting-function version, checkpoint, or numerical setting is used.
observed_winners <- setNames(
  criterion_winners_overall$key,
  criterion_winners_overall$criterion
)
expected_winners <- c(
  AIC = "MCNFAC_g2_q1",
  BIC = "MCNFAC_g1_q1",
  EDC = "MCNFAC_g2_q2"
)
if(!identical(observed_winners[names(expected_winners)], expected_winners)){
  stop(
    "The criterion winners do not match the verified manuscript results.\n",
    "Observed: ",
    paste(names(observed_winners), observed_winners, sep = "=", collapse = ", "),
    "\nExpected: ",
    paste(names(expected_winners), expected_winners, sep = "=", collapse = ", "),
    "\nDo not combine checkpoints from different code, data, or settings."
  )
}

###############################################################
# 3.1 TWO-COMPONENT SENSITIVITY ANALYSIS
###############################################################

aic_obs_file <- file.path(
  confirmation_dir,
  "AIC_MCNFAC_g2_q1",
  "observation_results.csv"
)
edc_obs_file <- file.path(
  confirmation_dir,
  "EDC_MCNFAC_g2_q2",
  "observation_results.csv"
)

if(!file.exists(aic_obs_file) || !file.exists(edc_obs_file)){
  stop("The AIC/EDC observation-level files needed for sensitivity analysis are missing.")
}

aic_obs <- read.csv(aic_obs_file, stringsAsFactors = FALSE)
edc_obs <- read.csv(edc_obs_file, stringsAsFactors = FALSE)

paired <- merge(
  aic_obs,
  edc_obs,
  by = "RowID",
  suffixes = c("_AIC", "_EDC"),
  all = FALSE,
  sort = TRUE
)

if(nrow(paired) != 61L){
  stop("Expected 61 paired observation-level results, but found ",
       nrow(paired), ".")
}

# Align the two arbitrary mixture-component labels by choosing the permutation
# with the larger MAP agreement.  For g=2, the only alternatives are the
# identity and the label swap 1 <-> 2.
agreement_identity <- mean(
  paired$MAP_cluster_AIC == paired$MAP_cluster_EDC
)
agreement_swapped <- mean(
  paired$MAP_cluster_AIC == (3L - paired$MAP_cluster_EDC)
)

if(agreement_swapped > agreement_identity){
  paired$MAP_cluster_EDC_aligned <- 3L - paired$MAP_cluster_EDC
  edc_label_alignment <- "labels_swapped"
}else{
  paired$MAP_cluster_EDC_aligned <- paired$MAP_cluster_EDC
  edc_label_alignment <- "identity"
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y){
  tab <- table(x, y)
  n <- sum(tab)
  if(n < 2L) return(NA_real_)

  index <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- 0.5 * (row_pairs + col_pairs)
  denominator <- maximum - expected

  if(abs(denominator) < .Machine$double.eps){
    return(if(identical(as.integer(x), as.integer(y))) 1 else 0)
  }
  (index - expected) / denominator
}

map_agreement_n <- sum(
  paired$MAP_cluster_AIC == paired$MAP_cluster_EDC_aligned
)
map_agreement_rate <- map_agreement_n / nrow(paired)
map_ari <- adjusted_rand_index(
  paired$MAP_cluster_AIC,
  paired$MAP_cluster_EDC_aligned
)

aic_flags <- as.logical(paired$flag_contaminated_0.5_AIC)
edc_flags <- as.logical(paired$flag_contaminated_0.5_EDC)

sensitivity_summary <- data.frame(
  quantity = c(
    "EDC label alignment",
    "Number of paired observations",
    "Matching MAP assignments",
    "MAP agreement proportion",
    "Adjusted Rand index",
    "AIC model contamination flags",
    "EDC model contamination flags",
    "Common contamination flags",
    "Pearson correlation of contamination probabilities",
    "Spearman correlation of contamination probabilities"
  ),
  value = c(
    edc_label_alignment,
    nrow(paired),
    map_agreement_n,
    sprintf("%.12f", map_agreement_rate),
    sprintf("%.12f", map_ari),
    sum(aic_flags),
    sum(edc_flags),
    sum(aic_flags & edc_flags),
    sprintf(
      "%.12f",
      cor(
        paired$contamination_probability_AIC,
        paired$contamination_probability_EDC,
        method = "pearson"
      )
    ),
    sprintf(
      "%.12f",
      cor(
        paired$contamination_probability_AIC,
        paired$contamination_probability_EDC,
        method = "spearman"
      )
    )
  ),
  stringsAsFactors = FALSE
)

map_crossclassification <- as.data.frame.matrix(
  table(
    AIC_MCNFAC_g2_q1 = paired$MAP_cluster_AIC,
    EDC_MCNFAC_g2_q2 = paired$MAP_cluster_EDC_aligned
  )
)
map_crossclassification$AIC_MCNFAC_g2_q1 <- row.names(map_crossclassification)
map_crossclassification <- map_crossclassification[
  , c("AIC_MCNFAC_g2_q1", setdiff(names(map_crossclassification),
                                   "AIC_MCNFAC_g2_q1")),
  drop = FALSE
]
row.names(map_crossclassification) <- NULL

write.csv(
  sensitivity_summary,
  file.path(release_dir, "07_two_component_sensitivity_summary.csv"),
  row.names = FALSE
)
write.csv(
  map_crossclassification,
  file.path(release_dir, "08_MAP_crossclassification_AIC_vs_EDC.csv"),
  row.names = FALSE
)

###############################################################
# 3.2 COPY THE AUTHORITATIVE OUTPUTS INTO ONE RELEASE DIRECTORY
###############################################################

copy_required <- function(source, destination = basename(source)){
  if(!file.exists(source)) stop("Required output is missing: ", source)
  target <- file.path(release_dir, destination)

  if(dir.exists(source)){
    if(!dir.exists(target)){
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
    }
    source_files <- list.files(
      source,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE,
      include.dirs = FALSE
    )
    relative_files <- substring(source_files, nchar(source) + 2L)
    target_files <- file.path(target, relative_files)
    target_dirs <- unique(dirname(target_files))
    invisible(vapply(
      target_dirs,
      dir.create,
      logical(1),
      recursive = TRUE,
      showWarnings = FALSE
    ))
    copied <- file.copy(source_files, target_files, overwrite = TRUE)
    if(length(copied) != length(source_files) || any(!copied)){
      stop("Could not copy the complete directory ", source, " to ", target)
    }
  }else{
    copied <- file.copy(source, target, overwrite = TRUE)
    if(!isTRUE(copied)) stop("Could not copy ", source, " to ", target)
  }
  invisible(target)
}

copy_required(file.path(base_dir, "01_censoring_summary.csv"))
copy_required(file.path(base_dir, "02_preprocessing_constants.csv"))
copy_required(file.path(confirmation_dir, "01_targeted_endpoint_audit.csv"),
              "03_targeted_endpoint_audit.csv")
copy_required(file.path(confirmation_dir, "02_iteration_trace.csv"),
              "04_iteration_trace.csv")
copy_required(file.path(confirmation_dir, "03_final_all_candidate_models.csv"),
              "05_final_all_candidate_models.csv")
copy_required(file.path(confirmation_dir, "04_final_criterion_winners_overall.csv"),
              "06_final_criterion_winners_overall.csv")
copy_required(file.path(confirmation_dir, "05_final_criterion_winners_by_family.csv"),
              "09_final_criterion_winners_by_family.csv")
copy_required(file.path(confirmation_dir, "06_affine_invariance_check.csv"),
              "10_affine_invariance_check.csv")

for(selected_dir in c(
  "AIC_MCNFAC_g2_q1",
  "BIC_MCNFAC_g1_q1",
  "EDC_MCNFAC_g2_q2"
)){
  copy_required(file.path(confirmation_dir, selected_dir), selected_dir)
}

copy_required(base_rds, "LakeMichigan_base_analysis.rds")
copy_required(confirmation_rds, "LakeMichigan_final_confirmation.rds")

###############################################################
# 3.3 FINAL HUMAN-READABLE SUMMARY AND SESSION INFORMATION
###############################################################

winner_lines <- apply(
  criterion_winners_overall,
  1,
  function(z){
    paste0(
      z[["criterion"]], ": ", z[["key"]],
      " (value = ", format(as.numeric(z[["criterion_value"]]), digits = 12),
      ")"
    )
  }
)

formal_max_affine_difference <- max(
  abs(scale_check_table$difference[scale_check_table$formal_check]),
  na.rm = TRUE
)

summary_lines <- c(
  "FINAL LAKE MICHIGAN MFA-CR / MCNFA-CR ANALYSIS",
  "================================================",
  "",
  paste0("Analysis date: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Complete-case sample size: 61",
  "Candidate models: 16",
  paste0("Admissible converged candidates: ", sum(formal_rows)),
  "",
  "Overall criterion winners:",
  paste0("  ", winner_lines),
  "",
  paste0(
    "Maximum admissible affine-likelihood discrepancy: ",
    format(formal_max_affine_difference, scientific = TRUE, digits = 8)
  ),
  paste0(
    "AIC/EDC two-component MAP agreement: ", map_agreement_n,
    "/", nrow(paired), " (", sprintf("%.2f", 100 * map_agreement_rate), "%)"
  ),
  paste0("AIC/EDC adjusted Rand index: ", sprintf("%.6f", map_ari)),
  paste0(
    "Contamination flags (AIC, EDC, common): ",
    sum(aic_flags), ", ", sum(edc_flags), ", ", sum(aic_flags & edc_flags)
  ),
  "",
  "Interpretive constraint:",
  "The information criteria disagree. No criterion-independent final model",
  "is declared. MCNFA-CR(2,1) is used only for AIC-conditional exploratory",
  "two-component interpretation, with sensitivity to MCNFA-CR(2,2) and",
  "explicit acknowledgement of the BIC-selected one-component model."
)

writeLines(summary_lines, file.path(release_dir, "FINAL_ANALYSIS_SUMMARY.txt"))
writeLines(
  capture.output(sessionInfo()),
  file.path(release_dir, "R_SESSION_INFO.txt")
)

manifest <- data.frame(
  file = sort(list.files(release_dir, recursive = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(release_dir, "FINAL_FILE_MANIFEST.csv"),
  row.names = FALSE
)

cat("\n=====================================================\n")
cat("FINAL RELEASE CREATED SUCCESSFULLY\n")
cat("=====================================================\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nFinal consolidated output directory:", release_dir, "\n")
cat("=====================================================\n")

###############################################################
# 4. RECREATE ALL REPORTED TABLES AND FIGURES
###############################################################

cat("\n\n=====================================================\n")
cat("STAGE 4 OF 4: PUBLICATION TABLES AND FIGURES\n")
cat("=====================================================\n")

Sys.setenv(
  LAKE_DATA_FILE = "737176_v3_lake_michigan_chemistry.csv",
  LAKE_OUTPUT_DIR = file.path("publication_outputs", "figures")
)
source("03_LakeMichigan_pairwise_figure.R", chdir = FALSE)
Sys.unsetenv(c("LAKE_DATA_FILE", "LAKE_OUTPUT_DIR"))

exit_status <- system2(
  command = file.path(R.home("bin"), "Rscript"),
  args = c(
    shQuote("04_LakeMichigan_publication_outputs.R"),
    shQuote("MCNFAC_LakeMichigan_FINAL_RELEASE"),
    shQuote("publication_outputs")
  )
)
if(!identical(exit_status, 0L)) stop("Publication-output script failed.")

if(!file.exists(file.path("publication_outputs", "VALIDATION_REPORT.txt"))){
  stop("Stage 4 did not create publication_outputs/VALIDATION_REPORT.txt")
}

cat("\nAll Lake Michigan publication outputs were created successfully.\n")
