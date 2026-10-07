# Reproduce every Lake Michigan table and fitted-model figure reported in the
# manuscript and Supplementary Material from the authoritative CSV outputs.
#
# Fast archived-results run:
#   Rscript 04_LakeMichigan_publication_outputs.R published_results
#
# After a full refit, the script automatically uses:
#   MCNFAC_LakeMichigan_FINAL_RELEASE

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if(length(file_arg) == 1L){
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  setwd(dirname(script_path))
}

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^--", args)]

default_results <- if(dir.exists("MCNFAC_LakeMichigan_FINAL_RELEASE")) {
  "MCNFAC_LakeMichigan_FINAL_RELEASE"
} else {
  "published_results"
}
results_dir <- if(length(args) >= 1L) args[[1L]] else default_results
output_dir <- if(length(args) >= 2L) args[[2L]] else "publication_outputs"
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if(length(missing_packages) > 0L){
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

data_file <- "737176_v3_lake_michigan_chemistry.csv"
if(!file.exists(data_file)) stop("Missing data file: ", data_file)
if(!dir.exists(results_dir)) stop("Missing results directory: ", results_dir)

read_required <- function(...){
  path <- file.path(...)
  if(!file.exists(path)) stop("Missing required result: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

model_label <- function(key){
  family <- ifelse(grepl("^MCNFAC", key), "MCNFA-CR", "MFA-CR")
  g <- sub(".*_g([0-9]+)_q.*", "\\1", key)
  q <- sub(".*_q([0-9]+)$", "\\1", key)
  paste0(family, "(", g, ", ", q, ")")
}

response_names <- c("SRP", "TDP", "PP", "Chl", "PC", "PN")
covariate_names <- c("Lat", "Long")
analysis_names <- c(response_names, covariate_names)
detection_limits <- c(SRP = 0.50, TDP = 1.00, PP = 0.50,
                      Chl = 0.50, PC = 0.50, PN = 0.10)

raw <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
lake_2017 <- raw[raw$Year == 2017, , drop = FALSE]
lake <- lake_2017[complete.cases(lake_2017[, analysis_names]), , drop = FALSE]
row.names(lake) <- NULL
lake$RowID <- seq_len(nrow(lake))
if(nrow(lake_2017) != 73L || nrow(lake) != 61L){
  stop("Expected 73 records in 2017 and 61 complete analysis records.")
}

final_models <- read_required(results_dir, "03_final_all_candidate_models.csv")
winners <- read_required(results_dir, "04_final_criterion_winners_overall.csv")
affine <- read_required(results_dir, "06_affine_invariance_check.csv")
aic_dir <- file.path(results_dir, "AIC_MCNFAC_g2_q1")
edc_dir <- file.path(results_dir, "EDC_MCNFAC_g2_q2")
aic_obs <- read_required(aic_dir, "observation_results.csv")
edc_obs <- read_required(edc_dir, "observation_results.csv")
aic_components <- read_required(aic_dir, "component_summary.csv")
edc_components <- read_required(edc_dir, "component_summary.csv")

if(nrow(final_models) != 16L || anyDuplicated(final_models$key)){
  stop("The candidate-model audit must contain 16 unique specifications.")
}
expected_winners <- c(AIC = "MCNFAC_g2_q1", BIC = "MCNFAC_g1_q1",
                      EDC = "MCNFAC_g2_q2")
observed_winners <- setNames(winners$key, winners$criterion)
if(!identical(observed_winners[names(expected_winners)], expected_winners)){
  stop("Criterion winners differ from the verified analysis.")
}
if(nrow(aic_obs) != 61L || nrow(edc_obs) != 61L){
  stop("Selected-fit observation files must each contain 61 rows.")
}

# -----------------------------------------------------------------------------
# Supplementary Tables D.1--D.5: descriptive data summaries
# -----------------------------------------------------------------------------

missing_table <- data.frame(
  Variable = analysis_names,
  Records = nrow(lake_2017),
  Missing = vapply(lake_2017[, analysis_names, drop = FALSE],
                   function(x) sum(is.na(x)), integer(1)),
  `Missing (%)` = 100 * vapply(lake_2017[, analysis_names, drop = FALSE],
                               function(x) mean(is.na(x)), numeric(1)),
  check.names = FALSE
)
write.csv(missing_table, file.path(table_dir, "Table_D1_missingness.csv"),
          row.names = FALSE)

descriptive <- do.call(rbind, lapply(analysis_names, function(v){
  x <- lake[[v]]
  data.frame(
    Variable = v, Minimum = min(x), Q1 = unname(quantile(x, 0.25)),
    Median = median(x), Mean = mean(x), Q3 = unname(quantile(x, 0.75)),
    Maximum = max(x), SD = sd(x), check.names = FALSE
  )
}))
write.csv(descriptive, file.path(table_dir, "Table_D2_descriptive_statistics.csv"),
          row.names = FALSE)

central_moment <- function(x, order) mean((x - mean(x))^order)
distributional <- do.call(rbind, lapply(response_names, function(v){
  x <- lake[[v]]
  m2 <- central_moment(x, 2L)
  data.frame(
    Variable = v, Mean = mean(x), SD = sd(x),
    Skewness = central_moment(x, 3L) / m2^(3/2),
    Kurtosis = central_moment(x, 4L) / m2^2
  )
}))
write.csv(distributional,
          file.path(table_dir, "Table_D3_distributional_summaries.csv"),
          row.names = FALSE)

preprocessing <- read_required(results_dir, "02_preprocessing_constants.csv")
response_preprocessing <- preprocessing[
  preprocessing$variable %in% response_names,
  c("variable", "response_center", "response_scale", "detection_limit"),
  drop = FALSE
]
response_preprocessing$standardized_limit <-
  (response_preprocessing$detection_limit - response_preprocessing$response_center) /
  response_preprocessing$response_scale
write.csv(response_preprocessing,
          file.path(table_dir, "Table_D4_response_standardization.csv"),
          row.names = FALSE)

spatial_correlations <- data.frame(
  Response = response_names,
  Latitude = vapply(response_names,
                    function(v) cor(lake[[v]], lake$Lat), numeric(1)),
  Longitude = vapply(response_names,
                     function(v) cor(lake[[v]], lake$Long), numeric(1))
)
write.csv(spatial_correlations,
          file.path(table_dir, "Table_D5_response_spatial_correlations.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Supplementary Tables D.6--D.8: numerical and fitted-model diagnostics
# -----------------------------------------------------------------------------

# The confirmation-stage object contains a logical `formal_check` column,
# whereas the archived CSV intentionally retains only the likelihood audit.
# Recover the formal comparison set from the authoritative candidate table
# when that auxiliary column is absent.
if("formal_check" %in% names(affine)){
  formal_affine_rows <- !is.na(affine$formal_check) & affine$formal_check
} else {
  formal_keys <- final_models$key[
    !is.na(final_models$converged) & final_models$converged &
      !is.na(final_models$admissible) & final_models$admissible
  ]
  formal_affine_rows <- affine$key %in% formal_keys
}
affine_formal <- affine[formal_affine_rows, , drop = FALSE]
if(nrow(affine_formal) == 0L){
  stop("No admissible converged candidate is represented in the affine check.")
}
if(any(!affine_formal$passed_1e_6)){
  stop("At least one admissible candidate failed the 1e-6 affine check.")
}
affine_table <- data.frame(
  Model = model_label(affine_formal$key),
  `Transformed logLik` = affine_formal$transformed_original_logLik,
  `Direct logLik` = affine_formal$directly_evaluated_original_logLik,
  `Absolute difference` = abs(affine_formal$difference),
  check.names = FALSE
)
write.csv(affine_table,
          file.path(table_dir, "Table_D6_affine_likelihood_verification.csv"),
          row.names = FALSE)

high_component <- final_models[final_models$g %in% c(3L, 4L), , drop = FALSE]
high_component$family_order <- match(high_component$model, c("MFAC", "MCNFAC"))
high_component <- high_component[order(high_component$family_order,
                                       high_component$g, high_component$q),
                                 , drop = FALSE]
size_table <- data.frame(
  Model = model_label(high_component$key),
  `Minimum effective size` = high_component$min_effective_component_size,
  `Minimum MAP size` = high_component$min_MAP_component_size,
  check.names = FALSE
)
write.csv(size_table,
          file.path(table_dir, "Table_D7_high_component_sizes.csv"),
          row.names = FALSE)

edc_table <- data.frame(
  Component = edc_components$component,
  pi = edc_components$pi,
  nu = edc_components$nu,
  eta = edc_components$eta,
  `MAP size` = edc_components$MAP_cluster_size,
  check.names = FALSE
)
write.csv(edc_table,
          file.path(table_dir, "Table_D8_EDC_component_summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Main-paper Tables 2--5
# -----------------------------------------------------------------------------

issue_labels <- function(reason){
  labels <- character(0)
  if(grepl("uniqueness", reason, ignore.case = TRUE)) labels <- c(labels, "Uniqueness")
  if(grepl("contamination proportion", reason, ignore.case = TRUE)) {
    labels <- c(labels, "contamination")
  }
  if(grepl("inflation factor", reason, ignore.case = TRUE)) labels <- c(labels, "inflation")
  paste(labels, collapse = ", ")
}

table2 <- data.frame(
  Model = model_label(high_component$key),
  Status = ifelse(high_component$terminal_class == "converged_nonregular",
                  "NR", "PB"),
  `LogLik at termination` = high_component$logLik,
  `Boundary issue(s)` = vapply(high_component$admissibility_reason,
                               issue_labels, character(1)),
  check.names = FALSE
)
write.csv(table2, file.path(table_dir, "Table_2_convergence_diagnostics.csv"),
          row.names = FALSE)

winner_rows <- final_models[match(winners$key, final_models$key), , drop = FALSE]
table3 <- data.frame(
  Criterion = winners$criterion,
  `Selected model` = model_label(winners$key),
  logLik = winner_rows$logLik, k = winner_rows$k,
  AIC = winner_rows$AIC, BIC = winner_rows$BIC, EDC = winner_rows$EDC,
  check.names = FALSE
)
write.csv(table3, file.path(table_dir, "Table_3_criterion_winners.csv"),
          row.names = FALSE)

paired <- merge(
  aic_obs[, c("RowID", "MAP_cluster", "contamination_probability",
              "flag_contaminated_0.5")],
  edc_obs[, c("RowID", "MAP_cluster", "contamination_probability",
              "flag_contaminated_0.5")],
  by = "RowID", suffixes = c("_AIC", "_EDC"), sort = TRUE
)
identity_agreement <- mean(paired$MAP_cluster_AIC == paired$MAP_cluster_EDC)
swap_agreement <- mean(paired$MAP_cluster_AIC == (3L - paired$MAP_cluster_EDC))
if(swap_agreement > identity_agreement){
  paired$MAP_cluster_EDC_aligned <- 3L - paired$MAP_cluster_EDC
} else {
  paired$MAP_cluster_EDC_aligned <- paired$MAP_cluster_EDC
}
cross <- addmargins(table(paired$MAP_cluster_AIC,
                          paired$MAP_cluster_EDC_aligned))
table4 <- data.frame(
  AIC_component = row.names(cross),
  as.data.frame.matrix(cross),
  row.names = NULL,
  check.names = FALSE
)
write.csv(table4, file.path(table_dir, "Table_4_MAP_crossclassification.csv"),
          row.names = FALSE)

lake_aic <- merge(lake, aic_obs[, c("RowID", "MAP_cluster")],
                  by = "RowID", sort = TRUE)
table5 <- do.call(rbind, lapply(response_names, function(v){
  x1 <- lake_aic[lake_aic$MAP_cluster == 1L, v]
  x2 <- lake_aic[lake_aic$MAP_cluster == 2L, v]
  data.frame(Response = v, Component1_Mean = mean(x1), Component1_SD = sd(x1),
             Component2_Mean = mean(x2), Component2_SD = sd(x2))
}))
write.csv(table5, file.path(table_dir, "Table_5_AIC_MAP_component_profiles.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Main-paper Figures 2 and 3
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Main-paper Figures 2 and 3
# -----------------------------------------------------------------------------

effects <- data.frame(
  Response = factor(response_names, levels = response_names),
  Standardized_difference = vapply(response_names, function(v) {
    means <- tapply(lake_aic[[v]], lake_aic$MAP_cluster, mean)
    (means[["2"]] - means[["1"]]) / stats::sd(lake_aic[[v]])
  }, numeric(1))
)

utils::write.csv(
  effects,
  file.path(table_dir, "Figure_2_standardized_differences_data.csv"),
  row.names = FALSE
)

effects$Direction <- factor(
  ifelse(
    effects$Standardized_difference >= 0,
    "Positive",
    "Negative"
  ),
  levels = c("Negative", "Positive")
)

effect_colours <- c(
  "Negative" = "#6C55C8",
  "Positive" = "#079A9A"
)

p_effect <- ggplot2::ggplot(
  effects,
  ggplot2::aes(
    x = Standardized_difference,
    y = Response,
    colour = Direction
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = Standardized_difference,
      yend = Response
    ),
    linewidth = 3,
    lineend = "round",
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = 0),
    shape = 21,
    size = 3.8,
    stroke = 1,
    colour = "#9AA7B5",
    fill = "white",
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Direction),
    shape = 21,
    size = 6,
    stroke = 0.9,
    colour = "white",
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.2f", Standardized_difference),
      hjust = ifelse(Standardized_difference >= 0, -0.25, 1.25)
    ),
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "#718096"
  ) +
  ggplot2::scale_colour_manual(values = effect_colours) +
  ggplot2::scale_fill_manual(values = effect_colours) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.16, 0.18))
  ) +
  ggplot2::labs(
    title = "Between-Component Effect Sizes",
    subtitle = "Component 2 relative to Component 1",
    x = "Standardized mean difference",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      colour = "#E7ECF2",
      linewidth = 0.45
    ),
    panel.grid.major.x = ggplot2::element_line(
      colour = "#E7ECF2",
      linewidth = 0.45
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      colour = "#243247"
    ),
    axis.text.x = ggplot2::element_text(
      colour = "#718096"
    ),
    axis.title.x = ggplot2::element_text(
      face = "bold",
      colour = "#40516A",
      margin = ggplot2::margin(t = 10)
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      colour = "#1E293B",
      size = 14
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      colour = "#718096",
      size = 10
    ),
    plot.margin = ggplot2::margin(12, 30, 12, 20)
  ) +
  ggplot2::coord_cartesian(clip = "off")

ggplot2::ggsave(
  file.path(figure_dir, "Figure_2_component_mean_differences.pdf"),
  p_effect,
  width = 7.2,
  height = 4.8
)

ggplot2::ggsave(
  file.path(figure_dir, "Figure_2_component_mean_differences.png"),
  p_effect,
  width = 7.2,
  height = 4.8,
  dpi = 400,
  bg = "white"
)

correlation_from_sigma <- function(path) {
  sigma <- as.matrix(
    utils::read.csv(path, check.names = FALSE)
  )
  scale_values <- sqrt(diag(sigma))
  sigma / outer(scale_values, scale_values)
}

correlation_rows <- vector("list", 2)

for (component in 1:2) {
  corr <- correlation_from_sigma(
    file.path(
      aic_dir,
      sprintf("Sigma_component_%d.csv", component)
    )
  )
  
  dimnames(corr) <- list(response_names, response_names)
  
  grid <- expand.grid(
    Row = response_names,
    Column = response_names,
    stringsAsFactors = FALSE
  )
  
  grid$Correlation <- mapply(
    function(r, cc) corr[r, cc],
    grid$Row,
    grid$Column
  )
  
  grid$Component <- paste("Component", component)
  correlation_rows[[component]] <- grid
}

correlation_data <- do.call(rbind, correlation_rows)

correlation_data$Row <- factor(
  correlation_data$Row,
  levels = rev(response_names)
)

correlation_data$Column <- factor(
  correlation_data$Column,
  levels = response_names
)

utils::write.csv(
  correlation_data,
  file.path(table_dir, "Figure_3_fitted_correlations_data.csv"),
  row.names = FALSE
)

p_corr <- ggplot2::ggplot(
  correlation_data,
  ggplot2::aes(
    x = Column,
    y = Row,
    fill = Correlation
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.2f", Correlation)
    ),
    size = 3.1
  ) +
  ggplot2::facet_wrap(
    ~Component,
    nrow = 1
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = ggplot2::element_text(
      face = "bold"
    )
  )

ggplot2::ggsave(
  file.path(figure_dir, "Figure_3_fitted_correlations.pdf"),
  p_corr,
  width = 10.0,
  height = 4.9
)

ggplot2::ggsave(
  file.path(figure_dir, "Figure_3_fitted_correlations.png"),
  p_corr,
  width = 10.0,
  height = 4.9,
  dpi = 400,
  bg = "white"
)

# -----------------------------------------------------------------------------
# Reproducibility checks and concise report
# -----------------------------------------------------------------------------

choose2 <- function(x) x * (x - 1) / 2
adjusted_rand_index <- function(x, y){
  tab <- table(x, y); n <- sum(tab)
  observed <- sum(choose2(tab))
  rows <- sum(choose2(rowSums(tab))); cols <- sum(choose2(colSums(tab)))
  expected <- rows * cols / choose2(n)
  maximum <- 0.5 * (rows + cols)
  (observed - expected) / (maximum - expected)
}

# Guard against accidentally combining outputs from different historical fits.
# These tolerances allow ordinary floating-point variation but reject the older
# 19/42 solution that is incompatible with the final model-comparison run.
if(!identical(as.integer(aic_components$MAP_cluster_size), c(30L, 31L))){
  stop("AIC component sizes do not match the verified final fit (30, 31).")
}
if(max(abs(aic_components$pi - c(0.5049345644, 0.4950654356))) > 1e-5 ||
   max(abs(aic_components$nu - c(0.1533196377, 0.0837363177))) > 1e-5 ||
   max(abs(aic_components$eta - c(8.5225990602, 3.9707166744))) > 1e-4){
  stop("AIC component parameters do not match the verified final fit.")
}
map_matches <- sum(paired$MAP_cluster_AIC == paired$MAP_cluster_EDC_aligned)
map_ari <- adjusted_rand_index(
  paired$MAP_cluster_AIC, paired$MAP_cluster_EDC_aligned
)
if(map_matches != 57L || abs(map_ari - 0.7507616214) > 1e-8){
  stop("AIC/EDC MAP sensitivity results do not match the verified analysis.")
}

report <- c(
  "LAKE MICHIGAN PUBLICATION-OUTPUT VALIDATION",
  "===========================================",
  paste("Results source:", normalizePath(results_dir, mustWork = TRUE)),
  "2017 records: 73",
  "Complete analysis records: 61",
  paste("AIC MCNFA-CR(2,1) MAP sizes:",
        paste(aic_components$MAP_cluster_size, collapse = ", ")),
  paste("AIC MCNFA-CR(2,1) mixing proportions:",
        paste(sprintf("%.6f", aic_components$pi), collapse = ", ")),
  paste("AIC MCNFA-CR(2,1) contamination proportions:",
        paste(sprintf("%.6f", aic_components$nu), collapse = ", ")),
  paste("AIC MCNFA-CR(2,1) inflation factors:",
        paste(sprintf("%.6f", aic_components$eta), collapse = ", ")),
  paste("Mean maximum posterior probability:",
        sprintf("%.12f", mean(aic_obs$maximum_cluster_probability))),
  paste("AIC/EDC MAP matches:",
        map_matches, "/61"),
  paste("AIC/EDC ARI:", sprintf("%.12f", map_ari)),
  paste("Contamination flags (AIC, EDC, common):",
        sum(paired$flag_contaminated_0.5_AIC),
        sum(paired$flag_contaminated_0.5_EDC),
        sum(paired$flag_contaminated_0.5_AIC &
              paired$flag_contaminated_0.5_EDC)),
  paste("Maximum formal affine discrepancy:",
        format(max(abs(affine_formal$difference)), scientific = TRUE)),
  "All structural checks: PASSED"
)
writeLines(report, file.path(output_dir, "VALIDATION_REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))
cat(paste(report, collapse = "\n"), "\n")
cat("Publication outputs written to:", output_dir, "\n")
