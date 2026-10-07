# ==============================================================================
# Publication-quality exploratory pairwise plot: Lake Michigan data (2017)
#
# Outputs:
#   fig/Figure_LakeMichigan_Pairwise.pdf
#   fig/Figure_LakeMichigan_Pairwise.png
#   fig/LakeMichigan_pairwise_correlations.csv
#
# Important:
#   This figure is descriptive. It uses the recorded response values for the
#   61 complete observations and does not adjust the Pearson correlations or
#   marginal densities for censoring. Dashed lines identify detection limits
#   for responses that contain censored observations.
# ==============================================================================

rm(list = ls())

# ---- 1. Packages -------------------------------------------------------------

required_packages <- c("ggplot2", "GGally", "rlang")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# ---- 2. User settings --------------------------------------------------------

# Put the CSV file in the working directory, or replace this path explicitly.
data_file <- Sys.getenv(
  "LAKE_DATA_FILE",
  unset = "737176_v3_lake_michigan_chemistry.csv"
)
output_directory <- Sys.getenv(
  "LAKE_OUTPUT_DIR",
  unset = file.path("publication_outputs", "figures")
)

response_names <- c("SRP", "TDP", "PP", "Chl", "PC", "PN")
covariate_names <- c("Lat", "Long")
analysis_names <- c(response_names, covariate_names)

# Detection limits used in the MCNFA-CR real-data analysis.
detection_limits <- c(
  SRP = 0.50,
  TDP = 1.00,
  PP  = 0.50,
  Chl = 0.50,
  PC  = 0.50,
  PN  = 0.10
)

# Expected counts provide safeguards against inadvertently plotting a different
# subset or a modified version of the dataset.
expected_2017_records <- 73L
expected_complete_records <- 61L
expected_censored_counts <- c(
  SRP = 43L,
  TDP = 16L,
  PP  = 1L,
  Chl = 10L,
  PC  = 0L,
  PN  = 0L
)

# Colour-blind-friendly colours.
point_colour <- "#1F4E79"
detection_limit_colour <- "#B2182B"

# ---- 3. Read and validate the data ------------------------------------------

if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlace it in the working directory or edit 'data_file'.",
    call. = FALSE
  )
}

lake_raw <- utils::read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_names <- c("Year", analysis_names)
absent_names <- setdiff(required_names, names(lake_raw))

if (length(absent_names) > 0L) {
  stop(
    "The following required columns are absent from the data: ",
    paste(absent_names, collapse = ", "),
    call. = FALSE
  )
}

lake_2017 <- lake_raw[lake_raw$Year == 2017, , drop = FALSE]

if (nrow(lake_2017) != expected_2017_records) {
  stop(
    "Expected ", expected_2017_records,
    " records for 2017, but found ", nrow(lake_2017), ".",
    call. = FALSE
  )
}

# Coerce only the variables used in the analysis and reject invalid entries.
for (variable in analysis_names) {
  original_values <- lake_2017[[variable]]
  numeric_values <- suppressWarnings(as.numeric(original_values))
  
  newly_missing <- is.na(numeric_values) & !is.na(original_values)
  if (any(newly_missing)) {
    stop(
      "Non-numeric values were found in column '", variable, "'.",
      call. = FALSE
    )
  }
  
  lake_2017[[variable]] <- numeric_values
}

complete_indicator <- stats::complete.cases(lake_2017[, analysis_names])
lake_complete <- lake_2017[complete_indicator, analysis_names, drop = FALSE]

if (nrow(lake_complete) != expected_complete_records) {
  stop(
    "Expected ", expected_complete_records,
    " complete records, but found ", nrow(lake_complete), ".",
    call. = FALSE
  )
}

response_data <- lake_complete[, response_names, drop = FALSE]

# Censoring is defined strictly as a recorded value below the detection limit;
# a value equal to its detection limit is not counted as censored.
censored_indicator <- sweep(
  as.matrix(response_data),
  MARGIN = 2L,
  STATS = detection_limits[response_names],
  FUN = "<"
)
censored_counts <- colSums(censored_indicator)

if (!identical(
  unname(as.integer(censored_counts)),
  unname(as.integer(expected_censored_counts[response_names]))
)) {
  stop(
    "The censoring counts do not match the expected values. Found: ",
    paste(
      paste0(names(censored_counts), "=", censored_counts),
      collapse = ", "
    ),
    call. = FALSE
  )
}

# ---- 4. Correlation matrix ---------------------------------------------------

correlation_matrix <- stats::cor(
  response_data,
  method = "pearson",
  use = "everything"
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  round(correlation_matrix, 6L),
  file = file.path(
    output_directory,
    "LakeMichigan_pairwise_correlations.csv"
  ),
  row.names = TRUE
)

# ---- 5. Custom GGally panels -------------------------------------------------

mapping_variable <- function(mapping_element) {
  rlang::as_name(rlang::get_expr(mapping_element))
}

lower_scatter_panel <- function(data, mapping, ...) {
  x_name <- mapping_variable(mapping$x)
  y_name <- mapping_variable(mapping$y)
  
  panel <- ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_point(
      shape = 21,
      size = 1.55,
      stroke = 0.25,
      colour = point_colour,
      fill = ggplot2::alpha(point_colour, 0.55),
      alpha = 0.80,
      na.rm = TRUE
    )
  
  if (censored_counts[[x_name]] > 0L) {
    panel <- panel +
      ggplot2::geom_vline(
        xintercept = detection_limits[[x_name]],
        colour = detection_limit_colour,
        linetype = "22",
        linewidth = 0.40,
        alpha = 0.80
      )
  }
  
  if (censored_counts[[y_name]] > 0L) {
    panel <- panel +
      ggplot2::geom_hline(
        yintercept = detection_limits[[y_name]],
        colour = detection_limit_colour,
        linetype = "22",
        linewidth = 0.40,
        alpha = 0.80
      )
  }
  
  panel
}

diagonal_density_panel <- function(data, mapping, ...) {
  x_name <- mapping_variable(mapping$x)
  
  panel <- ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_density(
      colour = point_colour,
      fill = point_colour,
      linewidth = 0.55,
      alpha = 0.25,
      adjust = 1,
      na.rm = TRUE
    )
  
  if (censored_counts[[x_name]] > 0L) {
    panel <- panel +
      ggplot2::geom_vline(
        xintercept = detection_limits[[x_name]],
        colour = detection_limit_colour,
        linetype = "22",
        linewidth = 0.45,
        alpha = 0.85
      )
  }
  
  panel
}

upper_correlation_panel <- function(data, mapping, ...) {
  x_values <- rlang::eval_tidy(mapping$x, data = data)
  y_values <- rlang::eval_tidy(mapping$y, data = data)
  
  correlation <- stats::cor(
    x_values,
    y_values,
    method = "pearson",
    use = "complete.obs"
  )
  
  label_size <- 3.75 + 1.20 * abs(correlation)
  
  ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::annotate(
      geom = "text",
      x = mean(range(x_values, finite = TRUE)),
      y = mean(range(y_values, finite = TRUE)),
      label = sprintf("r = %.3f", correlation),
      colour = if (correlation >= 0) point_colour else detection_limit_colour,
      family = "sans",
      fontface = "bold",
      size = label_size
    ) +
    ggplot2::theme_void()
}

# ---- 6. Construct the figure -------------------------------------------------

pairwise_plot <- GGally::ggpairs(
  data = response_data,
  columns = seq_along(response_names),
  upper = list(continuous = upper_correlation_panel),
  lower = list(continuous = lower_scatter_panel),
  diag = list(continuous = diagonal_density_panel),
  axisLabels = "show",
  columnLabels = response_names,
  labeller = "label_value",
  progress = FALSE
) +
  ggplot2::theme_bw(base_size = 9.5, base_family = "sans") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      colour = "#E6E6E6",
      linewidth = 0.25
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(
      colour = "#6E6E6E",
      linewidth = 0.40,
      fill = NA
    ),
    strip.background = ggplot2::element_rect(
      fill = "#EAF0F6",
      colour = "#6E6E6E",
      linewidth = 0.40
    ),
    strip.text = ggplot2::element_text(
      colour = "#1A1A1A",
      face = "bold",
      size = 10
    ),
    axis.title = ggplot2::element_text(size = 9.5),
    axis.text = ggplot2::element_text(
      colour = "#333333",
      size = 7.5
    ),
    plot.margin = ggplot2::margin(7, 9, 7, 7)
  )

# ---- 7. Export ---------------------------------------------------------------

pdf_file <- file.path(
  output_directory,
  "Figure_LakeMichigan_Pairwise.pdf"
)
png_file <- file.path(
  output_directory,
  "Figure_LakeMichigan_Pairwise.png"
)

ggplot2::ggsave(
  filename = pdf_file,
  plot = pairwise_plot,
  device = "pdf",
  width = 8.5,
  height = 8.5,
  units = "in",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = png_file,
  plot = pairwise_plot,
  device = "png",
  width = 8.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# ---- 8. Console report -------------------------------------------------------

cat("\nLake Michigan pairwise figure completed successfully.\n")
cat("2017 records:", nrow(lake_2017), "\n")
cat("Complete analysis records:", nrow(lake_complete), "\n")
cat(
  "Censored response counts:",
  paste(
    paste0(names(censored_counts), "=", censored_counts),
    collapse = ", "
  ),
  "\n"
)
cat("\nPearson correlation matrix (descriptive; unadjusted for censoring):\n")
print(round(correlation_matrix, 3L))
cat("\nFiles created:\n")
cat("1.", pdf_file, "\n")
cat("2.", png_file, "\n")
cat(
  "3.",
  file.path(output_directory, "LakeMichigan_pairwise_correlations.csv"),
  "\n"
)


