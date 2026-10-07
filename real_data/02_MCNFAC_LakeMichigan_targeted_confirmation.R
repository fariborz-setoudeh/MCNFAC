###############################################################
# FINAL TARGETED CONFIRMATION ANALYSIS
# Lake Michigan 2017: MFA-CR and MCNFA-CR
#
# Purpose
# -------
# This script continues only the six candidates that did not converge in the
# preceding analysis:
#   MFA-CR:   (g,q) = (3,2), (4,1), (4,2)
#   MCNFA-CR: (g,q) = (3,2), (4,1), (4,2)
#
# It does not overwrite the preceding results. It:
#   1. uses genuinely different initial partitions;
#   2. screens all starts and retains their terminal parameters;
#   3. continues three complementary starts per candidate;
#   4. diagnoses every terminal solution, even without convergence;
#   5. distinguishes formal IC values from diagnostic terminal IC values;
#   6. combines the confirmation fits with the preceding fitted candidates;
#   7. exports downstream results for the AIC-, BIC-, and EDC-selected fits;
#   8. independently checks response/covariate affine invariance by evaluating
#      the back-transformed fit on the original response and covariate scales.
###############################################################

rm(list = ls())

###############################################################
# 0. PATHS AND USER SETTINGS
###############################################################

FUNCTION_FILE <- "MCNFAC_Functions.R"
DATA_FILE <- "737176_v3_lake_michigan_chemistry.csv"
BASE_DIR <- "MCNFAC_LakeMichigan_2017_Lin_revision_v2"
BASE_RDS <- file.path(BASE_DIR, "LakeMichigan_2017_complete_analysis.rds")
OUT_DIR <- "MCNFAC_LakeMichigan_final_confirmation"

if(!file.exists(FUNCTION_FILE)) stop("Cannot find ", FUNCTION_FILE)
if(!file.exists(DATA_FILE)) stop("Cannot find ", DATA_FILE)
if(!file.exists(BASE_RDS)) stop("Cannot find ", BASE_RDS)
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(FUNCTION_FILE)

required_pkgs <- c("mvtnorm", "MomTrunc")
for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE))
    stop("Required package is not installed: ", pkg)
}

set.seed(20260920)

# Stage 1 deliberately explores several partitions. Stage 2 continues three
# complementary endpoints: the best likelihood, the best regular endpoint,
# and the best component-balanced endpoint.
N_SCREEN_STARTS <- 12L
SCREEN_ITER <- 250L
N_CONTINUE <- 3L
FINAL_TOTAL_ITER <- 3000L
TOL <- 1e-6
CONSECUTIVE_CONVERGENCE <- 3L
TRACE_EVERY <- 25L

D_FLOOR <- 1e-4
NU_LOWER <- 1e-5
NU_UPPER <- 1 - 1e-5
ETA_LOWER <- 1 + 1e-4
ETA_UPPER <- 1e3
MONOTONE_TOL <- 1e-5
MAX_CORR_CONDITION <- 1e8

# A persistent boundary is declared only when the same strong degeneracy
# pattern is observed at six consecutive 50-iteration checks after iteration
# 500. The endpoint remains in the audit table and is never treated as a
# converged maximum.
BOUNDARY_CHECK_EVERY <- 50L
BOUNDARY_CHECKS_REQUIRED <- 6L
BOUNDARY_CHECK_START <- 500L

TARGET_GRID <- expand.grid(
  model = c("MFAC", "MCNFAC"),
  g = c(3L, 4L),
  q = c(1L, 2L),
  stringsAsFactors = FALSE
)
TARGET_GRID <- TARGET_GRID[
  (TARGET_GRID$g == 3L & TARGET_GRID$q == 2L) |
    TARGET_GRID$g == 4L,
  , drop = FALSE
]

###############################################################
# 1. LOAD AND VERIFY THE PRECEDING ANALYSIS
###############################################################

base <- readRDS(BASE_RDS)
needed <- c(
  "lake", "response_names", "detection_limits", "censor_matrix",
  "l_mat", "u_mat", "X", "response_center", "response_scale",
  "observed_loglikelihood_jacobian", "X_center", "X_scale",
  "model_table", "fits"
)
missing_base <- setdiff(needed, names(base))
if(length(missing_base))
  stop("The base RDS is missing: ", paste(missing_base, collapse = ", "))

lake <- base$lake
response_names <- base$response_names
DL <- base$detection_limits
cens_raw <- base$censor_matrix
l_mat <- base$l_mat
u_mat <- base$u_mat
X <- base$X
y_center <- base$response_center
y_scale <- base$response_scale
loglik_jacobian <- base$observed_loglikelihood_jacobian
x_center <- base$X_center
x_scale <- base$X_scale
fits <- base$fits

if(nrow(lake) != 61L) stop("Expected 61 complete observations.")
if(!all(dim(l_mat) == dim(u_mat))) stop("Invalid likelihood-bound dimensions.")
if(nrow(X) != nrow(l_mat)) stop("X and response bounds have different rows.")

to_original_loglik <- function(ll_standardized){
  ll_standardized - loglik_jacobian
}

count_params_mcnfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p + 2)
}

count_params_mfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p)
}

parameter_count <- function(model, g, q){
  if(model == "MCNFAC")
    count_params_mcnfac(g, ncol(l_mat), ncol(X), q)
  else
    count_params_mfac(g, ncol(l_mat), ncol(X), q)
}

information_criteria <- function(logLik, k, n = nrow(X)){
  c(
    AIC = -2 * logLik + 2 * k,
    BIC = -2 * logLik + k * log(n),
    EDC = -2 * logLik + k * 0.2 * sqrt(n)
  )
}

###############################################################
# 2. CONSTRAINTS AND TERMINAL DIAGNOSTICS
###############################################################

project_parameters <- function(params, model){
  for(i in seq_len(params$g)){
    params$D[[i]] <- diag(pmax(diag(params$D[[i]]), D_FLOOR))
  }
  if(model == "MCNFAC"){
    params$nu <- pmin(pmax(params$nu, NU_LOWER), NU_UPPER)
    params$eta <- pmin(pmax(params$eta, ETA_LOWER), ETA_UPPER)
  }else{
    params$nu <- rep(1e-6, params$g)
    params$eta <- rep(1 + 1e-4, params$g)
  }
  rmfac_check_params(params)
  params
}

diagnose_endpoint <- function(params, estep, model, q){
  p <- params$p
  d <- params$d
  minimum_required <- max(p + 1L, d + q + 1L)
  effective_n <- colSums(estep$zhat)
  map_n <- tabulate(
    max.col(estep$zhat, ties.method = "first"),
    nbins = params$g
  )

  min_D <- Inf
  min_cov_eigen <- Inf
  max_corr_condition <- 0

  for(i in seq_len(params$g)){
    Sigma_i <- rmfac_Sigma_i(params$B[[i]], params$D[[i]])
    Sigma_i <- (Sigma_i + t(Sigma_i)) / 2
    ev_sigma <- eigen(Sigma_i, symmetric = TRUE,
                      only.values = TRUE)$values
    min_D <- min(min_D, diag(params$D[[i]]))
    min_cov_eigen <- min(min_cov_eigen, ev_sigma)

    component_sd <- sqrt(pmax(diag(Sigma_i), .Machine$double.eps))
    Corr_i <- Sigma_i / tcrossprod(component_sd)
    ev_corr <- eigen((Corr_i + t(Corr_i))/2, symmetric = TRUE,
                     only.values = TRUE)$values
    if(all(is.finite(ev_corr)) && min(ev_corr) > 0){
      max_corr_condition <- max(
        max_corr_condition,
        max(ev_corr) / min(ev_corr)
      )
    }else{
      max_corr_condition <- Inf
    }
  }

  reasons <- character(0)
  if(min(effective_n) < minimum_required)
    reasons <- c(reasons, sprintf("effective component size < %d", minimum_required))
  if(min(map_n) < minimum_required)
    reasons <- c(reasons, sprintf("MAP component size < %d", minimum_required))
  if(!is.finite(min_cov_eigen) || min_cov_eigen <= 0)
    reasons <- c(reasons, "non-positive component covariance")
  if(min_D <= D_FLOOR * (1 + 1e-6))
    reasons <- c(reasons, "uniqueness at numerical floor")
  if(!is.finite(max_corr_condition) ||
     max_corr_condition >= MAX_CORR_CONDITION)
    reasons <- c(reasons, "ill-conditioned component correlation matrix")

  if(model == "MCNFAC"){
    if(any(params$nu <= NU_LOWER * (1 + 1e-6)) ||
       any(params$nu >= NU_UPPER * (1 - 1e-6)))
      reasons <- c(reasons, "contamination proportion at numerical boundary")
    if(any(params$eta <= ETA_LOWER * (1 + 1e-6)) ||
       any(params$eta >= ETA_UPPER * (1 - 1e-6)))
      reasons <- c(reasons, "inflation factor at numerical boundary")
  }

  strong_boundary <-
    min_D <= D_FLOOR * (1 + 1e-6) &&
    (min(map_n) < minimum_required ||
       !is.finite(max_corr_condition) ||
       max_corr_condition >= MAX_CORR_CONDITION)

  list(
    admissible = length(reasons) == 0L,
    reason = if(length(reasons)) paste(unique(reasons), collapse = "; ") else "none",
    strong_boundary = strong_boundary,
    minimum_required_component_size = minimum_required,
    min_mixing_proportion = min(params$pi),
    min_effective_component_size = min(effective_n),
    min_MAP_component_size = min(map_n),
    min_uniqueness = min_D,
    min_covariance_eigenvalue = min_cov_eigen,
    max_correlation_condition_number = max_corr_condition,
    min_contamination_proportion = if(model == "MCNFAC") min(params$nu) else NA_real_,
    max_contamination_proportion = if(model == "MCNFAC") max(params$nu) else NA_real_,
    min_inflation_factor = if(model == "MCNFAC") min(params$eta) else NA_real_,
    max_inflation_factor = if(model == "MCNFAC") max(params$eta) else NA_real_
  )
}

###############################################################
# 3. GENUINELY DIVERSE INITIAL PARTITIONS
###############################################################

Y_init <- rmfac_init_impute_bounds(l_mat, u_mat)
Y_partition <- scale(Y_init)
Y_partition <- as.matrix(Y_partition)
if(any(!is.finite(Y_partition)))
  stop("Non-finite values in the partitioning matrix.")

balanced_random_partition <- function(n, g){
  sample(rep(seq_len(g), length.out = n))
}

# Some otherwise valid k-means or hierarchical starts can contain a singleton
# cluster.  Such a partition cannot initialize a within-component covariance.
# Repair only the starting partition (not the fitted model) by transferring
# observations from the largest donor cluster until every component contains
# at least two observations.
repair_partition <- function(cluster, g, Y = Y_partition, minimum_size = 2L){
  cluster <- as.integer(cluster)
  if(length(cluster) != nrow(Y) || any(!cluster %in% seq_len(g)))
    stop("Invalid proposed partition.")
  if(nrow(Y) < g * minimum_size)
    stop("The sample is too small for the requested minimum initial cluster size.")

  repeat{
    counts <- tabulate(cluster, nbins = g)
    deficient <- which(counts < minimum_size)
    if(!length(deficient)) break

    target <- deficient[1L]
    donors <- which(counts > minimum_size)
    if(!length(donors))
      stop("Unable to repair the proposed initial partition.")
    donor <- donors[which.max(counts[donors])]
    donor_rows <- which(cluster == donor)

    if(counts[target] > 0L){
      target_center <- colMeans(Y[cluster == target, , drop = FALSE])
      distances <- rowSums(
        sweep(Y[donor_rows, , drop = FALSE], 2L, target_center, "-")^2
      )
      move_row <- donor_rows[which.min(distances)]
    }else{
      donor_center <- colMeans(Y[donor_rows, , drop = FALSE])
      distances <- rowSums(
        sweep(Y[donor_rows, , drop = FALSE], 2L, donor_center, "-")^2
      )
      move_row <- donor_rows[which.max(distances)]
    }
    cluster[move_row] <- target
  }

  cluster
}

make_partition <- function(start_id, g){
  n <- nrow(Y_partition)

  if(start_id == 1L){
    cl <- kmeans(Y_partition, centers = g, nstart = 50,
                 iter.max = 100)$cluster
    cl <- repair_partition(cl, g)
    return(list(cluster = cl, type = "kmeans_50"))
  }

  if(start_id >= 2L && start_id <= 5L){
    distinct_rows <- !duplicated(as.data.frame(Y_partition))
    candidate_rows <- which(distinct_rows)
    if(length(candidate_rows) < g)
      stop("Fewer distinct initialization rows than components.")
    centers <- Y_partition[sample(candidate_rows, g), , drop = FALSE]
    cl <- kmeans(Y_partition, centers = centers, nstart = 1,
                 iter.max = 100)$cluster
    cl <- repair_partition(cl, g)
    return(list(cluster = cl, type = paste0("kmeans_random_centers_", start_id - 1L)))
  }

  if(start_id == 6L){
    cl <- cutree(hclust(dist(Y_partition), method = "ward.D2"), k = g)
    cl <- repair_partition(cl, g)
    return(list(cluster = cl, type = "hierarchical_ward"))
  }

  if(start_id == 7L){
    cl <- cutree(hclust(dist(Y_partition), method = "complete"), k = g)
    cl <- repair_partition(cl, g)
    return(list(cluster = cl, type = "hierarchical_complete"))
  }

  cl <- repair_partition(balanced_random_partition(n, g), g)
  list(cluster = cl, type = paste0("balanced_random_", start_id - 7L))
}

parameters_from_partition <- function(cluster, model, g, q, start_id){
  n <- nrow(Y_init)
  p <- ncol(Y_init)
  d <- ncol(X)

  if(length(cluster) != n || any(!cluster %in% seq_len(g)))
    stop("Invalid initial partition.")
  if(any(tabulate(cluster, nbins = g) < 2L))
    stop("An initial component contains fewer than two observations.")

  pi_start <- tabulate(cluster, nbins = g) / n
  beta <- vector("list", g)
  B <- vector("list", g)
  D <- vector("list", g)

  for(i in seq_len(g)){
    idx <- which(cluster == i)
    beta[[i]] <- rmfac_init_beta_cluster(Y_init[idx, , drop = FALSE],
                                         X[idx, , drop = FALSE])
    fa <- rmfac_init_FA(Y_init[idx, , drop = FALSE], q)
    B[[i]] <- fa$B
    D[[i]] <- diag(pmax(diag(fa$D), D_FLOOR))
  }

  contamination_grid <- data.frame(
    nu = c(0.02, 0.05, 0.10, 0.20),
    eta = c(2, 5, 10, 20)
  )
  grid_row <- 1L + ((start_id - 1L) %% nrow(contamination_grid))

  params <- list(
    g = g, p = p, q = q, d = d,
    pi = pi_start,
    beta = beta,
    B = B,
    D = D,
    nu = if(model == "MCNFAC") rep(contamination_grid$nu[grid_row], g)
      else rep(1e-6, g),
    eta = if(model == "MCNFAC") rep(contamination_grid$eta[grid_row], g)
      else rep(1 + 1e-4, g)
  )

  project_parameters(params, model)
}

###############################################################
# 4. ONE MONITORED AECM SEGMENT
###############################################################

run_segment <- function(params, model, q, max_additional_iterations,
                        initial_iteration = 0L, start_id = NA_integer_,
                        start_type = NA_character_){
  start_time <- proc.time()[["elapsed"]]
  trace_rows <- list()
  trace_index <- 1L
  convergence_streak <- 0L
  boundary_streak <- 0L
  converged <- FALSE
  persistent_boundary <- FALSE
  error_message <- ""
  minimum_increment <- Inf

  result <- tryCatch({
    params <- project_parameters(params, model)
    ll_old <- rmfac_loglikelihood(l_mat, u_mat, X, params)
    if(!is.finite(ll_old)) stop("Initial log-likelihood is not finite.")

    for(local_iter in seq_len(max_additional_iterations)){
      total_iter <- initial_iteration + local_iter

      estep <- rmfac_E_step(l_mat, u_mat, X, params)
      factor_moments <- rmfac_factor_posteriors(params, estep, X)
      estep$fhat <- factor_moments$fhat
      if(!is.null(factor_moments$Fhat)) estep$Fhat <- factor_moments$Fhat

      if(model == "MCNFAC"){
        params <- rmfac_CM1_update(
          params, estep, X,
          eta_eps = ETA_LOWER - 1
        )
        params <- project_parameters(params, model)
      }else{
        pi_new <- pmax(colMeans(estep$zhat), 1e-12)
        params$pi <- pi_new / sum(pi_new)
        params$nu <- rep(1e-6, params$g)
        params$eta <- rep(1 + 1e-4, params$g)
      }

      params_old <- params
      params <- rmfac_CM2_update_beta(params, estep, X)
      params <- rmfac_CM2_update_B(params_old, params, estep, X)
      params <- rmfac_CM2_update_D(
        params_old, params, estep, X,
        d_floor = D_FLOOR
      )
      params <- project_parameters(params, model)

      ll_new <- rmfac_loglikelihood(l_mat, u_mat, X, params)
      if(!is.finite(ll_new)) stop("Non-finite log-likelihood during iteration.")

      increment <- ll_new - ll_old
      minimum_increment <- min(minimum_increment, increment)
      if(increment < -MONOTONE_TOL)
        stop(sprintf("Log-likelihood decreased by %.6e at iteration %d.",
                     increment, total_iter))

      relative_increment <- abs(increment) / (abs(ll_old) + 1e-10)
      if(relative_increment < TOL){
        convergence_streak <- convergence_streak + 1L
      }else{
        convergence_streak <- 0L
      }
      ll_old <- ll_new

      do_trace <- total_iter == 1L ||
        total_iter %% TRACE_EVERY == 0L ||
        convergence_streak >= CONSECUTIVE_CONVERGENCE

      if(do_trace){
        endpoint_estep <- rmfac_E_step(l_mat, u_mat, X, params)
        endpoint_diag <- diagnose_endpoint(params, endpoint_estep, model, q)
        trace_rows[[trace_index]] <- data.frame(
          model = model, g = params$g, q = q,
          start_id = start_id, start_type = start_type,
          iteration = total_iter,
          logLik_standardized = ll_new,
          logLik = to_original_loglik(ll_new),
          increment = increment,
          relative_increment = relative_increment,
          min_effective_component_size = endpoint_diag$min_effective_component_size,
          min_MAP_component_size = endpoint_diag$min_MAP_component_size,
          min_uniqueness = endpoint_diag$min_uniqueness,
          max_correlation_condition_number = endpoint_diag$max_correlation_condition_number,
          diagnostic_reason = endpoint_diag$reason,
          stringsAsFactors = FALSE
        )
        trace_index <- trace_index + 1L

        # A visible heartbeat is important because the higher-order candidates
        # can require several hours.  This output also distinguishes a slow
        # active fit from a stalled R session.
        if(total_iter %% TRACE_EVERY == 0L){
          cat(
            "  ", model, "_g", params$g, "_q", q,
            " | start ", start_id,
            " | iter ", total_iter,
            " | logLik ", sprintf("%.6f", to_original_loglik(ll_new)),
            " | rel.change ", sprintf("%.3e", relative_increment),
            " | elapsed ", sprintf("%.1f min", (proc.time()[["elapsed"]] - start_time)/60),
            "\n", sep = ""
          )
          flush.console()
        }
      }

      if(total_iter >= BOUNDARY_CHECK_START &&
         total_iter %% BOUNDARY_CHECK_EVERY == 0L){
        boundary_estep <- rmfac_E_step(l_mat, u_mat, X, params)
        boundary_diag <- diagnose_endpoint(params, boundary_estep, model, q)
        if(isTRUE(boundary_diag$strong_boundary)){
          boundary_streak <- boundary_streak + 1L
        }else{
          boundary_streak <- 0L
        }
        if(boundary_streak >= BOUNDARY_CHECKS_REQUIRED){
          persistent_boundary <- TRUE
          break
        }
      }

      if(convergence_streak >= CONSECUTIVE_CONVERGENCE){
        converged <- TRUE
        break
      }
    }

    final_iteration <- initial_iteration + local_iter
    final_estep <- rmfac_E_step(l_mat, u_mat, X, params)
    final_diag <- diagnose_endpoint(params, final_estep, model, q)

    list(
      params = params,
      estep = final_estep,
      diagnostics = final_diag,
      converged = converged,
      persistent_boundary = persistent_boundary,
      iteration = final_iteration,
      logLik_standardized = ll_old,
      logLik = to_original_loglik(ll_old),
      minimum_increment = minimum_increment,
      trace = if(length(trace_rows)) do.call(rbind, trace_rows) else NULL,
      error = ""
    )
  }, error = function(e){
    error_message <<- conditionMessage(e)
    NULL
  })

  elapsed <- proc.time()[["elapsed"]] - start_time

  if(is.null(result)){
    return(list(
      params = params, estep = NULL, diagnostics = NULL,
      converged = FALSE, persistent_boundary = FALSE,
      iteration = initial_iteration,
      logLik_standardized = NA_real_, logLik = NA_real_,
      minimum_increment = if(is.finite(minimum_increment)) minimum_increment else NA_real_,
      trace = if(length(trace_rows)) do.call(rbind, trace_rows) else NULL,
      error = error_message, elapsed_seconds = elapsed,
      start_id = start_id, start_type = start_type
    ))
  }

  result$elapsed_seconds <- elapsed
  result$start_id <- start_id
  result$start_type <- start_type
  result
}

endpoint_class <- function(result){
  if(nzchar(result$error)) return("numerical_failure")
  if(isTRUE(result$converged) && isTRUE(result$diagnostics$admissible))
    return("converged_admissible")
  if(isTRUE(result$converged)) return("converged_nonregular")
  if(isTRUE(result$persistent_boundary)) return("persistent_boundary")
  if(!is.null(result$diagnostics) && !isTRUE(result$diagnostics$admissible))
    return("iteration_limit_nonregular")
  "iteration_limit_regular"
}

endpoint_row <- function(result, stage){
  diag <- result$diagnostics
  data.frame(
    stage = stage,
    start_id = result$start_id,
    start_type = result$start_type,
    terminal_class = endpoint_class(result),
    converged = isTRUE(result$converged),
    admissible = !is.null(diag) && isTRUE(diag$admissible),
    persistent_boundary = isTRUE(result$persistent_boundary),
    iteration = result$iteration,
    logLik_standardized = result$logLik_standardized,
    logLik = result$logLik,
    minimum_increment = result$minimum_increment,
    elapsed_seconds = result$elapsed_seconds,
    minimum_required_component_size = if(is.null(diag)) NA_integer_ else diag$minimum_required_component_size,
    min_mixing_proportion = if(is.null(diag)) NA_real_ else diag$min_mixing_proportion,
    min_effective_component_size = if(is.null(diag)) NA_real_ else diag$min_effective_component_size,
    min_MAP_component_size = if(is.null(diag)) NA_integer_ else diag$min_MAP_component_size,
    min_uniqueness = if(is.null(diag)) NA_real_ else diag$min_uniqueness,
    min_covariance_eigenvalue = if(is.null(diag)) NA_real_ else diag$min_covariance_eigenvalue,
    max_correlation_condition_number = if(is.null(diag)) NA_real_ else diag$max_correlation_condition_number,
    min_contamination_proportion = if(is.null(diag)) NA_real_ else diag$min_contamination_proportion,
    max_contamination_proportion = if(is.null(diag)) NA_real_ else diag$max_contamination_proportion,
    min_inflation_factor = if(is.null(diag)) NA_real_ else diag$min_inflation_factor,
    max_inflation_factor = if(is.null(diag)) NA_real_ else diag$max_inflation_factor,
    diagnostic_reason = if(is.null(diag)) NA_character_ else diag$reason,
    error = result$error,
    stringsAsFactors = FALSE
  )
}

###############################################################
# 5. SCREEN, SELECT COMPLEMENTARY ENDPOINTS, AND CONTINUE
###############################################################

choose_continuations <- function(results, number_to_keep = N_CONTINUE){
  eligible <- which(vapply(
    results,
    function(x) is.finite(x$logLik_standardized) && !nzchar(x$error),
    logical(1)
  ))
  if(!length(eligible)) return(integer(0))

  chosen <- integer(0)

  # 1. Largest attained likelihood.
  ll <- vapply(results[eligible], `[[`, numeric(1), "logLik_standardized")
  chosen <- c(chosen, eligible[which.max(ll)])

  # 2. Largest likelihood among endpoints currently satisfying diagnostics.
  regular <- eligible[vapply(
    results[eligible],
    function(x) !is.null(x$diagnostics) && isTRUE(x$diagnostics$admissible),
    logical(1)
  )]
  if(length(regular)){
    regular_ll <- vapply(results[regular], `[[`, numeric(1), "logLik_standardized")
    chosen <- c(chosen, regular[which.max(regular_ll)])
  }

  # 3. Most component-balanced endpoint, with likelihood breaking ties.
  balance <- vapply(
    results[eligible],
    function(x) if(is.null(x$diagnostics)) -Inf else
      x$diagnostics$min_effective_component_size,
    numeric(1)
  )
  best_balance <- max(balance)
  balance_candidates <- eligible[balance == best_balance]
  balance_ll <- vapply(results[balance_candidates], `[[`, numeric(1),
                       "logLik_standardized")
  chosen <- c(chosen, balance_candidates[which.max(balance_ll)])

  # Fill remaining places by attained likelihood, without duplicates.
  ordered <- eligible[order(
    vapply(results[eligible], `[[`, numeric(1), "logLik_standardized"),
    decreasing = TRUE
  )]
  chosen <- unique(c(chosen, ordered))
  head(chosen, number_to_keep)
}

checkpoint_file <- file.path(OUT_DIR, "confirmation_checkpoint.rds")
confirmation_fits <- list()
all_endpoint_rows <- list()
all_trace_rows <- list()

if(file.exists(checkpoint_file)){
  checkpoint <- try(readRDS(checkpoint_file), silent = TRUE)
  if(!inherits(checkpoint, "try-error") && is.list(checkpoint)){
    if(!is.null(checkpoint$confirmation_fits))
      confirmation_fits <- checkpoint$confirmation_fits
    if(!is.null(checkpoint$fits)) fits <- checkpoint$fits
    if(!is.null(checkpoint$endpoint_rows))
      all_endpoint_rows <- checkpoint$endpoint_rows
    if(!is.null(checkpoint$trace_rows))
      all_trace_rows <- checkpoint$trace_rows
    cat("Resuming from checkpoint with", length(confirmation_fits),
        "completed targeted candidate(s).\n")
  }else{
    warning("The confirmation checkpoint could not be read and is ignored.")
  }
}

endpoint_counter <- length(all_endpoint_rows) + 1L
trace_counter <- length(all_trace_rows) + 1L

for(candidate_index in seq_len(nrow(TARGET_GRID))){
  model <- TARGET_GRID$model[candidate_index]
  g <- TARGET_GRID$g[candidate_index]
  q <- TARGET_GRID$q[candidate_index]
  key <- paste0(model, "_g", g, "_q", q)

  if(key %in% names(confirmation_fits)){
    cat("\nUsing completed targeted candidate from checkpoint:", key, "\n")
    next
  }

  cat("\n=====================================================\n")
  cat("TARGETED CONFIRMATION:", key, "\n")
  cat("=====================================================\n")

  screen_results <- vector("list", N_SCREEN_STARTS)
  for(start_id in seq_len(N_SCREEN_STARTS)){
    set.seed(20260920 + 10000L*candidate_index + start_id)
    partition_type <- paste0("start_", start_id)
    screen_results[[start_id]] <- tryCatch({
      partition <- make_partition(start_id, g)
      partition_type <- partition$type
      params <- parameters_from_partition(
        partition$cluster, model, g, q, start_id
      )
      run_segment(
        params = params,
        model = model,
        q = q,
        max_additional_iterations = SCREEN_ITER,
        initial_iteration = 0L,
        start_id = start_id,
        start_type = partition$type
      )
    }, error = function(e){
      # A single unusable initialization is evidence about that start only; it
      # must not terminate the remaining starts or candidates.
      list(
        params = NULL, estep = NULL, diagnostics = NULL,
        converged = FALSE, persistent_boundary = FALSE,
        iteration = 0L,
        logLik_standardized = NA_real_, logLik = NA_real_,
        minimum_increment = NA_real_, trace = NULL,
        error = paste("initialization failure:", conditionMessage(e)),
        elapsed_seconds = 0,
        start_id = start_id, start_type = partition_type
      )
    })

    row <- endpoint_row(screen_results[[start_id]], "screen")
    row$model <- model; row$g <- g; row$q <- q; row$key <- key
    all_endpoint_rows[[endpoint_counter]] <- row
    endpoint_counter <- endpoint_counter + 1L

    if(!is.null(screen_results[[start_id]]$trace)){
      all_trace_rows[[trace_counter]] <- screen_results[[start_id]]$trace
      trace_counter <- trace_counter + 1L
    }

    cat(
      "Screen", start_id, "of", N_SCREEN_STARTS,
      "|", screen_results[[start_id]]$start_type,
      "| class:", endpoint_class(screen_results[[start_id]]),
      "| logLik:", round(screen_results[[start_id]]$logLik, 4), "\n"
    )
  }

  continuation_indices <- choose_continuations(screen_results)
  if(!length(continuation_indices)){
    warning("No usable screen endpoint for ", key)
    next
  }

  continuation_results <- vector("list", length(continuation_indices))
  for(h in seq_along(continuation_indices)){
    idx <- continuation_indices[h]
    screen_fit <- screen_results[[idx]]
    additional_iterations <- max(0L, FINAL_TOTAL_ITER - screen_fit$iteration)

    continuation_results[[h]] <- run_segment(
      params = screen_fit$params,
      model = model,
      q = q,
      max_additional_iterations = additional_iterations,
      initial_iteration = screen_fit$iteration,
      start_id = screen_fit$start_id,
      start_type = screen_fit$start_type
    )

    row <- endpoint_row(continuation_results[[h]], "continuation")
    row$model <- model; row$g <- g; row$q <- q; row$key <- key
    all_endpoint_rows[[endpoint_counter]] <- row
    endpoint_counter <- endpoint_counter + 1L

    if(!is.null(continuation_results[[h]]$trace)){
      all_trace_rows[[trace_counter]] <- continuation_results[[h]]$trace
      trace_counter <- trace_counter + 1L
    }

    cat(
      "Continuation", h, "of", length(continuation_indices),
      "| original start", idx,
      "| class:", endpoint_class(continuation_results[[h]]),
      "| iterations:", continuation_results[[h]]$iteration,
      "| logLik:", round(continuation_results[[h]]$logLik, 4), "\n"
    )
  }

  valid_continuations <- Filter(
    function(x) is.finite(x$logLik_standardized) && !nzchar(x$error),
    continuation_results
  )
  if(!length(valid_continuations)){
    warning("All continuation endpoints failed for ", key)
    next
  }

  regular_converged <- Filter(
    function(x) isTRUE(x$converged) && isTRUE(x$diagnostics$admissible),
    valid_continuations
  )
  selection_pool <- if(length(regular_converged))
    regular_converged else valid_continuations
  best_index <- which.max(vapply(
    selection_pool, `[[`, numeric(1), "logLik_standardized"
  ))
  best <- selection_pool[[best_index]]

  k <- parameter_count(model, g, q)
  terminal_ic <- information_criteria(best$logLik, k)
  formal <- isTRUE(best$converged) && isTRUE(best$diagnostics$admissible)

  confirmation_fits[[key]] <- list(
    model = model, g = g, q = q,
    converged = isTRUE(best$converged),
    admissible = isTRUE(best$diagnostics$admissible),
    terminal_class = endpoint_class(best),
    admissibility_reason = best$diagnostics$reason,
    params = best$params,
    estep = best$estep,
    logLik_standardized = best$logLik_standardized,
    logLik = best$logLik,
    k = k,
    AIC = if(formal) terminal_ic["AIC"] else NA_real_,
    BIC = if(formal) terminal_ic["BIC"] else NA_real_,
    EDC = if(formal) terminal_ic["EDC"] else NA_real_,
    terminal_AIC = terminal_ic["AIC"],
    terminal_BIC = terminal_ic["BIC"],
    terminal_EDC = terminal_ic["EDC"],
    iter = best$iteration,
    diagnostics = best$diagnostics,
    start_id = best$start_id,
    start_type = best$start_type
  )

  fits[[key]] <- confirmation_fits[[key]]
  saveRDS(
    list(
      confirmation_fits = confirmation_fits,
      fits = fits,
      endpoint_rows = all_endpoint_rows,
      trace_rows = all_trace_rows,
      settings = list(
        N_SCREEN_STARTS = N_SCREEN_STARTS,
        SCREEN_ITER = SCREEN_ITER,
        N_CONTINUE = N_CONTINUE,
        FINAL_TOTAL_ITER = FINAL_TOTAL_ITER,
        TOL = TOL
      )
    ),
    checkpoint_file
  )
}

endpoint_table <- do.call(rbind, all_endpoint_rows)
endpoint_table <- endpoint_table[, c(
  "key", "model", "g", "q", setdiff(names(endpoint_table),
    c("key", "model", "g", "q"))
)]
write.csv(endpoint_table,
          file.path(OUT_DIR, "01_targeted_endpoint_audit.csv"),
          row.names = FALSE)

if(length(all_trace_rows)){
  trace_table <- do.call(rbind, all_trace_rows)
  write.csv(trace_table,
            file.path(OUT_DIR, "02_iteration_trace.csv"),
            row.names = FALSE)
}

###############################################################
# 6. BUILD THE FINAL ALL-CANDIDATE TABLE
###############################################################

candidate_row <- function(key, fit){
  diag <- fit$diagnostics
  model <- fit$model
  g <- fit$g
  q <- fit$q
  k <- if(!is.null(fit$k)) fit$k else parameter_count(model, g, q)
  ll <- fit$logLik
  ll_std <- fit$logLik_standardized
  converged <- isTRUE(fit$converged)
  admissible <- isTRUE(fit$admissible)
  formal <- converged && admissible && is.finite(ll)
  terminal_ic <- if(is.finite(ll)) information_criteria(ll, k) else
    c(AIC = NA_real_, BIC = NA_real_, EDC = NA_real_)

  data.frame(
    key = key, model = model, g = g, q = q,
    converged = converged,
    admissible = admissible,
    terminal_class = if(!is.null(fit$terminal_class)) fit$terminal_class else
      if(converged && admissible) "converged_admissible" else
        if(converged) "converged_nonregular" else "not_converged",
    admissibility_reason = if(!is.null(fit$admissibility_reason))
      fit$admissibility_reason else if(!is.null(diag)) diag$reason else NA_character_,
    iteration = if(!is.null(fit$iter)) fit$iter else NA_integer_,
    logLik_standardized = ll_std,
    logLik = ll,
    k = k,
    AIC = if(formal) terminal_ic["AIC"] else NA_real_,
    BIC = if(formal) terminal_ic["BIC"] else NA_real_,
    EDC = if(formal) terminal_ic["EDC"] else NA_real_,
    terminal_AIC_diagnostic = terminal_ic["AIC"],
    terminal_BIC_diagnostic = terminal_ic["BIC"],
    terminal_EDC_diagnostic = terminal_ic["EDC"],
    minimum_required_component_size = if(is.null(diag)) NA_integer_ else diag$minimum_required_component_size,
    min_mixing_proportion = if(is.null(diag)) NA_real_ else diag$min_mixing_proportion,
    min_effective_component_size = if(is.null(diag)) NA_real_ else diag$min_effective_component_size,
    min_MAP_component_size = if(is.null(diag)) NA_integer_ else diag$min_MAP_component_size,
    min_uniqueness = if(is.null(diag)) NA_real_ else diag$min_uniqueness,
    min_covariance_eigenvalue = if(is.null(diag)) NA_real_ else diag$min_covariance_eigenvalue,
    max_correlation_condition_number = if(is.null(diag)) NA_real_ else diag$max_correlation_condition_number,
    min_contamination_proportion = if(is.null(diag)) NA_real_ else diag$min_contamination_proportion,
    max_contamination_proportion = if(is.null(diag)) NA_real_ else diag$max_contamination_proportion,
    min_inflation_factor = if(is.null(diag)) NA_real_ else diag$min_inflation_factor,
    max_inflation_factor = if(is.null(diag)) NA_real_ else diag$max_inflation_factor,
    stringsAsFactors = FALSE
  )
}

all_keys <- unlist(lapply(c("MFAC", "MCNFAC"), function(model){
  unlist(lapply(1:4, function(g){
    paste0(model, "_g", g, "_q", 1:2)
  }))
}))

final_rows <- lapply(all_keys, function(key){
  fit <- fits[[key]]
  if(is.null(fit)) stop("Missing fit object for ", key)
  candidate_row(key, fit)
})
final_table <- do.call(rbind, final_rows)

for(criterion in c("AIC", "BIC", "EDC")){
  usable <- final_table$converged & final_table$admissible &
    is.finite(final_table[[criterion]])
  final_table[[paste0("delta_", criterion)]] <- NA_real_
  final_table[[paste0("rank_", criterion)]] <- NA_integer_
  if(any(usable)){
    final_table[[paste0("delta_", criterion)]][usable] <-
      final_table[[criterion]][usable] - min(final_table[[criterion]][usable])
    final_table[[paste0("rank_", criterion)]][usable] <-
      rank(final_table[[criterion]][usable], ties.method = "min")
  }
}

final_table <- final_table[order(final_table$model, final_table$g,
                                 final_table$q), , drop = FALSE]
row.names(final_table) <- NULL
write.csv(final_table,
          file.path(OUT_DIR, "03_final_all_candidate_models.csv"),
          row.names = FALSE)

select_winner <- function(criterion, family = NULL){
  eligible <- final_table$converged & final_table$admissible &
    is.finite(final_table[[criterion]])
  if(!is.null(family)) eligible <- eligible & final_table$model == family
  tab <- final_table[eligible, , drop = FALSE]
  if(!nrow(tab)) stop("No eligible model for ", criterion,
                      if(!is.null(family)) paste0(" in ", family) else "")
  tab[which.min(tab[[criterion]]), , drop = FALSE]
}

criterion_winners_overall <- do.call(rbind, lapply(
  c("AIC", "BIC", "EDC"),
  function(criterion){
    winner <- select_winner(criterion)
    data.frame(
      criterion = criterion,
      key = winner$key,
      family = winner$model,
      g = winner$g,
      q = winner$q,
      logLik = winner$logLik,
      k = winner$k,
      criterion_value = winner[[criterion]],
      stringsAsFactors = FALSE
    )
  }
))

criterion_winners_by_family <- do.call(rbind, lapply(
  c("AIC", "BIC", "EDC"),
  function(criterion){
    do.call(rbind, lapply(c("MFAC", "MCNFAC"), function(family){
      winner <- select_winner(criterion, family)
      data.frame(
        criterion = criterion,
        family = family,
        key = winner$key,
        g = winner$g,
        q = winner$q,
        logLik = winner$logLik,
        k = winner$k,
        criterion_value = winner[[criterion]],
        stringsAsFactors = FALSE
      )
    }))
  }
))

write.csv(criterion_winners_overall,
          file.path(OUT_DIR, "04_final_criterion_winners_overall.csv"),
          row.names = FALSE)
write.csv(criterion_winners_by_family,
          file.path(OUT_DIR, "05_final_criterion_winners_by_family.csv"),
          row.names = FALSE)

###############################################################
# 7. INDEPENDENT ORIGINAL-SCALE LIKELIHOOD CHECK
###############################################################

Y_raw <- as.matrix(lake[, response_names, drop = FALSE])
l_raw <- Y_raw
u_raw <- Y_raw
for(k in seq_along(response_names)){
  idx <- which(cens_raw[, k])
  if(length(idx)){
    l_raw[idx, k] <- -Inf
    u_raw[idx, k] <- DL[k]
  }
}
X_raw <- cbind(Intercept = 1, Lat = lake$Lat, Long = lake$Long)

backtransform_fit <- function(params){
  transformed <- unserialize(serialize(params, NULL))
  response_scale_matrix <- diag(y_scale)

  for(i in seq_len(params$g)){
    beta_std <- matrix(params$beta[[i]], nrow = params$p,
                       ncol = params$d, byrow = TRUE)
    beta_raw <- beta_std
    beta_raw[, 2] <- y_scale * beta_std[, 2] / x_scale[1]
    beta_raw[, 3] <- y_scale * beta_std[, 3] / x_scale[2]
    beta_raw[, 1] <- y_center + y_scale * beta_std[, 1] -
      beta_raw[, 2] * x_center[1] -
      beta_raw[, 3] * x_center[2]
    transformed$beta[[i]] <- as.vector(t(beta_raw))
    transformed$B[[i]] <- response_scale_matrix %*% params$B[[i]]
    transformed$D[[i]] <- response_scale_matrix %*% params$D[[i]] %*%
      response_scale_matrix
  }
  rmfac_check_params(transformed)
  transformed
}

scale_check_rows <- list()
scale_check_counter <- 1L
for(key in all_keys){
  fit <- fits[[key]]
  if(is.null(fit$params) || !is.finite(fit$logLik_standardized)) next
  raw_params <- backtransform_fit(fit$params)
  direct_raw_ll <- rmfac_loglikelihood(l_raw, u_raw, X_raw, raw_params)
  transformed_ll <- to_original_loglik(fit$logLik_standardized)
  difference <- direct_raw_ll - transformed_ll
  scale_check_rows[[scale_check_counter]] <- data.frame(
    key = key,
    logLik_standardized = fit$logLik_standardized,
    transformed_original_logLik = transformed_ll,
    directly_evaluated_original_logLik = direct_raw_ll,
    difference = difference,
    passed_1e_6 = is.finite(difference) && abs(difference) <= 1e-6,
    stringsAsFactors = FALSE
  )
  scale_check_counter <- scale_check_counter + 1L
}
scale_check_table <- do.call(rbind, scale_check_rows)

# Attach the final convergence/admissibility classification.  All endpoints
# remain visible in the audit table, but only admissible converged solutions
# are relevant to the formal likelihood-equivalence check reported in the
# manuscript.
scale_check_table$converged <- final_table$converged[
  match(scale_check_table$key, final_table$key)
]
scale_check_table$admissible <- final_table$admissible[
  match(scale_check_table$key, final_table$key)
]
scale_check_table$formal_check <-
  scale_check_table$converged & scale_check_table$admissible

write.csv(scale_check_table,
          file.path(OUT_DIR, "06_affine_invariance_check.csv"),
          row.names = FALSE)

formal_scale_check <- scale_check_table$formal_check
if(!any(formal_scale_check))
  stop("No admissible converged candidate is available for the affine check.")

if(any(!scale_check_table$passed_1e_6[formal_scale_check]))
  warning("At least one admissible affine-invariance likelihood check failed.")

###############################################################
# 8. CRITERION-SPECIFIC DOWNSTREAM SUMMARIES
###############################################################

extract_observation_results <- function(fit, criterion){
  estep <- fit$estep
  if(is.null(estep)) estep <- rmfac_E_step(l_mat, u_mat, X, fit$params)
  z <- estep$zhat
  cluster <- max.col(z, ties.method = "first")

  out <- data.frame(
    RowID = lake$RowID,
    Date_UTC = lake$Date_UTC,
    Time_UTC = lake$Time_UTC,
    Site = lake$Site,
    Lat = lake$Lat,
    Long = lake$Long,
    criterion = criterion,
    model = fit$model,
    g = fit$g,
    q = fit$q,
    MAP_cluster = cluster,
    maximum_cluster_probability = apply(z, 1, max),
    stringsAsFactors = FALSE
  )
  for(i in seq_len(ncol(z))) out[[paste0("posterior_cluster_", i)]] <- z[, i]

  if(fit$model == "MCNFAC"){
    contamination_by_component <- 1 - estep$vhat
    out$contamination_probability <- rowSums(z * contamination_by_component)
    out$flag_contaminated_0.5 <- out$contamination_probability > 0.5
  }
  out
}

export_selected_fit <- function(criterion, key){
  fit <- fits[[key]]
  selected_dir <- file.path(OUT_DIR, paste0(criterion, "_", key))
  if(!dir.exists(selected_dir)) dir.create(selected_dir, recursive = TRUE)

  obs <- extract_observation_results(fit, criterion)
  write.csv(obs, file.path(selected_dir, "observation_results.csv"),
            row.names = FALSE)

  map_cluster <- obs$MAP_cluster
  component_summary <- data.frame(
    component = seq_len(fit$g),
    pi = fit$params$pi,
    nu = if(fit$model == "MCNFAC") fit$params$nu else NA_real_,
    eta = if(fit$model == "MCNFAC") fit$params$eta else NA_real_,
    MAP_cluster_size = tabulate(map_cluster, nbins = fit$g),
    stringsAsFactors = FALSE
  )
  write.csv(component_summary,
            file.path(selected_dir, "component_summary.csv"),
            row.names = FALSE)

  raw_params <- backtransform_fit(fit$params)
  beta_rows <- list()
  beta_counter <- 1L
  for(i in seq_len(fit$g)){
    beta_matrix <- matrix(raw_params$beta[[i]], nrow = fit$params$p,
                          ncol = fit$params$d, byrow = TRUE)
    colnames(beta_matrix) <- c("Intercept", "Lat", "Long")
    for(r in seq_along(response_names)){
      for(cc in seq_len(ncol(beta_matrix))){
        beta_rows[[beta_counter]] <- data.frame(
          component = i,
          response = response_names[r],
          covariate = colnames(beta_matrix)[cc],
          beta = beta_matrix[r, cc],
          stringsAsFactors = FALSE
        )
        beta_counter <- beta_counter + 1L
      }
    }
    write.csv(raw_params$B[[i]],
              file.path(selected_dir, sprintf("B_component_%d.csv", i)),
              row.names = FALSE)
    write.csv(raw_params$D[[i]],
              file.path(selected_dir, sprintf("D_component_%d.csv", i)),
              row.names = FALSE)
    write.csv(rmfac_Sigma_i(raw_params$B[[i]], raw_params$D[[i]]),
              file.path(selected_dir, sprintf("Sigma_component_%d.csv", i)),
              row.names = FALSE)
  }
  write.csv(do.call(rbind, beta_rows),
            file.path(selected_dir, "beta_original_coordinates.csv"),
            row.names = FALSE)
}

for(i in seq_len(nrow(criterion_winners_overall))){
  export_selected_fit(
    criterion_winners_overall$criterion[i],
    criterion_winners_overall$key[i]
  )
}

###############################################################
# 9. SAVE COMPLETE CONFIRMATION OBJECT AND FINAL REPORT
###############################################################

saveRDS(
  list(
    settings = list(
      N_SCREEN_STARTS = N_SCREEN_STARTS,
      SCREEN_ITER = SCREEN_ITER,
      N_CONTINUE = N_CONTINUE,
      FINAL_TOTAL_ITER = FINAL_TOTAL_ITER,
      TOL = TOL,
      D_FLOOR = D_FLOOR,
      NU_LOWER = NU_LOWER,
      NU_UPPER = NU_UPPER,
      ETA_LOWER = ETA_LOWER,
      ETA_UPPER = ETA_UPPER,
      MAX_CORR_CONDITION = MAX_CORR_CONDITION
    ),
    confirmation_fits = confirmation_fits,
    fits = fits,
    endpoint_table = endpoint_table,
    trace_table = if(exists("trace_table")) trace_table else NULL,
    final_table = final_table,
    criterion_winners_overall = criterion_winners_overall,
    criterion_winners_by_family = criterion_winners_by_family,
    scale_check_table = scale_check_table
  ),
  file.path(OUT_DIR, "LakeMichigan_final_confirmation.rds")
)

cat("\n=====================================================\n")
cat("FINAL TARGETED CONFIRMATION COMPLETED\n")
cat("=====================================================\n")
cat("\nCriterion-specific overall winners:\n")
print(criterion_winners_overall)
cat("\nTargeted candidate outcomes:\n")
print(final_table[final_table$key %in% names(confirmation_fits),
                  c("key", "terminal_class", "converged", "admissible",
                    "iteration", "logLik", "AIC", "BIC", "EDC",
                    "terminal_AIC_diagnostic", "terminal_BIC_diagnostic",
                    "terminal_EDC_diagnostic", "admissibility_reason")])
cat("\nMaximum absolute affine-invariance discrepancy among admissible fits:",
    max(abs(scale_check_table$difference[scale_check_table$formal_check]),
        na.rm = TRUE), "\n")
cat("\nResults saved in:", OUT_DIR, "\n")
cat("=====================================================\n")
