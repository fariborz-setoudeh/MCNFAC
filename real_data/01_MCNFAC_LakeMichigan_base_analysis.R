###############################################################
# MCNFAC / MFAC REAL-DATA ANALYSIS
# Lake Michigan water chemistry, Year 2017
#
# Responses (Y): SRP, TDP, PP, Chl, PC, PN
# Responses: standardized for numerical stability; reported parameters and
# log-likelihoods are transformed back to the original measurement scale
# Covariates (X): intercept plus standardized Latitude and Longitude
#
# Detection limits supplied by Prof. Tsung-I Lin:
#   SRP = 0.50, TDP = 1.00, PP = 0.50,
#   Chl = 0.50, PC = 0.50, PN = 0.10
#
# IMPORTANT:
# - Observations STRICTLY BELOW the detection limit are treated
#   as left-censored. This reproduces Lin's 2017 counts:
#   SRP=43, TDP=16, PP=1, Chl=10, PC=0, PN=0.
# - The same MCNFAC/MFAC fitting engine used in the simulation
#   study is retained here.
# - Standardizing the responses is an invertible affine reparameterization.
#   Detection limits are transformed by the same map, and the likelihood is
#   adjusted by the exact Jacobian before AIC/BIC/EDC are calculated.
# - Latitude and longitude are standardized only for numerical stability;
#   an intercept is included, so this does not alter the fitted regression
#   family. Coefficients are exported in the original coordinate system.
###############################################################

rm(list = ls())

###############################################################
# 0. USER PATHS
###############################################################

FUNCTION_FILE <- "MCNFAC_Functions.R"
DATA_FILE     <- "737176_v3_lake_michigan_chemistry.csv"
OUT_DIR       <- "MCNFAC_LakeMichigan_2017_Lin_revision_v2"

if(!file.exists(FUNCTION_FILE))
  stop("Cannot find ", FUNCTION_FILE)

if(!file.exists(DATA_FILE))
  stop("Cannot find ", DATA_FILE)

if(!dir.exists(OUT_DIR))
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(FUNCTION_FILE)

required_pkgs <- c("mvtnorm", "MomTrunc")
for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE))
    stop("Required package is not installed: ", pkg)
}

# ggplot2 is used only for figures; the numerical analysis does
# not depend on it.
HAS_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)

set.seed(20260807)

###############################################################
# 1. FITTING SETTINGS
###############################################################

MAX_ITER <- 600
TOL      <- 1e-6
N_INIT   <- 10
RESUME_FROM_CHECKPOINT <- TRUE

# Scale-aware numerical safeguards. Because responses are standardized,
# D_FLOOR has the same interpretation for all six responses. A solution that
# finishes on one of these safeguards is retained in the audit trail but is
# not admitted to the information-criterion comparison.
D_FLOOR       <- 1e-4
NU_LOWER      <- 1e-5
NU_UPPER      <- 1 - 1e-5
ETA_LOWER     <- 1 + 1e-4
ETA_UPPER     <- 1e3
MONOTONE_TOL  <- 1e-5
MAX_CORR_CONDITION <- 1e8

# Main candidate grid. With n=61 and p=6, keeping q modest is
# deliberate. q=1,2 is adequate for the primary analysis.
G_GRID <- 1:4
Q_GRID <- 1:2

# BIC is used only to define provisional family representatives for
# downstream summaries.  AIC, BIC, and EDC winners are all reported, and no
# final working model is declared automatically when they disagree.
PROVISIONAL_CRITERION <- "BIC"

###############################################################
# 2. PARAMETER COUNTS AND INFORMATION CRITERIA
###############################################################

count_params_mcnfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p + 2)
}

count_params_mfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p)
}

info_criteria <- function(logLik, k, n){
  if(!is.finite(logLik)) stop("logLik must be finite.")
  if(k < 0) stop("k must be non-negative.")
  if(n <= 0) stop("n must be positive.")
  
  data.frame(
    AIC = -2 * logLik + 2 * k,
    BIC = -2 * logLik + k * log(n),
    EDC = -2 * logLik + k * 0.2 * sqrt(n)
  )
}

get_map_class <- function(estep){
  if(!is.null(estep$zhat))
    return(apply(estep$zhat, 1, which.max))
  if(!is.null(estep$z_ij))
    return(apply(estep$z_ij, 1, which.max))
  NULL
}

###############################################################
# 3. READ AND VERIFY THE 2017 DATA
###############################################################

raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

required_cols <- c(
  "Year", "Date_UTC", "Time_UTC", "Site", "Lat", "Long",
  "Depth_Site", "Depth_Smp",
  "SRP", "TDP", "PP", "Chl", "PC", "PN"
)

missing_cols <- setdiff(required_cols, names(raw))
if(length(missing_cols) > 0){
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

response_names <- c("SRP", "TDP", "PP", "Chl", "PC", "PN")
covariate_names <- c("Lat", "Long")
analysis_vars <- c(response_names, covariate_names)

lake <- raw[raw$Year == 2017, , drop = FALSE]
lake <- lake[complete.cases(lake[, analysis_vars]), , drop = FALSE]
row.names(lake) <- NULL

# Prof. Lin's example contains exactly 61 complete observations.
if(nrow(lake) != 61){
  stop(
    "The 2017 complete-case sample should contain 61 observations, but found ",
    nrow(lake), ". Check the dataset/version and filtering."
  )
}

lake$RowID <- seq_len(nrow(lake))

Y_raw <- as.matrix(lake[, response_names, drop = FALSE])
storage.mode(Y_raw) <- "double"

DL <- c(
  SRP = 0.50,
  TDP = 1.00,
  PP  = 0.50,
  Chl = 0.50,
  PC  = 0.50,
  PN  = 0.10
)

# Strict inequality is intentional. For PP, one observation is
# exactly 0.50; Lin's reported count treats it as observed, not censored.
cens_raw <- sweep(Y_raw, 2, DL, FUN = "<")

censor_counts <- colSums(cens_raw)
expected_counts <- c(SRP = 43, TDP = 16, PP = 1, Chl = 10, PC = 0, PN = 0)

cat("\n===== 2017 DATA CHECK =====\n")
cat("n =", nrow(lake), "\n")
cat("\nCensoring counts (< detection limit):\n")
print(censor_counts)

if(!identical(as.integer(censor_counts), as.integer(expected_counts))){
  stop(
    "Censoring counts do not reproduce Prof. Lin's 2017 example.\n",
    "Observed: ", paste(censor_counts, collapse = ", "), "\n",
    "Expected: ", paste(expected_counts, collapse = ", ")
  )
}

censor_summary <- data.frame(
  variable = response_names,
  detection_limit = as.numeric(DL),
  n_censored = as.integer(censor_counts),
  percent_censored = 100 * as.numeric(censor_counts) / nrow(lake),
  stringsAsFactors = FALSE
)

write.csv(
  censor_summary,
  file.path(OUT_DIR, "01_censoring_summary.csv"),
  row.names = FALSE
)

###############################################################
# 4. BUILD THE CENSORED RESPONSE MATRICES
###############################################################

# Yc_raw is the observable representation of the response:
# exact observations retain their value; left-censored values
# are represented by the corresponding detection limit.
Yc_raw <- Y_raw
for(k in seq_along(response_names)){
  Yc_raw[cens_raw[, k], k] <- DL[k]
}

# Centers and scales are computed from the observable representation, never
# from the unobserved below-limit values.
y_center <- colMeans(Yc_raw)
y_scale  <- apply(Yc_raw, 2, sd)

if(any(!is.finite(y_scale)) || any(y_scale <= 0))
  stop("Invalid response scaling constants.")

Y_scaled <- sweep(Yc_raw, 2, y_center, FUN = "-")
Y_scaled <- sweep(Y_scaled, 2, y_scale, FUN = "/")

DL_scaled <- (DL - y_center) / y_scale

# Likelihood matrices on the standardized response scale. Exact observations
# have lower=upper; a left-censored response contributes (-Inf, scaled DL].
l_mat <- Y_scaled
u_mat <- Y_scaled

for(k in seq_along(response_names)){
  idx <- which(cens_raw[, k])
  if(length(idx) > 0){
    l_mat[idx, k] <- -Inf
    u_mat[idx, k] <- DL_scaled[k]
  }
}

colnames(l_mat) <- response_names
colnames(u_mat) <- response_names

# For a partially censored observation, only exactly observed coordinates
# contribute a density Jacobian. This constant converts the standardized-scale
# observed log-likelihood to the original response scale.
loglik_jacobian <- sum(
  sweep(!cens_raw, 2, log(y_scale), FUN = "*")
)

to_original_loglik <- function(loglik_standardized){
  loglik_standardized - loglik_jacobian
}

###############################################################
# 5. COVARIATE MATRIX
###############################################################

# Latitude and longitude are centered and scaled for numerical stability.  An
# intercept is essential: without it, centering the spatial covariates would
# force every component regression mean to zero at the average location.
X_raw <- as.matrix(lake[, covariate_names, drop = FALSE])
storage.mode(X_raw) <- "double"

x_center <- colMeans(X_raw)
x_scale  <- apply(X_raw, 2, sd)

if(any(!is.finite(x_scale)) || any(x_scale <= 0))
  stop("Invalid covariate scaling constants.")

X_spatial <- scale(X_raw, center = x_center, scale = x_scale)
X_spatial <- as.matrix(X_spatial)
X <- cbind(
  Intercept = 1,
  Lat_std = X_spatial[, 1],
  Long_std = X_spatial[, 2]
)
storage.mode(X) <- "double"

# Save preprocessing constants for complete reproducibility.
preprocess_constants <- data.frame(
  variable = c(response_names, covariate_names),
  response_center = c(y_center, NA_real_, NA_real_),
  response_scale = c(y_scale, NA_real_, NA_real_),
  covariate_center = c(rep(NA_real_, length(response_names)), x_center),
  covariate_scale = c(rep(NA_real_, length(response_names)), x_scale),
  detection_limit = c(as.numeric(DL), NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
write.csv(
  preprocess_constants,
  file.path(OUT_DIR, "02_preprocessing_constants.csv"),
  row.names = FALSE
)

# Stable initialization and parameter projection for the constrained numerical
# search. Projection is used only for explicit compactness safeguards.
initialize_stable <- function(g, q){
  params <- rmfac_initialize(l_mat, u_mat, X, g = g, q = q)
  for(i in seq_len(g)){
    params$D[[i]] <- diag(pmax(diag(params$D[[i]]), D_FLOOR))
  }
  params$nu <- pmin(pmax(params$nu, NU_LOWER), NU_UPPER)
  params$eta <- pmin(pmax(params$eta, ETA_LOWER), ETA_UPPER)
  rmfac_check_params(params)
  params
}

project_mcnfac_parameters <- function(params){
  params$nu <- pmin(pmax(params$nu, NU_LOWER), NU_UPPER)
  params$eta <- pmin(pmax(params$eta, ETA_LOWER), ETA_UPPER)
  for(i in seq_len(params$g)){
    params$D[[i]] <- diag(pmax(diag(params$D[[i]]), D_FLOOR))
  }
  params
}

###############################################################
# 6. FIT MCNFAC
# Improved real-data initialization:
# MFAC-anchored contamination starts + ordinary starts
###############################################################

fit_mcnfac_once_legacy <- function(
    l_mat, u_mat, X, g, q,
    n_init = N_INIT,
    max_iter = MAX_ITER,
    tol = TOL,
    mfac_anchor = NULL
){
  
  n <- nrow(X)
  p <- ncol(l_mat)
  d <- ncol(X)
  
  best_fit <- NULL
  best_ll  <- -Inf
  all_init <- vector("list", n_init)
  
  total_start <- proc.time()[["elapsed"]]
  
  #############################################################
  # MFAC-anchored contamination starting values
  #############################################################
  
  anchor_grid <- data.frame(
    nu = c(
      1e-6,
      0.01,
      0.05,
      0.05,
      0.10,
      0.20
    ),
    eta = c(
      1 + 1e-4,
      2,
      2,
      5,
      5,
      10
    )
  )
  
  valid_anchor <-
    !is.null(mfac_anchor) &&
    !inherits(mfac_anchor, "fit_error") &&
    isTRUE(mfac_anchor$converged) &&
    !is.null(mfac_anchor$params) &&
    is.finite(mfac_anchor$logLik)
  
  if(valid_anchor){
    n_anchor <- min(nrow(anchor_grid), n_init)
  }else{
    n_anchor <- 0
  }
  
  #############################################################
  # Initializations
  #############################################################
  
  for(init in seq_len(n_init)){
    
    start_type <- "standard"
    start_nu   <- NA_real_
    start_eta  <- NA_real_
    
    # ----------------------------------------------------------
    # MFAC-anchored starts
    # ----------------------------------------------------------
    
    if(valid_anchor && init <= n_anchor){
      
      # Deep copy of fitted MFAC parameters
      params <- unserialize(
        serialize(
          mfac_anchor$params,
          NULL
        )
      )
      
      start_nu  <- anchor_grid$nu[init]
      start_eta <- anchor_grid$eta[init]
      
      params$nu  <- rep(start_nu, g)
      params$eta <- rep(start_eta, g)
      
      start_type <- "MFAC_anchor"
      
    }else{
      
      # --------------------------------------------------------
      # Ordinary simulation-style initialization
      # --------------------------------------------------------
      
      params <- initialize_stable(g = g, q = q)
    }
    
    rmfac_check_params(params)
    
    ll_old <- -Inf
    ll_new <- NA_real_
    converged <- FALSE
    
    init_start <- proc.time()[["elapsed"]]
    
    #############################################################
    # MCNFAC iteration
    #############################################################
    
    for(iter in seq_len(max_iter)){
      
      ## E-step
      
      estep <- rmfac_E_step(
        l_mat,
        u_mat,
        X,
        params
      )
      
      fact <- rmfac_factor_posteriors(
        params,
        estep,
        X
      )
      
      estep$fhat <- fact$fhat
      
      if(!is.null(fact$Fhat)){
        estep$Fhat <- fact$Fhat
      }
      
      ## Cycle 1
      
      params <- rmfac_CM1_update(
        params,
        estep,
        X
      )
      
      ## Cycle 2
      
      params_old <- params
      
      params <- rmfac_CM2_update_beta(
        params,
        estep,
        X
      )
      
      params <- rmfac_CM2_update_B(
        params_old,
        params,
        estep,
        X
      )
      
      params <- rmfac_CM2_update_D(
        params_old,
        params,
        estep,
        X
      )
      
      rmfac_check_params(params)
      
      ## Observed-data log-likelihood
      
      ll_new <- rmfac_loglikelihood(
        l_mat,
        u_mat,
        X,
        params
      )
      
      if(!is.finite(ll_new)){
        break
      }
      
      delta <- ll_new - ll_old
      
      if(is.finite(ll_old)){
        
        if(delta < -1e-7){
          warning(
            sprintf(
              paste0(
                "MCNFAC g=%d q=%d init=%d: ",
                "logLik decreased by %.6e at iter %d."
              ),
              g,
              q,
              init,
              delta,
              iter
            )
          )
        }
        
        rel_change <-
          abs(delta) /
          (abs(ll_old) + 1e-10)
        
        if(rel_change < tol){
          converged <- TRUE
          break
        }
      }
      
      ll_old <- ll_new
    }
    
    init_elapsed <-
      proc.time()[["elapsed"]] -
      init_start
    
    #############################################################
    # Save initialization information
    #############################################################
    
    all_init[[init]] <- data.frame(
      model = "MCNFAC",
      g = g,
      q = q,
      initialization = init,
      start_type = start_type,
      start_nu = start_nu,
      start_eta = start_eta,
      converged = converged,
      iter = iter,
      logLik =
        if(is.finite(ll_new))
          ll_new
      else
        NA_real_,
      elapsed_seconds = init_elapsed,
      stringsAsFactors = FALSE
    )
    
    #############################################################
    # Only converged solutions compete for best fit
    #############################################################
    
    if(
      !isTRUE(converged) ||
      !is.finite(ll_new)
    ){
      next
    }
    
    k <- count_params_mcnfac(
      g,
      p,
      d,
      q
    )
    
    IC <- info_criteria(
      ll_new,
      k,
      n
    )
    
    estep_final <- rmfac_E_step(
      l_mat,
      u_mat,
      X,
      params
    )
    
    zhat <- get_map_class(
      estep_final
    )
    
    entropy <- -mean(
      rowSums(
        estep_final$zhat *
          log(estep_final$zhat + 1e-12)
      )
    )
    
    current_fit <- list(
      model = "MCNFAC",
      g = g,
      q = q,
      logLik = ll_new,
      AIC = IC$AIC,
      BIC = IC$BIC,
      EDC = IC$EDC,
      k = k,
      entropy = entropy,
      iter = iter,
      converged = converged,
      cpu_time = init_elapsed,
      params = params,
      zhat = zhat,
      estep = estep_final,
      initialization = init,
      start_type = start_type,
      start_nu = start_nu,
      start_eta = start_eta
    )
    
    if(ll_new > best_ll){
      best_ll <- ll_new
      best_fit <- current_fit
    }
  }
  
  total_elapsed <-
    proc.time()[["elapsed"]] -
    total_start
  
  init_table <- do.call(
    rbind,
    all_init
  )
  
  #############################################################
  # If no initialization converged
  #############################################################
  
  if(is.null(best_fit)){
    
    return(
      list(
        model = "MCNFAC",
        g = g,
        q = q,
        converged = FALSE,
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        EDC = NA_real_,
        k = count_params_mcnfac(g, p, d, q),
        entropy = NA_real_,
        iter = NA_integer_,
        cpu_time = NA_real_,
        total_elapsed = total_elapsed,
        initialization = NA_integer_,
        params = NULL,
        zhat = NULL,
        estep = NULL,
        init_table = init_table
      )
    )
  }
  
  best_fit$total_elapsed <- total_elapsed
  best_fit$init_table <- init_table
  
  best_fit
}

###############################################################
# 7. FIT MFAC: SAME NESTED-GAUSSIAN IMPLEMENTATION AS SIMULATION
###############################################################

fit_mfac_once_legacy <- function(
    l_mat, u_mat, X, g, q,
    n_init = N_INIT,
    max_iter = MAX_ITER,
    tol = TOL
){
  n <- nrow(X)
  p <- ncol(l_mat)
  d <- ncol(X)
  
  best_fit <- NULL
  best_ll  <- -Inf
  all_init <- vector("list", n_init)
  
  total_start <- proc.time()[["elapsed"]]
  
  for(init in seq_len(n_init)){
    params <- initialize_stable(g = g, q = q)
    
    fixed_nu  <- rep(1e-6, g)
    fixed_eta <- rep(1 + 1e-4, g)
    params$nu  <- fixed_nu
    params$eta <- fixed_eta
    
    ll_old <- -Inf
    ll_new <- NA_real_
    converged <- FALSE
    
    init_start <- proc.time()[["elapsed"]]
    
    for(iter in seq_len(max_iter)){
      estep <- rmfac_E_step(l_mat, u_mat, X, params)
      
      fact <- rmfac_factor_posteriors(params, estep, X)
      estep$fhat <- fact$fhat
      if(!is.null(fact$Fhat)) estep$Fhat <- fact$Fhat
      
      pi_new <- colMeans(estep$zhat)
      pi_new <- pmax(pi_new, 1e-12)
      params$pi <- pi_new / sum(pi_new)
      
      params$nu  <- fixed_nu
      params$eta <- fixed_eta
      
      params_old <- params
      params <- rmfac_CM2_update_beta(params, estep, X)
      params <- rmfac_CM2_update_B(params_old, params, estep, X)
      params <- rmfac_CM2_update_D(params_old, params, estep, X)
      
      params$nu  <- fixed_nu
      params$eta <- fixed_eta
      
      rmfac_check_params(params)
      
      ll_new <- rmfac_loglikelihood(l_mat, u_mat, X, params)
      if(!is.finite(ll_new)) break
      
      delta <- ll_new - ll_old
      
      if(is.finite(ll_old)){
        if(delta < -1e-7){
          warning(sprintf(
            "MFAC g=%d q=%d init=%d: logLik decreased by %.6e at iter %d.",
            g, q, init, delta, iter
          ))
        }
        
        rel_change <- abs(delta) / (abs(ll_old) + 1e-10)
        if(rel_change < tol){
          converged <- TRUE
          break
        }
      }
      
      ll_old <- ll_new
    }
    
    init_elapsed <- proc.time()[["elapsed"]] - init_start
    
    all_init[[init]] <- data.frame(
      model = "MFAC",
      g = g,
      q = q,
      initialization = init,
      start_type = "standard",
      start_nu = NA_real_,
      start_eta = NA_real_,
      converged = converged,
      iter = iter,
      logLik =
        if(is.finite(ll_new))
          ll_new
      else
        NA_real_,
      elapsed_seconds = init_elapsed,
      stringsAsFactors = FALSE
    )
    
    if(!isTRUE(converged) || !is.finite(ll_new)) next
    
    k <- count_params_mfac(g, p, d, q)
    IC <- info_criteria(ll_new, k, n)
    
    estep_final <- rmfac_E_step(l_mat, u_mat, X, params)
    zhat <- get_map_class(estep_final)
    
    entropy <- -mean(
      rowSums(estep_final$zhat * log(estep_final$zhat + 1e-12))
    )
    
    current_fit <- list(
      model = "MFAC",
      g = g, q = q,
      logLik = ll_new,
      AIC = IC$AIC, BIC = IC$BIC, EDC = IC$EDC,
      k = k,
      entropy = entropy,
      iter = iter,
      converged = converged,
      cpu_time = init_elapsed,
      params = params,
      zhat = zhat,
      estep = estep_final,
      initialization = init
    )
    
    if(ll_new > best_ll){
      best_ll <- ll_new
      best_fit <- current_fit
    }
  }
  
  total_elapsed <- proc.time()[["elapsed"]] - total_start
  init_table <- do.call(rbind, all_init)
  
  if(is.null(best_fit)){
    return(list(
      model = "MFAC", g = g, q = q,
      converged = FALSE,
      logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, EDC = NA_real_,
      k = count_params_mfac(g, p, d, q),
      entropy = NA_real_, iter = NA_integer_,
      cpu_time = NA_real_, total_elapsed = total_elapsed,
      initialization = NA_integer_, params = NULL,
      zhat = NULL, estep = NULL, init_table = init_table
    ))
  }
  
  best_fit$total_elapsed <- total_elapsed
  best_fit$init_table <- init_table
  best_fit
}

fit_model_once <- function(
    model,
    l_mat,
    u_mat,
    X,
    g,
    q,
    n_init = N_INIT,
    mfac_anchor = NULL
){
  
  if(model == "MCNFAC"){
    
    return(
      fit_mcnfac_once(
        l_mat = l_mat,
        u_mat = u_mat,
        X = X,
        g = g,
        q = q,
        n_init = n_init,
        mfac_anchor = mfac_anchor
      )
    )
  }
  
  if(model == "MFAC"){
    
    return(
      fit_mfac_once(
        l_mat = l_mat,
        u_mat = u_mat,
        X = X,
        g = g,
        q = q,
        n_init = n_init
      )
    )
  }
  
  stop("Unknown model: ", model)
}

###############################################################
# 7A. ROBUST REAL-DATA FITTING WRAPPER
#
# These definitions replace the legacy wrappers above. The underlying E-step
# and CM updates in MCNFAC_Functions.R are unchanged. The wrapper adds:
#   (i) error isolation for every initialization;
#   (ii) standardized-scale compactness safeguards;
#   (iii) monotonicity checks;
#   (iv) predeclared admissibility diagnostics; and
#   (v) selection of the best admissible, not merely largest, likelihood.
###############################################################

assess_solution <- function(params, estep, model, q){
  n <- nrow(estep$zhat)
  p <- params$p
  d <- params$d
  min_component_n <- max(p + 1L, d + q + 1L)
  
  effective_n <- colSums(estep$zhat)
  map_n <- tabulate(max.col(estep$zhat, ties.method = "first"),
                    nbins = params$g)
  
  min_D <- Inf
  min_cov_eigen <- Inf
  max_corr_condition <- 0
  
  for(i in seq_len(params$g)){
    D_i <- params$D[[i]]
    Sigma_i <- rmfac_Sigma_i(params$B[[i]], D_i)
    Sigma_i <- (Sigma_i + t(Sigma_i)) / 2
    ev_sigma <- eigen(Sigma_i, symmetric = TRUE, only.values = TRUE)$values
    min_D <- min(min_D, diag(D_i))
    min_cov_eigen <- min(min_cov_eigen, ev_sigma)
    
    s <- sqrt(pmax(diag(Sigma_i), .Machine$double.eps))
    Corr_i <- Sigma_i / tcrossprod(s)
    ev_corr <- eigen((Corr_i + t(Corr_i))/2,
                     symmetric = TRUE, only.values = TRUE)$values
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
  if(min(effective_n) < min_component_n)
    reasons <- c(reasons, sprintf("effective component size < %d", min_component_n))
  if(min(map_n) < min_component_n)
    reasons <- c(reasons, sprintf("MAP component size < %d", min_component_n))
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
  
  list(
    admissible = length(reasons) == 0L,
    reason = if(length(reasons)) paste(unique(reasons), collapse = "; ") else "none",
    minimum_required_component_size = min_component_n,
    min_effective_component_size = min(effective_n),
    min_MAP_component_size = min(map_n),
    min_uniqueness = min_D,
    min_covariance_eigenvalue = min_cov_eigen,
    max_correlation_condition_number = max_corr_condition,
    min_mixing_proportion = min(params$pi),
    min_contamination_proportion = if(model == "MCNFAC") min(params$nu) else NA_real_,
    max_contamination_proportion = if(model == "MCNFAC") max(params$nu) else NA_real_,
    min_inflation_factor = if(model == "MCNFAC") min(params$eta) else NA_real_,
    max_inflation_factor = if(model == "MCNFAC") max(params$eta) else NA_real_
  )
}

make_start_record <- function(model, g, q, init, start_type,
                              start_nu = NA_real_, start_eta = NA_real_,
                              converged = FALSE, admissible = FALSE,
                              iter = NA_integer_, ll_std = NA_real_,
                              ll_original = NA_real_, min_delta = NA_real_,
                              elapsed = NA_real_, reason = NA_character_,
                              error = "", diagnostics = NULL){
  data.frame(
    model = model,
    g = g,
    q = q,
    initialization = init,
    start_type = start_type,
    start_nu = start_nu,
    start_eta = start_eta,
    converged = converged,
    admissible = admissible,
    iter = iter,
    logLik_standardized = ll_std,
    logLik = ll_original,
    min_logLik_increment = min_delta,
    elapsed_seconds = elapsed,
    min_effective_component_size = if(is.null(diagnostics)) NA_real_ else diagnostics$min_effective_component_size,
    min_MAP_component_size = if(is.null(diagnostics)) NA_real_ else diagnostics$min_MAP_component_size,
    min_uniqueness = if(is.null(diagnostics)) NA_real_ else diagnostics$min_uniqueness,
    max_correlation_condition_number = if(is.null(diagnostics)) NA_real_ else diagnostics$max_correlation_condition_number,
    admissibility_reason = reason,
    error = error,
    stringsAsFactors = FALSE
  )
}

fit_single_start <- function(model, l_mat, u_mat, X, g, q, init,
                             start_type = "standard", start_nu = NA_real_,
                             start_eta = NA_real_, anchor = NULL,
                             max_iter = MAX_ITER, tol = TOL){
  start_time <- proc.time()[["elapsed"]]
  
  tryCatch({
    if(start_type == "MFAC_anchor"){
      params <- unserialize(serialize(anchor$params, NULL))
      params$nu <- rep(start_nu, g)
      params$eta <- rep(start_eta, g)
      params <- project_mcnfac_parameters(params)
    }else{
      params <- initialize_stable(g, q)
    }
    
    if(model == "MFAC"){
      fixed_nu <- rep(1e-6, g)
      fixed_eta <- rep(1 + 1e-4, g)
      params$nu <- fixed_nu
      params$eta <- fixed_eta
    }
    
    rmfac_check_params(params)
    ll_old <- rmfac_loglikelihood(l_mat, u_mat, X, params)
    if(!is.finite(ll_old)) stop("Initial log-likelihood is not finite.")
    
    converged <- FALSE
    min_delta <- Inf
    
    for(iter in seq_len(max_iter)){
      estep <- rmfac_E_step(l_mat, u_mat, X, params)
      fact <- rmfac_factor_posteriors(params, estep, X)
      estep$fhat <- fact$fhat
      if(!is.null(fact$Fhat)) estep$Fhat <- fact$Fhat
      
      if(model == "MCNFAC"){
        params <- rmfac_CM1_update(params, estep, X,
                                   eta_eps = ETA_LOWER - 1)
        params <- project_mcnfac_parameters(params)
      }else{
        pi_new <- pmax(colMeans(estep$zhat), 1e-12)
        params$pi <- pi_new / sum(pi_new)
        params$nu <- fixed_nu
        params$eta <- fixed_eta
      }
      
      params_old <- params
      params <- rmfac_CM2_update_beta(params, estep, X)
      params <- rmfac_CM2_update_B(params_old, params, estep, X)
      params <- rmfac_CM2_update_D(params_old, params, estep, X,
                                   d_floor = D_FLOOR)
      
      if(model == "MCNFAC"){
        params <- project_mcnfac_parameters(params)
      }else{
        params$nu <- fixed_nu
        params$eta <- fixed_eta
      }
      
      rmfac_check_params(params)
      ll_new <- rmfac_loglikelihood(l_mat, u_mat, X, params)
      if(!is.finite(ll_new)) stop("Non-finite log-likelihood during iteration.")
      
      delta <- ll_new - ll_old
      min_delta <- min(min_delta, delta)
      if(delta < -MONOTONE_TOL){
        stop(sprintf("Log-likelihood decreased by %.6e at iteration %d.",
                     delta, iter))
      }
      
      rel_change <- abs(delta) / (abs(ll_old) + 1e-10)
      ll_old <- ll_new
      if(rel_change < tol){
        converged <- TRUE
        break
      }
    }
    
    elapsed <- proc.time()[["elapsed"]] - start_time
    ll_std <- ll_old
    ll_original <- to_original_loglik(ll_std)
    
    if(!converged){
      record <- make_start_record(
        model, g, q, init, start_type, start_nu, start_eta,
        converged = FALSE, admissible = FALSE, iter = iter,
        ll_std = ll_std, ll_original = ll_original,
        min_delta = min_delta, elapsed = elapsed,
        reason = "maximum iterations reached", error = ""
      )
      return(list(fit = NULL, record = record))
    }
    
    estep_final <- rmfac_E_step(l_mat, u_mat, X, params)
    diagnostics <- assess_solution(params, estep_final, model, q)
    k <- if(model == "MCNFAC")
      count_params_mcnfac(g, ncol(l_mat), ncol(X), q) else
        count_params_mfac(g, ncol(l_mat), ncol(X), q)
    IC <- info_criteria(ll_original, k, nrow(X))
    entropy <- -mean(rowSums(
      estep_final$zhat * log(estep_final$zhat + 1e-12)
    ))
    
    fit <- list(
      model = model, g = g, q = q,
      converged = TRUE,
      admissible = diagnostics$admissible,
      admissibility_reason = diagnostics$reason,
      logLik_standardized = ll_std,
      logLik = ll_original,
      AIC = IC$AIC, BIC = IC$BIC, EDC = IC$EDC,
      k = k, entropy = entropy, iter = iter,
      cpu_time = elapsed, params = params,
      zhat = get_map_class(estep_final), estep = estep_final,
      initialization = init, start_type = start_type,
      start_nu = start_nu, start_eta = start_eta,
      diagnostics = diagnostics,
      min_logLik_increment = min_delta
    )
    
    record <- make_start_record(
      model, g, q, init, start_type, start_nu, start_eta,
      converged = TRUE, admissible = diagnostics$admissible, iter = iter,
      ll_std = ll_std, ll_original = ll_original,
      min_delta = min_delta, elapsed = elapsed,
      reason = diagnostics$reason, error = "", diagnostics = diagnostics
    )
    list(fit = fit, record = record)
  }, error = function(e){
    elapsed <- proc.time()[["elapsed"]] - start_time
    record <- make_start_record(
      model, g, q, init, start_type, start_nu, start_eta,
      elapsed = elapsed, reason = "numerical failure",
      error = conditionMessage(e)
    )
    list(fit = NULL, record = record)
  })
}

fit_candidate_robust <- function(model, l_mat, u_mat, X, g, q,
                                 n_init = N_INIT, mfac_anchor = NULL){
  total_start <- proc.time()[["elapsed"]]
  start_results <- vector("list", n_init)
  
  anchor_grid <- data.frame(
    nu = c(NU_LOWER, 0.01, 0.05, 0.05, 0.10, 0.20),
    eta = c(ETA_LOWER, 2, 2, 5, 5, 10)
  )
  valid_anchor <- model == "MCNFAC" && !is.null(mfac_anchor) &&
    isTRUE(mfac_anchor$converged) && isTRUE(mfac_anchor$admissible) &&
    !is.null(mfac_anchor$params)
  n_anchor <- if(valid_anchor) min(nrow(anchor_grid), n_init) else 0L
  
  for(init in seq_len(n_init)){
    if(valid_anchor && init <= n_anchor){
      start_results[[init]] <- fit_single_start(
        model, l_mat, u_mat, X, g, q, init,
        start_type = "MFAC_anchor",
        start_nu = anchor_grid$nu[init],
        start_eta = anchor_grid$eta[init],
        anchor = mfac_anchor
      )
    }else{
      start_results[[init]] <- fit_single_start(
        model, l_mat, u_mat, X, g, q, init,
        start_type = "standard"
      )
    }
  }
  
  init_table <- do.call(rbind, lapply(start_results, `[[`, "record"))
  converged_fits <- Filter(
    function(x) !is.null(x) && isTRUE(x$converged),
    lapply(start_results, `[[`, "fit")
  )
  admissible_fits <- Filter(
    function(x) isTRUE(x$admissible), converged_fits
  )
  
  total_elapsed <- proc.time()[["elapsed"]] - total_start
  pool <- if(length(admissible_fits)) admissible_fits else converged_fits
  
  if(!length(pool)){
    k <- if(model == "MCNFAC")
      count_params_mcnfac(g, ncol(l_mat), ncol(X), q) else
        count_params_mfac(g, ncol(l_mat), ncol(X), q)
    return(list(
      model = model, g = g, q = q,
      converged = FALSE, admissible = FALSE,
      admissibility_reason = "no converged initialization",
      logLik_standardized = NA_real_, logLik = NA_real_,
      AIC = NA_real_, BIC = NA_real_, EDC = NA_real_, k = k,
      entropy = NA_real_, iter = NA_integer_, cpu_time = NA_real_,
      total_elapsed = total_elapsed, initialization = NA_integer_,
      params = NULL, zhat = NULL, estep = NULL,
      init_table = init_table,
      n_starts = n_init,
      n_converged_starts = sum(init_table$converged),
      n_admissible_starts = sum(init_table$admissible)
    ))
  }
  
  ll_values <- vapply(pool, function(x) x$logLik, numeric(1))
  best_fit <- pool[[which.max(ll_values)]]
  best_fit$total_elapsed <- total_elapsed
  best_fit$init_table <- init_table
  best_fit$n_starts <- n_init
  best_fit$n_converged_starts <- sum(init_table$converged)
  best_fit$n_admissible_starts <- sum(init_table$admissible)
  
  if(!length(admissible_fits)){
    best_fit$admissible <- FALSE
    best_fit$admissibility_reason <- paste0(
      "no admissible initialization; best converged start: ",
      best_fit$admissibility_reason
    )
  }
  best_fit
}

fit_mfac_once <- function(l_mat, u_mat, X, g, q, n_init = N_INIT,
                          max_iter = MAX_ITER, tol = TOL){
  fit_candidate_robust("MFAC", l_mat, u_mat, X, g, q, n_init, NULL)
}

fit_mcnfac_once <- function(l_mat, u_mat, X, g, q, n_init = N_INIT,
                            max_iter = MAX_ITER, tol = TOL,
                            mfac_anchor = NULL){
  fit_candidate_robust("MCNFAC", l_mat, u_mat, X, g, q, n_init, mfac_anchor)
}

###############################################################
# 8. FIT ALL CANDIDATE MODELS
#
# IMPORTANT:
# - MFAC candidates are fitted first.
# - For each MCNFAC (g,q), the corresponding fitted MFAC (g,q)
#   is supplied as an initialization anchor.
# - MCNFAC still also uses ordinary rmfac_initialize() starts.
###############################################################

models <- c("MFAC", "MCNFAC")

checkpoint_file <- file.path(OUT_DIR, "candidate_fit_checkpoint.rds")
fits <- list()
if(RESUME_FROM_CHECKPOINT && file.exists(checkpoint_file)){
  checkpoint_object <- try(readRDS(checkpoint_file), silent = TRUE)
  if(!inherits(checkpoint_object, "try-error") && is.list(checkpoint_object)){
    fits <- checkpoint_object
    cat("Loaded", length(fits), "completed candidate(s) from checkpoint.\n")
  }else{
    warning("The candidate checkpoint could not be read and will be ignored.")
  }
}

summary_rows <- list()
init_rows <- list()

row_counter  <- 1
init_counter <- 1


cat("\n=============================================\n")
cat("FITTING LAKE MICHIGAN 2017 DATA\n")
cat("=============================================\n")

cat(
  "Models:",
  paste(models, collapse = ", "),
  "\n"
)

cat(
  "g grid:",
  paste(G_GRID, collapse = ", "),
  "\n"
)

cat(
  "q grid:",
  paste(Q_GRID, collapse = ", "),
  "\n"
)

cat(
  "Initializations per candidate:",
  N_INIT,
  "\n"
)

cat(
  "MAX_ITER =",
  MAX_ITER,
  "| TOL =",
  TOL,
  "\n\n"
)


###############################################################
# Fit candidates
###############################################################

for(model in models){
  
  for(g in G_GRID){
    
    for(q in Q_GRID){
      
      #########################################################
      # Reproducible candidate-specific seed
      #########################################################
      
      candidate_seed <-
        20260807 +
        1000 * g +
        10 * q
      
      set.seed(candidate_seed)
      
      
      #########################################################
      # Candidate name
      #########################################################
      
      key <- paste0(
        model,
        "_g",
        g,
        "_q",
        q
      )
      
      
      cat("\n---------------------------------------------\n")
      cat("Fitting", key, "\n")
      cat("---------------------------------------------\n")
      
      candidate_already_saved <- key %in% names(fits)
      if(candidate_already_saved){
        fit <- fits[[key]]
        cat("Using completed candidate from checkpoint.\n")
      }
      
      
      #########################################################
      # For MCNFAC:
      # retrieve corresponding MFAC(g,q) fit as anchor
      #########################################################
      
      mfac_anchor <- NULL
      
      if(!candidate_already_saved && model == "MCNFAC"){
        
        mfac_key <- paste0(
          "MFAC_g",
          g,
          "_q",
          q
        )
        
        if(mfac_key %in% names(fits)){
          
          candidate_mfac <- fits[[mfac_key]]
          
          if(
            !inherits(candidate_mfac, "fit_error") &&
            isTRUE(candidate_mfac$converged) &&
            isTRUE(candidate_mfac$admissible) &&
            !is.null(candidate_mfac$params) &&
            is.finite(candidate_mfac$logLik)
          ){
            
            mfac_anchor <- candidate_mfac
            
            cat(
              "Using",
              mfac_key,
              "as MCNFAC initialization anchor.\n"
            )
            
            cat(
              "MFAC anchor logLik =",
              round(candidate_mfac$logLik, 6),
              "\n"
            )
            
          }else{
            
            cat(
              "Corresponding MFAC fit is unavailable, inadmissible, or did not converge.\n",
              "MCNFAC will use standard initializations only.\n"
            )
          }
          
        }else{
          
          cat(
            "Corresponding MFAC candidate not found.\n",
            "MCNFAC will use standard initializations only.\n"
          )
        }
      }
      
      
      #########################################################
      # Fit candidate
      #########################################################
      
      if(!candidate_already_saved){
        fit <- tryCatch(
          
          fit_model_once(
            
            model = model,
            
            l_mat = l_mat,
            
            u_mat = u_mat,
            
            X = X,
            
            g = g,
            
            q = q,
            
            n_init = N_INIT,
            
            mfac_anchor = mfac_anchor
            
          ),
          
          error = function(e){
            
            structure(
              
              list(
                
                model = model,
                
                g = g,
                
                q = q,
                
                converged = FALSE,
                
                message = conditionMessage(e)
                
              ),
              
              class = "fit_error"
              
            )
          }
        )
        
        
        #########################################################
        # Store fit
        #########################################################
        
        fits[[key]] <- fit
        saveRDS(fits, checkpoint_file)
      }
      
      
      #########################################################
      # Fit failed
      #########################################################
      
      if(inherits(fit, "fit_error")){
        
        summary_rows[[row_counter]] <- data.frame(
          
          model = model,
          
          g = g,
          
          q = q,
          
          converged = FALSE,
          
          admissible = FALSE,
          
          admissibility_reason = fit$message,
          
          initialization = NA_integer_,
          
          start_type = NA_character_,
          
          start_nu = NA_real_,
          
          start_eta = NA_real_,
          
          iter = NA_integer_,
          
          logLik = NA_real_,
          
          logLik_standardized = NA_real_,
          
          k =
            if(model == "MCNFAC"){
              
              count_params_mcnfac(
                g,
                ncol(l_mat),
                ncol(X),
                q
              )
              
            }else{
              
              count_params_mfac(
                g,
                ncol(l_mat),
                ncol(X),
                q
              )
            },
          
          AIC = NA_real_,
          
          BIC = NA_real_,
          
          EDC = NA_real_,
          
          entropy = NA_real_,
          
          winning_init_seconds = NA_real_,
          
          total_elapsed_seconds = NA_real_,
          
          n_starts = NA_integer_,
          
          n_converged_starts = NA_integer_,
          
          n_admissible_starts = NA_integer_,
          
          minimum_required_component_size = NA_integer_,
          
          min_mixing_proportion = NA_real_,
          
          min_effective_component_size = NA_real_,
          
          min_MAP_component_size = NA_integer_,
          
          min_uniqueness = NA_real_,
          
          min_covariance_eigenvalue = NA_real_,
          
          max_correlation_condition_number = NA_real_,
          
          min_contamination_proportion = NA_real_,
          
          max_contamination_proportion = NA_real_,
          
          min_inflation_factor = NA_real_,
          
          max_inflation_factor = NA_real_,
          
          error = fit$message,
          
          stringsAsFactors = FALSE
          
        )
        
        
      }else{
        
        #######################################################
        # Successful / returned fit
        #######################################################
        
        summary_rows[[row_counter]] <- data.frame(
          
          model = fit$model,
          
          g = fit$g,
          
          q = fit$q,
          
          converged = fit$converged,
          
          admissible = isTRUE(fit$admissible),
          
          admissibility_reason = if(!is.null(fit$admissibility_reason))
            fit$admissibility_reason else NA_character_,
          
          initialization = fit$initialization,
          
          start_type =
            if(!is.null(fit$start_type))
              fit$start_type
          else
            if(model == "MFAC")
              "standard"
          else
            NA_character_,
          
          start_nu =
            if(!is.null(fit$start_nu))
              fit$start_nu
          else
            NA_real_,
          
          start_eta =
            if(!is.null(fit$start_eta))
              fit$start_eta
          else
            NA_real_,
          
          iter = fit$iter,
          
          logLik = fit$logLik,
          
          logLik_standardized = fit$logLik_standardized,
          
          k = fit$k,
          
          AIC = fit$AIC,
          
          BIC = fit$BIC,
          
          EDC = fit$EDC,
          
          entropy = fit$entropy,
          
          winning_init_seconds = fit$cpu_time,
          
          total_elapsed_seconds = fit$total_elapsed,
          
          n_starts = fit$n_starts,
          
          n_converged_starts = fit$n_converged_starts,
          
          n_admissible_starts = fit$n_admissible_starts,
          
          minimum_required_component_size = if(is.null(fit$diagnostics)) NA_integer_ else
            fit$diagnostics$minimum_required_component_size,
          
          min_mixing_proportion = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_mixing_proportion,
          
          min_effective_component_size = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_effective_component_size,
          
          min_MAP_component_size = if(is.null(fit$diagnostics)) NA_integer_ else
            fit$diagnostics$min_MAP_component_size,
          
          min_uniqueness = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_uniqueness,
          
          min_covariance_eigenvalue = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_covariance_eigenvalue,
          
          max_correlation_condition_number = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$max_correlation_condition_number,
          
          min_contamination_proportion = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_contamination_proportion,
          
          max_contamination_proportion = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$max_contamination_proportion,
          
          min_inflation_factor = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$min_inflation_factor,
          
          max_inflation_factor = if(is.null(fit$diagnostics)) NA_real_ else
            fit$diagnostics$max_inflation_factor,
          
          error = "",
          
          stringsAsFactors = FALSE
          
        )
        
        
        #######################################################
        # Save every initialization
        #######################################################
        
        if(!is.null(fit$init_table)){
          
          init_rows[[init_counter]] <-
            fit$init_table
          
          init_counter <-
            init_counter + 1
        }
      }
      
      
      #########################################################
      # Console summary for candidate
      #########################################################
      
      if(
        !inherits(fit, "fit_error") &&
        isTRUE(fit$converged)
      ){
        
        cat(
          "Finished:",
          key,
          "\n"
        )
        
        cat(
          "Best logLik =",
          round(fit$logLik, 6),
          "\n"
        )
        
        cat(
          "BIC =",
          round(fit$BIC, 6),
          "\n"
        )
        
        cat(
          "Admissible =", isTRUE(fit$admissible),
          "| converged starts =", fit$n_converged_starts, "of", fit$n_starts,
          "| admissible starts =", fit$n_admissible_starts, "\n"
        )
        
        if(!isTRUE(fit$admissible))
          cat("Diagnostic reason:", fit$admissibility_reason, "\n")
        
        cat(
          "Winning initialization =",
          fit$initialization,
          "\n"
        )
        
        if(model == "MCNFAC"){
          
          cat(
            "Winning start type =",
            fit$start_type,
            "\n"
          )
          
          if(
            identical(
              fit$start_type,
              "MFAC_anchor"
            )
          ){
            
            cat(
              "Starting nu =",
              fit$start_nu,
              "| starting eta =",
              fit$start_eta,
              "\n"
            )
          }
        }
        
      }else{
        
        cat(
          "Candidate did not produce a converged fit:",
          key,
          "\n"
        )
      }
      
      
      row_counter <-
        row_counter + 1
    }
  }
}


###############################################################
# Combine candidate summaries
###############################################################

model_table <- do.call(
  rbind,
  summary_rows
)

model_table$key <- paste0(model_table$model, "_g", model_table$g,
                          "_q", model_table$q)

for(criterion in c("AIC", "BIC", "EDC")){
  usable <- model_table$converged & model_table$admissible &
    is.finite(model_table[[criterion]])
  model_table[[paste0("delta_", criterion)]] <- NA_real_
  model_table[[paste0("rank_", criterion)]] <- NA_integer_
  if(any(usable)){
    model_table[[paste0("delta_", criterion)]][usable] <-
      model_table[[criterion]][usable] - min(model_table[[criterion]][usable])
    model_table[[paste0("rank_", criterion)]][usable] <-
      rank(model_table[[criterion]][usable], ties.method = "min")
  }
}

model_table <- model_table[
  order(
    model_table$model,
    model_table$BIC
  ),
  ,
  drop = FALSE
]

row.names(model_table) <- NULL


###############################################################
# Save candidate-model table
###############################################################

write.csv(
  
  model_table,
  
  file.path(
    OUT_DIR,
    "03_all_candidate_models.csv"
  ),
  
  row.names = FALSE
)


###############################################################
# Save all initialization results
###############################################################

if(length(init_rows) > 0){
  
  init_table_all <- do.call(
    rbind,
    init_rows
  )
  
  write.csv(
    
    init_table_all,
    
    file.path(
      OUT_DIR,
      "04_all_initializations.csv"
    ),
    
    row.names = FALSE
  )
}


###############################################################
# Print candidate results
###############################################################

cat("\n=============================================\n")
cat("ALL CANDIDATE MODELS\n")
cat("=============================================\n")

print(model_table)


###############################################################
# Additional paired likelihood diagnostic
###############################################################

cat("\n=============================================\n")
cat("PAIRWISE MCNFAC vs MFAC CHECK\n")
cat("=============================================\n")

for(g in G_GRID){
  
  for(q in Q_GRID){
    
    key_mfac <- paste0(
      "MFAC_g",
      g,
      "_q",
      q
    )
    
    key_mcnfac <- paste0(
      "MCNFAC_g",
      g,
      "_q",
      q
    )
    
    if(
      key_mfac %in% names(fits) &&
      key_mcnfac %in% names(fits)
    ){
      
      f0 <- fits[[key_mfac]]
      f1 <- fits[[key_mcnfac]]
      
      if(
        !inherits(f0, "fit_error") &&
        !inherits(f1, "fit_error") &&
        isTRUE(f0$converged) &&
        isTRUE(f1$converged) &&
        is.finite(f0$logLik) &&
        is.finite(f1$logLik)
      ){
        
        diff_ll <-
          f1$logLik -
          f0$logLik
        
        cat(
          "g =", g,
          "| q =", q,
          "| MFAC =", round(f0$logLik, 6),
          "| MCNFAC =", round(f1$logLik, 6),
          "| difference =", round(diff_ll, 6),
          "\n"
        )
        
        if(diff_ll < -1e-4){
          
          warning(
            sprintf(
              paste0(
                "MCNFAC still below MFAC for g=%d q=%d ",
                "(difference = %.6f)."
              ),
              g,
              q,
              diff_ll
            )
          )
        }
      }
    }
  }
}

###############################################################
# 9. REPORT CRITERION-SPECIFIC WINNERS
###############################################################

select_best_key <- function(model_name, criterion = "BIC"){
  tab <- model_table[
    model_table$model == model_name &
      model_table$converged &
      model_table$admissible &
      is.finite(model_table[[criterion]]),
    , drop = FALSE
  ]
  
  if(nrow(tab) == 0)
    stop("No converged admissible candidate model for ", model_name)
  
  best_row <- tab[which.min(tab[[criterion]]), , drop = FALSE]
  list(
    row = best_row,
    key = paste0(model_name, "_g", best_row$g, "_q", best_row$q)
  )
}

criterion_names <- c("AIC", "BIC", "EDC")

criterion_winners_by_family <- do.call(
  rbind,
  lapply(criterion_names, function(criterion){
    do.call(rbind, lapply(c("MFAC", "MCNFAC"), function(model_name){
      ans <- select_best_key(model_name, criterion)
      data.frame(
        criterion = criterion,
        family = model_name,
        key = ans$key,
        g = ans$row$g,
        q = ans$row$q,
        logLik = ans$row$logLik,
        k = ans$row$k,
        criterion_value = ans$row[[criterion]],
        stringsAsFactors = FALSE
      )
    }))
  })
)

select_overall <- function(criterion){
  tab <- model_table[
    model_table$converged & model_table$admissible &
      is.finite(model_table[[criterion]]), , drop = FALSE
  ]
  if(nrow(tab) == 0)
    stop("No converged admissible candidate for ", criterion, ".")
  ans <- tab[which.min(tab[[criterion]]), , drop = FALSE]
  data.frame(
    criterion = criterion,
    family = ans$model,
    key = ans$key,
    g = ans$g,
    q = ans$q,
    logLik = ans$logLik,
    k = ans$k,
    criterion_value = ans[[criterion]],
    stringsAsFactors = FALSE
  )
}

criterion_winners_overall <- do.call(
  rbind, lapply(criterion_names, select_overall)
)

write.csv(
  criterion_winners_by_family,
  file.path(OUT_DIR, "05_criterion_winners_by_family.csv"),
  row.names = FALSE
)
write.csv(
  criterion_winners_overall,
  file.path(OUT_DIR, "06_criterion_winners_overall.csv"),
  row.names = FALSE
)

criteria_agree_overall <- length(unique(criterion_winners_overall$key)) == 1
criteria_agree_by_family <- tapply(
  criterion_winners_by_family$key,
  criterion_winners_by_family$family,
  function(x) length(unique(x)) == 1
)

cat("\n===== CRITERION-SPECIFIC WINNERS WITHIN EACH FAMILY =====\n")
print(criterion_winners_by_family)
cat("\n===== CRITERION-SPECIFIC WINNERS OVERALL =====\n")
print(criterion_winners_overall)
if(!criteria_agree_overall){
  cat("\nAIC, BIC, and EDC do not select one common overall model.\n")
  cat("Accordingly, this script does not declare a unique final model.\n")
}

# BIC representatives are retained only to produce comparable downstream
# component and contamination summaries. They are not automatically declared
# to be the final scientific models.
best_mfac_info   <- select_best_key("MFAC", PROVISIONAL_CRITERION)
best_mcnfac_info <- select_best_key("MCNFAC", PROVISIONAL_CRITERION)

best_mfac   <- fits[[best_mfac_info$key]]
best_mcnfac <- fits[[best_mcnfac_info$key]]

best_table <- rbind(best_mfac_info$row, best_mcnfac_info$row)
row.names(best_table) <- NULL

write.csv(
  best_table,
  file.path(OUT_DIR, "07_provisional_BIC_family_representatives.csv"),
  row.names = FALSE
)

cat("\n===== PROVISIONAL FAMILY REPRESENTATIVES BY",
    PROVISIONAL_CRITERION, "=====\n")
print(best_table)

###############################################################
# 10. BASIC NESTED-MODEL SANITY CHECK
###############################################################

# For the same g and q, MCNFAC contains MFAC as a limiting/nested
# Gaussian case. A materially smaller optimized MCNFAC likelihood
# is therefore a warning sign about local maxima or convergence.
paired_ll <- merge(
  model_table[, c("model", "g", "q", "converged", "logLik")],
  model_table[, c("model", "g", "q", "converged", "logLik")],
  by = c("g", "q"), suffixes = c("_1", "_2")
)
paired_ll <- paired_ll[
  paired_ll$model_1 == "MCNFAC" & paired_ll$model_2 == "MFAC",
  , drop = FALSE
]
paired_ll$MCNFAC_minus_MFAC_logLik <- paired_ll$logLik_1 - paired_ll$logLik_2

write.csv(
  paired_ll,
  file.path(OUT_DIR, "08_paired_loglikelihood_check.csv"),
  row.names = FALSE
)

if(any(
  paired_ll$converged_1 & paired_ll$converged_2 &
  paired_ll$MCNFAC_minus_MFAC_logLik < -1e-4,
  na.rm = TRUE
)){
  warning(
    "At least one converged MCNFAC candidate has a lower log-likelihood than the corresponding MFAC candidate. ",
    "Inspect local maxima/initializations before publication."
  )
}

###############################################################
# 11. EXTRACT POSTERIORS AND CONTAMINATION INFORMATION
###############################################################

extract_subject_results <- function(fit, model_name){
  if(is.null(fit$estep))
    fit$estep <- rmfac_E_step(l_mat, u_mat, X, fit$params)
  
  estep <- fit$estep
  zmat <- estep$zhat
  map_cluster <- apply(zmat, 1, which.max)
  max_cluster_prob <- apply(zmat, 1, max)
  
  out <- data.frame(
    RowID = lake$RowID,
    Date_UTC = lake$Date_UTC,
    Time_UTC = lake$Time_UTC,
    Site = lake$Site,
    Depth_Site = lake$Depth_Site,
    Depth_Smp = lake$Depth_Smp,
    Lat = lake$Lat,
    Long = lake$Long,
    model = model_name,
    cluster = map_cluster,
    max_cluster_probability = max_cluster_prob,
    stringsAsFactors = FALSE
  )
  
  for(i in seq_len(ncol(zmat))){
    out[[paste0("posterior_cluster_", i)]] <- zmat[, i]
  }
  
  if(model_name == "MCNFAC"){
    # vhat_ij = posterior probability of the typical state within
    # component i. Hence 1-vhat_ij is contamination probability.
    contamination_by_component <- 1 - estep$vhat
    
    # Marginal contamination probability, averaging component-specific
    # contamination responsibilities with zhat_ij.
    marginal_contam_prob <- rowSums(zmat * contamination_by_component)
    
    # Contamination probability within the MAP-assigned component.
    map_contam_prob <- vapply(
      seq_len(nrow(zmat)),
      function(j) contamination_by_component[j, map_cluster[j]],
      numeric(1)
    )
    
    out$contamination_probability <- marginal_contam_prob
    out$map_component_contamination_probability <- map_contam_prob
    out$flag_contaminated_0.5 <- marginal_contam_prob > 0.5
    
    for(i in seq_len(ncol(zmat))){
      out[[paste0("contam_prob_component_", i)]] <- contamination_by_component[, i]
    }
  }
  
  out
}

subject_mfac <- extract_subject_results(best_mfac, "MFAC")
subject_mcnfac <- extract_subject_results(best_mcnfac, "MCNFAC")

write.csv(
  subject_mfac,
  file.path(OUT_DIR, "09_provisional_BIC_MFAC_observation_results.csv"),
  row.names = FALSE
)

write.csv(
  subject_mcnfac,
  file.path(OUT_DIR, "10_provisional_BIC_MCNFAC_observation_results.csv"),
  row.names = FALSE
)

###############################################################
# 12. PROVISIONAL BIC-MCNFAC PARAMETER SUMMARIES
###############################################################

mcnfac_component_summary <- data.frame(
  component = seq_len(best_mcnfac$g),
  pi = best_mcnfac$params$pi,
  nu = best_mcnfac$params$nu,
  eta = best_mcnfac$params$eta,
  MAP_cluster_size = as.integer(
    table(factor(subject_mcnfac$cluster, levels = seq_len(best_mcnfac$g)))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mcnfac_component_summary,
  file.path(OUT_DIR, "11_provisional_BIC_MCNFAC_component_summary.csv"),
  row.names = FALSE
)

# Regression coefficients on the original latitude/longitude scale, one row
# per response x covariate x component. The fitted likelihood uses centered
# and scaled spatial covariates, so the transformation below is algebraic and
# does not refit the model.
beta_rows <- list()
bc <- 1
for(i in seq_len(best_mcnfac$g)){
  beta_standardized_X <- matrix(
    best_mcnfac$params$beta[[i]],
    nrow = length(response_names),
    ncol = ncol(X),
    byrow = TRUE
  )
  
  beta_original_X <- beta_standardized_X
  beta_original_X[, 2] <-
    y_scale * beta_standardized_X[, 2] / x_scale[1]
  beta_original_X[, 3] <-
    y_scale * beta_standardized_X[, 3] / x_scale[2]
  beta_original_X[, 1] <- y_center +
    y_scale * beta_standardized_X[, 1] -
    beta_original_X[, 2] * x_center[1] -
    beta_original_X[, 3] * x_center[2]
  colnames(beta_original_X) <- c("Intercept", "Lat", "Long")
  
  for(r in seq_along(response_names)){
    for(cc in seq_len(ncol(X))){
      beta_rows[[bc]] <- data.frame(
        component = i,
        response = response_names[r],
        covariate = colnames(beta_original_X)[cc],
        beta = beta_original_X[r, cc],
        stringsAsFactors = FALSE
      )
      bc <- bc + 1
    }
  }
}
beta_table <- do.call(rbind, beta_rows)
write.csv(
  beta_table,
  file.path(OUT_DIR, "12_provisional_BIC_MCNFAC_beta_original_coordinates.csv"),
  row.names = FALSE
)

# Save B, D and Sigma for each component.
for(i in seq_len(best_mcnfac$g)){
  response_scale_matrix <- diag(y_scale)
  B_i <- response_scale_matrix %*% best_mcnfac$params$B[[i]]
  D_i <- response_scale_matrix %*% best_mcnfac$params$D[[i]] %*%
    response_scale_matrix
  Sigma_i <- rmfac_Sigma_i(B_i, D_i)
  
  write.csv(B_i,
            file.path(OUT_DIR, sprintf("13_B_component_%d.csv", i)),
            row.names = FALSE)
  write.csv(D_i,
            file.path(OUT_DIR, sprintf("14_D_component_%d.csv", i)),
            row.names = FALSE)
  write.csv(Sigma_i,
            file.path(OUT_DIR, sprintf("15_Sigma_component_%d.csv", i)),
            row.names = FALSE)
}

###############################################################
# 13. CONTAMINATED-OBSERVATION TABLE
###############################################################

contam_table <- subject_mcnfac[
  order(subject_mcnfac$contamination_probability, decreasing = TRUE),
  , drop = FALSE
]

# Attach the six original chemistry measurements and censoring
# indicators using RowID, because contam_table has been reordered.
contam_table <- cbind(
  contam_table,
  lake[contam_table$RowID, response_names, drop = FALSE]
)

for(v in response_names){
  contam_table[[paste0(v, "_censored")]] <-
    cens_raw[contam_table$RowID, v]
}

write.csv(
  contam_table,
  file.path(OUT_DIR, "16_provisional_BIC_MCNFAC_contamination_ranking.csv"),
  row.names = FALSE
)

###############################################################
# 14. SIMPLE PUBLICATION-READY FIGURES
###############################################################

if(HAS_GGPLOT2){
  library(ggplot2)
  
  # Spatial clustering plot
  p_cluster <- ggplot(
    subject_mcnfac,
    aes(x = Long, y = Lat)
  )  +
    geom_point(aes(shape = factor(cluster)), size = 3) +
    labs(
      x = "Longitude",
      y = "Latitude",
      shape = "MCNFAC cluster",
      title = "Lake Michigan 2017: MCNFAC clustering"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(OUT_DIR, "Figure1_MCNFAC_spatial_clusters.png"),
    p_cluster, width = 7, height = 5, dpi = 400
  )
  
  # Contamination-probability spatial plot
  p_contam <- ggplot(
    subject_mcnfac,
    aes(x = Long, y = Lat)
  ) +
    geom_point(
      aes(size = contamination_probability,
          shape = flag_contaminated_0.5),
      alpha = 0.85
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      size = "Posterior contamination\nprobability",
      shape = "Posterior > 0.5",
      title = "Lake Michigan 2017: MCNFAC contamination probabilities"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(OUT_DIR, "Figure2_MCNFAC_spatial_contamination.png"),
    p_contam, width = 7, height = 5, dpi = 400
  )
  
  # Information-criterion profiles by g and q. Showing all three criteria is
  # important here because agreement must be assessed rather than assumed.
  criterion_plot_data <- do.call(
    rbind,
    lapply(c("AIC", "BIC", "EDC"), function(criterion){
      tmp <- model_table[
        model_table$converged,
        c("model", "g", "q", "admissible", criterion)
      ]
      names(tmp)[5] <- "value"
      tmp$criterion <- criterion
      tmp
    })
  )
  
  p_criteria <- ggplot(
    criterion_plot_data,
    aes(x = g, y = value, linetype = factor(q),
        group = interaction(model, q))
  ) +
    geom_line(
      data = subset(criterion_plot_data, admissible),
      na.rm = TRUE
    ) +
    geom_point(aes(shape = admissible), size = 2.5, na.rm = TRUE) +
    facet_grid(criterion ~ model, scales = "free_y") +
    scale_x_continuous(breaks = G_GRID) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
    labs(
      x = "Number of mixture components (g)",
      y = "Criterion value",
      linetype = "q",
      shape = "Admissible",
      title = "Model-selection profiles for Lake Michigan 2017",
      subtitle = "All converged candidates are shown; lines join admissible fits"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(OUT_DIR, "Figure3_information_criterion_profiles.png"),
    p_criteria, width = 9, height = 9, dpi = 400
  )
}

###############################################################
# 15. SAVE COMPLETE R OBJECTS
###############################################################

saveRDS(
  list(
    lake = lake,
    response_names = response_names,
    detection_limits = DL,
    censor_matrix = cens_raw,
    l_mat = l_mat,
    u_mat = u_mat,
    X = X,
    response_center = y_center,
    response_scale = y_scale,
    observed_loglikelihood_jacobian = loglik_jacobian,
    X_center = x_center,
    X_scale = x_scale,
    numerical_safeguards = list(
      D_FLOOR = D_FLOOR,
      NU_LOWER = NU_LOWER,
      NU_UPPER = NU_UPPER,
      ETA_LOWER = ETA_LOWER,
      ETA_UPPER = ETA_UPPER,
      MONOTONE_TOL = MONOTONE_TOL,
      MAX_CORR_CONDITION = MAX_CORR_CONDITION
    ),
    model_table = model_table,
    criterion_winners_by_family = criterion_winners_by_family,
    criterion_winners_overall = criterion_winners_overall,
    criteria_agree_overall = criteria_agree_overall,
    criteria_agree_by_family = criteria_agree_by_family,
    best_table = best_table,
    fits = fits,
    best_mfac = best_mfac,
    best_mcnfac = best_mcnfac,
    subject_mfac = subject_mfac,
    subject_mcnfac = subject_mcnfac
  ),
  file.path(OUT_DIR, "LakeMichigan_2017_complete_analysis.rds")
)

###############################################################
# 16. FINAL CONSOLE REPORT
###############################################################

cat("\n\n=====================================================\n")
cat("LAKE MICHIGAN 2017 ANALYSIS COMPLETE\n")
cat("=====================================================\n")
cat("Sample size:", nrow(lake), "\n")
cat("Responses:", paste(response_names, collapse = ", "), "\n")
cat("Likelihood response scale: standardized for computation; logLik and information criteria transformed to the original response scale\n")
cat("Covariates: intercept plus standardized Lat and Long\n")
cat("\nCensoring counts:\n")
print(censor_summary)

cat("\nProvisional MFAC representative by", PROVISIONAL_CRITERION, ": g =",
    best_mfac$g, ", q =", best_mfac$q,
    ", logLik =", round(best_mfac$logLik, 4),
    ", BIC =", round(best_mfac$BIC, 4), "\n")

cat("Provisional MCNFAC representative by", PROVISIONAL_CRITERION, ": g =",
    best_mcnfac$g, ", q =", best_mcnfac$q,
    ", logLik =", round(best_mcnfac$logLik, 4),
    ", BIC =", round(best_mcnfac$BIC, 4), "\n")

if(criteria_agree_overall){
  cat("AIC, BIC, and EDC agree on the overall candidate:",
      criterion_winners_overall$key[1], "\n")
}else{
  cat("AIC, BIC, and EDC disagree; no unique final model is declared.\n")
}

cat("\nProvisional BIC-MCNFAC contamination parameters:\n")
print(mcnfac_component_summary)

cat("\nNumber of observations with marginal posterior contamination probability > 0.5:",
    sum(subject_mcnfac$flag_contaminated_0.5), "\n")

cat("\nTotal elapsed fitting time for provisional MFAC representative:",
    round(best_mfac$total_elapsed, 2), "seconds\n")
cat("Total elapsed fitting time for provisional MCNFAC representative:",
    round(best_mcnfac$total_elapsed, 2), "seconds\n")

cat("\nAll output saved in:", OUT_DIR, "\n")
cat("=====================================================\n")
