###############################################################
# MCNFAC / MFAC REAL-DATA ANALYSIS
# Lake Michigan water chemistry, Year 2017
#
# Responses (Y): SRP, TDP, PP, Chl, PC, PN
# Covariates (X): standardized Latitude and Longitude
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
# - Response standardization is based on the observed/censored
#   representation (censored values replaced by their detection
#   limit), not on the unobserved below-limit values.
###############################################################

rm(list = ls())

###############################################################
# 0. USER PATHS
###############################################################

FUNCTION_FILE <- "MCNFAC_Functions.R"
DATA_FILE     <- "737176_v3_lake_michigan_chemistry.csv"
OUT_DIR       <- "MCNFAC_LakeMichigan_2017_results"

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

# Main candidate grid. With n=61 and p=6, keeping q modest is
# deliberate. q=1,2 is adequate for the primary analysis.
G_GRID <- 1:4
Q_GRID <- 1:2

# Main model-selection criterion
SELECTION_CRITERION <- "BIC"

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

# Standardize responses for numerical stability and to prevent
# high-scale variables (especially PC) from dominating k-means
# initialization. Crucially, centers/scales are computed from
# Yc_raw, so actual below-DL values are not used as if observed.
y_center <- colMeans(Yc_raw)
y_scale  <- apply(Yc_raw, 2, sd)

if(any(!is.finite(y_scale)) || any(y_scale <= 0))
  stop("Invalid response scaling constants.")

Y_scaled <- sweep(Yc_raw, 2, y_center, FUN = "-")
Y_scaled <- sweep(Y_scaled, 2, y_scale, FUN = "/")

DL_scaled <- (DL - y_center) / y_scale

# Exact observations: lower = upper = observed standardized value.
# Left-censored observations: (-Inf, standardized detection limit].
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

###############################################################
# 5. COVARIATE MATRIX
###############################################################

# Lin's Lake Michigan suggestion specifies latitude and longitude
# as X. We standardize them for numerical stability. No artificial
# covariates are added.
X_raw <- as.matrix(lake[, covariate_names, drop = FALSE])
storage.mode(X_raw) <- "double"

x_center <- colMeans(X_raw)
x_scale  <- apply(X_raw, 2, sd)

if(any(!is.finite(x_scale)) || any(x_scale <= 0))
  stop("Invalid covariate scaling constants.")

X <- scale(X_raw, center = x_center, scale = x_scale)
X <- as.matrix(X)
colnames(X) <- c("Lat_std", "Long_std")

# Save preprocessing constants for complete reproducibility.
preprocess_constants <- data.frame(
  variable = c(response_names, covariate_names),
  center = c(y_center, x_center),
  scale = c(y_scale, x_scale),
  detection_limit = c(as.numeric(DL), NA_real_, NA_real_),
  stringsAsFactors = FALSE
)
write.csv(
  preprocess_constants,
  file.path(OUT_DIR, "02_preprocessing_constants.csv"),
  row.names = FALSE
)

###############################################################
# 6. FIT MCNFAC
# Improved real-data initialization:
# MFAC-anchored contamination starts + ordinary starts
###############################################################

fit_mcnfac_once <- function(
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
      
      params <- rmfac_initialize(
        l_mat,
        u_mat,
        X,
        g = g,
        q = q
      )
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

fit_mfac_once <- function(
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
    params <- rmfac_initialize(l_mat, u_mat, X, g = g, q = q)
    
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
# 8. FIT ALL CANDIDATE MODELS
#
# IMPORTANT:
# - MFAC candidates are fitted first.
# - For each MCNFAC (g,q), the corresponding fitted MFAC (g,q)
#   is supplied as an initialization anchor.
# - MCNFAC still also uses ordinary rmfac_initialize() starts.
###############################################################

models <- c("MFAC", "MCNFAC")

fits <- list()

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
      
      
      #########################################################
      # For MCNFAC:
      # retrieve corresponding MFAC(g,q) fit as anchor
      #########################################################
      
      mfac_anchor <- NULL
      
      if(model == "MCNFAC"){
        
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
              "Corresponding MFAC fit is unavailable or did not converge.\n",
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
      
      
      #########################################################
      # Fit failed
      #########################################################
      
      if(inherits(fit, "fit_error")){
        
        summary_rows[[row_counter]] <- data.frame(
          
          model = model,
          
          g = g,
          
          q = q,
          
          converged = FALSE,
          
          initialization = NA_integer_,
          
          start_type = NA_character_,
          
          start_nu = NA_real_,
          
          start_eta = NA_real_,
          
          iter = NA_integer_,
          
          logLik = NA_real_,
          
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
          
          k = fit$k,
          
          AIC = fit$AIC,
          
          BIC = fit$BIC,
          
          EDC = fit$EDC,
          
          entropy = fit$entropy,
          
          winning_init_seconds = fit$cpu_time,
          
          total_elapsed_seconds = fit$total_elapsed,
          
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
# 9. SELECT BEST MFAC AND MCNFAC MODELS
###############################################################

select_best_key <- function(model_name, criterion = "BIC"){
  tab <- model_table[
    model_table$model == model_name &
      model_table$converged &
      is.finite(model_table[[criterion]]),
    , drop = FALSE
  ]
  
  if(nrow(tab) == 0)
    stop("No converged candidate model for ", model_name)
  
  best_row <- tab[which.min(tab[[criterion]]), , drop = FALSE]
  list(
    row = best_row,
    key = paste0(model_name, "_g", best_row$g, "_q", best_row$q)
  )
}

best_mfac_info   <- select_best_key("MFAC", SELECTION_CRITERION)
best_mcnfac_info <- select_best_key("MCNFAC", SELECTION_CRITERION)

best_mfac   <- fits[[best_mfac_info$key]]
best_mcnfac <- fits[[best_mcnfac_info$key]]

best_table <- rbind(best_mfac_info$row, best_mcnfac_info$row)
row.names(best_table) <- NULL

write.csv(
  best_table,
  file.path(OUT_DIR, "05_best_models_by_BIC.csv"),
  row.names = FALSE
)

cat("\n===== BEST MODELS BY", SELECTION_CRITERION, "=====\n")
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
  file.path(OUT_DIR, "06_paired_loglikelihood_check.csv"),
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
  file.path(OUT_DIR, "07_best_MFAC_observation_results.csv"),
  row.names = FALSE
)

write.csv(
  subject_mcnfac,
  file.path(OUT_DIR, "08_best_MCNFAC_observation_results.csv"),
  row.names = FALSE
)

###############################################################
# 12. BEST-MCNFAC PARAMETER SUMMARIES
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
  file.path(OUT_DIR, "09_best_MCNFAC_component_summary.csv"),
  row.names = FALSE
)

# Regression coefficients, one row per response x covariate x component.
beta_rows <- list()
bc <- 1
for(i in seq_len(best_mcnfac$g)){
  beta_mat <- matrix(
    best_mcnfac$params$beta[[i]],
    nrow = length(response_names),
    ncol = ncol(X),
    byrow = TRUE
  )
  
  for(r in seq_along(response_names)){
    for(cc in seq_len(ncol(X))){
      beta_rows[[bc]] <- data.frame(
        component = i,
        response = response_names[r],
        covariate = colnames(X)[cc],
        beta = beta_mat[r, cc],
        stringsAsFactors = FALSE
      )
      bc <- bc + 1
    }
  }
}
beta_table <- do.call(rbind, beta_rows)
write.csv(
  beta_table,
  file.path(OUT_DIR, "10_best_MCNFAC_beta.csv"),
  row.names = FALSE
)

# Save B, D and Sigma for each component.
for(i in seq_len(best_mcnfac$g)){
  B_i <- best_mcnfac$params$B[[i]]
  D_i <- best_mcnfac$params$D[[i]]
  Sigma_i <- rmfac_Sigma_i(B_i, D_i)
  
  write.csv(B_i,
            file.path(OUT_DIR, sprintf("11_B_component_%d.csv", i)),
            row.names = FALSE)
  write.csv(D_i,
            file.path(OUT_DIR, sprintf("12_D_component_%d.csv", i)),
            row.names = FALSE)
  write.csv(Sigma_i,
            file.path(OUT_DIR, sprintf("13_Sigma_component_%d.csv", i)),
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
  file.path(OUT_DIR, "14_MCNFAC_contamination_ranking.csv"),
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
  
  # BIC profile by g and q
  p_bic <- ggplot(
    model_table[model_table$converged, ],
    aes(x = g, y = BIC, linetype = factor(q), group = interaction(model, q))
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~model, scales = "free_y") +
    scale_x_continuous(breaks = G_GRID) +
    labs(
      x = "Number of mixture components (g)",
      y = "BIC",
      linetype = "q",
      title = "Model-selection profiles for Lake Michigan 2017"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    file.path(OUT_DIR, "Figure3_BIC_profiles.png"),
    p_bic, width = 8, height = 5, dpi = 400
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
    X_center = x_center,
    X_scale = x_scale,
    model_table = model_table,
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
cat("Covariates: standardized Lat, Long\n")
cat("\nCensoring counts:\n")
print(censor_summary)

cat("\nBest MFAC by", SELECTION_CRITERION, ": g =",
    best_mfac$g, ", q =", best_mfac$q,
    ", logLik =", round(best_mfac$logLik, 4),
    ", BIC =", round(best_mfac$BIC, 4), "\n")

cat("Best MCNFAC by", SELECTION_CRITERION, ": g =",
    best_mcnfac$g, ", q =", best_mcnfac$q,
    ", logLik =", round(best_mcnfac$logLik, 4),
    ", BIC =", round(best_mcnfac$BIC, 4), "\n")

cat("\nBest MCNFAC contamination parameters:\n")
print(mcnfac_component_summary)

cat("\nNumber of observations with marginal posterior contamination probability > 0.5:",
    sum(subject_mcnfac$flag_contaminated_0.5), "\n")

cat("\nTotal elapsed fitting time for selected MFAC candidate:",
    round(best_mfac$total_elapsed, 2), "seconds\n")
cat("Total elapsed fitting time for selected MCNFAC candidate:",
    round(best_mcnfac$total_elapsed, 2), "seconds\n")

cat("\nAll output saved in:", OUT_DIR, "\n")
cat("=====================================================\n")
