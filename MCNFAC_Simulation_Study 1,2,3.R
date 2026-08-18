
###############################################################
# FINAL Simulation Study for MCNFAC following Prof. Lin's advice
#
# Purpose:
#   Compare MFAC vs MCNFAC for:
#     1) parameter recovery under censoring rates 0%,10%,20%,30%
#     2) recovery of true number of mixture components g
#     3) clustering performance
# Required local file:
#source("MCNFAC_Functions.R")
###############################################################

rm(list = ls())

###############################################################
# 0. Load functions and packages
###############################################################
if(!file.exists("MCNFAC_Functions.R"))
  stop("MCNFAC_Functions.R not found.")

source("MCNFAC_Functions.R")

required_pkgs <- c(
  "mvtnorm",
  "mclust",
  "gtools"
)
for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(mvtnorm)
  library(mclust)
  library(gtools)
  library(parallel)
})

set.seed(123)

###############################################################
# 1. USER SETTINGS
###############################################################

R_sim <- 100

studies_to_run <- c(
  "Study1_recovery",
  "Study2_component_selection",
  "Study3_clustering"
)

MAX_ITER <- 600
TOL <- 1e-6

# Parallel settings. Works on Windows using a PSOCK cluster.
# Leave one core free so the computer remains usable.
USE_PARALLEL <- TRUE
n_detect <- parallel::detectCores()

if(is.na(n_detect))
  n_detect <- 2

N_CORES <- max(1, n_detect - 1)

TASK_CHUNK_SIZE <- max(1, 2 * N_CORES)
cat("Parallel setting: USE_PARALLEL =", USE_PARALLEL, "| N_CORES =", N_CORES, "\n")

out_dir <- "MCNFAC_simulation2"
if(!dir.exists(out_dir)) dir.create(out_dir,
                                    recursive = TRUE,
                                    showWarnings = FALSE)

###############################################################
# 2. Parameter counts and information criteria
###############################################################

count_params_mcnfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p + 2)
}

count_params_mfac <- function(g, p, d, q){
  (g - 1) + g * (p*d + p*q - q*(q - 1)/2 + p)
}

info_criteria <- function(logLik, k, n){
  
  if(!is.finite(logLik))
    stop("logLik must be finite.")
  
  if(k < 0)
    stop("k must be non-negative.")
  
  if(n <= 0)
    stop("n must be positive.")
  
  data.frame(
    
    AIC = -2 * logLik + 2 * k,
    
    BIC = -2 * logLik + k * log(n),
    
    EDC = -2 * logLik + k * 0.2 * sqrt(n)
    
  )
}

###############################################################
# 3. Utility functions
###############################################################

make_X <- function(n){
  x4 <- rnorm(n)
  
  if(n > 1)  x4 <- as.numeric(scale(x4))
  
  X <- cbind(
    rbinom(n, 1, 0.5),
    rbinom(n, 1, 0.5),
    rbinom(n, 1, 0.6),
    x4
  )
  colnames(X) <- paste0("X", 1:4)
  X
}

sigma_from_params <- function(params, i) rmfac_Sigma_i(params$B[[i]],params$D[[i]])

get_map_class <- function(estep){
  if(!is.null(estep$zhat)) return(apply(estep$zhat, 1, which.max))
  if(!is.null(estep$z_ij)) return(apply(estep$z_ij, 1, which.max))
  NULL
}

misclass_rate <- function(true_z, pred_z){
  true_z <- as.integer(as.factor(true_z))
  pred_z <- as.integer(as.factor(pred_z))
  tab <- table(true_z, pred_z)
  K <- max(nrow(tab), ncol(tab))
  if(K>8) stop("Too many classes for exhaustive permutation.")
  tab_pad <- matrix(0, K, K)
  tab_pad[1:nrow(tab), 1:ncol(tab)] <- tab
  perms <- gtools::permutations(K, K)
  best_correct <- 0
  for(r in seq_len(nrow(perms))){
    pp <- perms[r, ]
    correct <- sum(tab_pad[cbind(1:K, pp)])
    if(correct > best_correct) best_correct <- correct
  }
  1 - best_correct / length(true_z)
}

bind_rows_fill <- function(x){
  x <- x[!sapply(x, is.null)]
  if(length(x) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(x, names)))
  x2 <- lapply(x, function(df){
    miss <- setdiff(all_names, names(df))
    for(m in miss) df[[m]] <- NA
    df[, all_names, drop = FALSE]
  })
  do.call(rbind, x2)
}

safe_mean <- function(x){
  
  x <- x[is.finite(x)]
  
  if(length(x)==0)
    return(NA_real_)
  
  mean(x)
  
}
###############################################################
# 3b. Parallel helper
###############################################################

# Runs a list of tasks either serially or in parallel, and saves progress
# after each chunk. This keeps the scientific design unchanged, but speeds
# up Monte Carlo replications and model fits across CPU cores.
process_tasks <- function(tasks, task_fun, progress_file){
  if(length(tasks) == 0) return(data.frame())
  
  old_results <- data.frame()
  
  if(file.exists(progress_file)){
    cat("\nExisting progress file found:", progress_file, "\n")
    old_results <- tryCatch(
      read.csv(progress_file),
      error=function(e){
        cat("Progress file corrupted. Starting from scratch.\n")
        data.frame()
      }
    )
  }
  
  if(nrow(old_results) > 0 &&
     all(c("rep", "n", "censor_target") %in% names(old_results))){
    
    done_keys <- unique(paste(old_results$rep,
                              old_results$n,
                              old_results$censor_target,
                              sep = "_"))
    
    task_keys <- sapply(tasks, function(tt){
      paste(tt$rep, tt$n, tt$censor, sep = "_")
    })
    
    tasks <- tasks[!(task_keys %in% done_keys)]
    
    cat("Completed tasks already found:", length(done_keys), "\n")
    cat("Remaining tasks to run:", length(tasks), "\n")
  }
  
  if(length(tasks) == 0){
    cat("All tasks already completed for this study.\n")
    return(old_results)
  }
  
  results <- list()
  if(nrow(old_results) > 0){
    results[[1]] <- old_results
    counter <- 2
  } else {
    counter <- 1
  }
  
  chunks <- split(tasks, ceiling(seq_along(tasks) / TASK_CHUNK_SIZE))
  
  for(cc in seq_along(chunks)){
    cat("\nProcessing chunk", cc, "of", length(chunks),
        "| tasks in chunk =", length(chunks[[cc]]), "\n")
    
    if(isTRUE(USE_PARALLEL) && exists("CL", envir = .GlobalEnv, inherits = FALSE) && !is.null(CL) && N_CORES > 1){
      chunk_res <- parallel::parLapplyLB(
        CL,
        chunks[[cc]],
        function(x){
          
          tryCatch(
            task_fun(x),
            error=function(e){
              
              data.frame(
                rep = x$rep,
                n = x$n,
                censor_target = x$censor,
                error = conditionMessage(e)
              )
              
            }
          )
          
        }
      )
    } else {
      chunk_res <- lapply(
        chunks[[cc]],
        function(task){
          
          tryCatch(
            task_fun(task),
            error=function(e){
              
              data.frame(
                rep=task$rep,
                n=task$n,
                censor_target=task$censor,
                error=conditionMessage(e)
              )
              
            }
          )
          
        }
      )
    }
    
    for(z in chunk_res){
      results[[counter]] <- z
      counter <- counter + 1
    }
    
    write.csv(bind_rows_fill(results), progress_file, row.names=FALSE)
  }
  
  bind_rows_fill(results)
}

###############################################################
# 4. Scenario definitions
###############################################################
###############################################################
# 4. Scenario definitions
###############################################################

make_scenario <- function(name){
  
  #############################################################
  # Scenario 1
  # Moderate contamination
  # (Baseline scenario used in Smoke Test G)
  #############################################################
  
  if(name == "moderate_g2_q1"){
    
    p <- 4
    q <- 1
    g <- 2
    d <- 4
    
    return(list(
      
      name = name,
      
      p = p,
      q = q,
      g = g,
      d = d,
      
      pi = c(0.55,0.45),
      
      beta = list(
        
        c(
          0.0,-0.5,-0.6,-0.4,
          -0.5,-1.0,-0.6,-0.7,
          -1.0,-0.8,-0.7,-0.9,
          -0.4,-0.6,-0.8,-0.5
        ),
        
        c(
          2.0,2.5,2.4,2.3,
          2.5,2.0,2.4,2.3,
          2.2,2.4,2.5,2.1,
          2.3,2.5,2.1,2.4
        )
        
      ),
      
      B = list(
        
        matrix(c(
          1.5,
          1.2,
          0.9,
          0.6
        ), ncol = 1),
        
        matrix(c(
          -0.2,
          0.2,
          0.8,
          1.4
        ), ncol = 1)
        
      ),
      
      D = list(
        
        diag(rep(0.05,4)),
        
        diag(rep(0.05,4))
        
      ),
      
      nu = c(0.10,0.20),
      
      eta = c(5,8)
      
    ))
  }
  
  #############################################################
  # Scenario 2
  # Heavy contamination only
  #############################################################
  
  if(name == "heavy_g2_q1"){
    
    scen <- make_scenario("moderate_g2_q1")
    
    scen$name <- name
    
    scen$nu  <- c(0.25,0.30)
    
    scen$eta <- c(15,20)
    
    return(scen)
    
  }
  
  #############################################################
  # Scenario 3
  # Complex latent structure
  #############################################################
  
  if(name == "complex_g3_q2"){
    
    p <- 6
    q <- 2
    g <- 3
    d <- 4
    
    return(list(
      
      name = name,
      
      p = p,
      q = q,
      g = g,
      d = d,
      
      pi = c(0.40,0.35,0.25),
      
      beta = list(
        
        ####################################################
        # Component 1
        ####################################################
        c(
          0.0,-0.5,-0.6,-0.4,
          -0.5,-1.0,-0.6,-0.7,
          -1.0,-0.8,-0.7,-0.9,
          -0.4,-0.6,-0.8,-0.5,
          -0.8,-0.9,-1.1,-0.7,
          -0.6,-0.5,-0.9,-1.0
        ),
        
        ####################################################
        # Component 2
        ####################################################
        c(
          2.0,2.5,2.4,2.3,
          2.5,2.0,2.4,2.3,
          2.2,2.4,2.5,2.1,
          2.3,2.5,2.1,2.4,
          2.6,2.3,2.4,2.5,
          2.4,2.6,2.3,2.2
        ),
        
        ####################################################
        # Component 3
        ####################################################
        c(
          4.2,4.5,4.4,4.3,
          4.6,4.2,4.5,4.4,
          4.3,4.6,4.2,4.5,
          4.5,4.3,4.6,4.2,
          4.4,4.7,4.5,4.3,
          4.6,4.4,4.3,4.7
        )
        
      ),
      
      B = list(
        
        matrix(c(
          1.60, 0.20,
          1.20, 0.40,
          0.80, 0.80,
          0.40, 1.20,
          0.20, 1.60,
          0.10, 1.80
        ), ncol = 2, byrow = TRUE),
        
        matrix(c(
          1.50, 0.30,
          1.10, 0.50,
          0.70, 0.90,
          0.30, 1.30,
          0.10, 1.60,
          0.20, 1.70
        ), ncol = 2, byrow = TRUE),
        
        matrix(c(
          1.60,0.20,
          1.20,0.40,
          0.80,0.80,
          0.40,1.20,
          0.20,1.60,
          0.00,1.80
        ), ncol=2, byrow=TRUE)
        
      ),
      
      D = list(
        
        diag(rep(0.05, 6)),
        
        diag(rep(0.05, 6)),
        
        diag(rep(0.05, 6))
        
      ),
      
      nu = c(0.10,0.20,0.15),
      
      eta = c(5,8,6)
      
    ))
  }
  
  stop("Unknown scenario: ", name)
  
}
###############################################################
# 5. Data generation from true MCNFAC
###############################################################

simulate_mcnfac_data <- function(n, scen, censor_rate){
  
  g <- scen$g
  p <- scen$p
  q <- scen$q
  
  X <- make_X(n)
  
  #########################################################
  # True component labels
  #########################################################
  
  z <- sample(
    1:g,
    size = n,
    replace = TRUE,
    prob = scen$pi
  )
  
  #########################################################
  # Generate latent variables
  #########################################################
  
  ystar <- matrix(0, n, p)
  
  Vtrue <- integer(n)
  
  Wtrue <- numeric(n)
  
  Ftrue <- matrix(0, n, q)
  
  for(j in seq_len(n)){
    
    i <- z[j]
    
    Vtrue[j] <- rbinom(
      1,
      1,
      1 - scen$nu[i]
    )
    
    Wtrue[j] <- rmfac_W_from_V(
      Vtrue[j],
      scen$eta[i]
    )
    
    mu <- rmfac_mu_ij(
      X[j, ],
      scen$beta[[i]],
      p
    )
    
    Ftrue[j, ] <- rnorm(
      q,
      mean = 0,
      sd = sqrt(1 / Wtrue[j])
    )
    
    eps <- as.numeric(
      
      mvtnorm::rmvnorm(
        
        1,
        
        sigma =
          (1 / Wtrue[j]) *
          rmfac_sym(scen$D[[i]])
        
      )
      
    )
    
    ystar[j, ] <-
      
      as.numeric(
        
        mu +
          scen$B[[i]] %*% Ftrue[j, ] +
          eps
        
      )
    
  }
  
  #########################################################
  # Build censoring
  #
  # Apply the requested censoring rate to every response
  # variable, so censor_actual is approximately censor_rate.
  # Odd variables: left-censored
  # Even variables: right-censored
  #########################################################
  
  l_mat <- ystar
  u_mat <- ystar
  Yc <- ystar
  
  cens <- matrix(
    0L,
    nrow = n,
    ncol = p
  )
  
  cutoff <- rep(NA_real_, p)
  
  if(censor_rate > 0){
    
    if(censor_rate >= 0.5){
      stop("censor_rate must be smaller than 0.5.")
    }
    
    for(k in seq_len(p)){
      
      if(k %% 2 == 1){
        
        # Left censoring for odd-numbered responses
        cutoff[k] <- as.numeric(
          quantile(
            ystar[, k],
            probs = censor_rate,
            type = 7
          )
        )
        
        cens_idx <- which(
          ystar[, k] <= cutoff[k]
        )
        
        l_mat[cens_idx, k] <- -Inf
        u_mat[cens_idx, k] <- cutoff[k]
        Yc[cens_idx, k] <- cutoff[k]
        cens[cens_idx, k] <- 1L
        
      }else{
        
        # Right censoring for even-numbered responses
        cutoff[k] <- as.numeric(
          quantile(
            ystar[, k],
            probs = 1 - censor_rate,
            type = 7
          )
        )
        
        cens_idx <- which(
          ystar[, k] >= cutoff[k]
        )
        
        l_mat[cens_idx, k] <- cutoff[k]
        u_mat[cens_idx, k] <- Inf
        Yc[cens_idx, k] <- cutoff[k]
        cens[cens_idx, k] <- 2L
      }
    }
  }
  
  #########################################################
  # True covariance matrices
  #########################################################
  
  Sigma.true <-
    
    lapply(
      
      seq_len(g),
      
      function(i){
        
        scen$B[[i]] %*%
          t(scen$B[[i]]) +
          scen$D[[i]]
        
      }
      
    )
  
  
  #########################################################
  # Actual censoring percentage
  #########################################################
  
  censor_actual <-
    
    mean(
      cens != 0
    )
  
  #########################################################
  # Return everything needed by all studies
  #########################################################
  
  list(
    
    X = X,
    
    ystar = ystar,
    
    Yc = Yc,
    
    l_mat = l_mat,
    
    u_mat = u_mat,
    
    cens = cens,
    
    cutoff = cutoff,
    
    z = z,
    
    V = Vtrue,
    
    W = Wtrue,
    
    F = Ftrue,
    
    Sigma.true = Sigma.true,
    
    censor_actual = censor_actual
    
  )
  
}
###############################################################
# 6. Fitting MCNFAC and MFAC
###############################################################
fit_mcnfac_once <- function(
    l_mat,
    u_mat,
    X,
    g,
    q,
    n_init = 3,
    max_iter = MAX_ITER,
    tol = TOL
){
  
  n <- nrow(X)
  p <- ncol(l_mat)
  d <- ncol(X)
  
  best_fit <- NULL
  best_ll  <- -Inf
  
  for(init in seq_len(n_init)){
    
    params <- rmfac_initialize(
      l_mat,
      u_mat,
      X,
      g = g,
      q = q
    )
    
    ll_old <- -Inf
    ll_new <- NA_real_
    
    converged <- FALSE
    
    
    time_start <- proc.time()[["elapsed"]]
    
    for(iter in seq_len(max_iter)){
      
      ## ----------------------------
      ## E-step
      ## ----------------------------
      
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
      
      if(!is.null(fact$Fhat))
        estep$Fhat <- fact$Fhat
      
      
      ## ----------------------------
      ## Cycle 1
      ## ----------------------------
      
      params <- rmfac_CM1_update(
        params,
        estep,
        X
      )
      
      
      ## ----------------------------
      ## Cycle 2
      ## ----------------------------
      
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
      
      
      ## ----------------------------
      ## Parameter validity check
      ## ----------------------------
      
      rmfac_check_params(params)
      
      ## ----------------------------
      ## Observed log-likelihood
      ## ----------------------------
      
      ll_new <- rmfac_loglikelihood(
        l_mat,
        u_mat,
        X,
        params
      )
      
      
      if(!is.finite(ll_new)){
        
        warning(
          sprintf(
            "Initialization %d: non-finite log-likelihood.",
            init
          )
        )
        
        break
        
      }
      
      
      ## ----------------------------
      ## Convergence check
      ## ----------------------------
      
      delta <- ll_new - ll_old
      
      if(is.finite(ll_old)){
        
        if(delta < -1e-7){
          
          warning(
            sprintf(
              "Initialization %d: logLik decreased by %.6e at iteration %d.",
              init,
              delta,
              iter
            )
          )
          
        }
        
        rel_change <- abs(delta) / (abs(ll_old) + 1e-10)
        
        if(rel_change < tol){
          
          converged <- TRUE
          break
          
        }
        
      }
      
      ll_old <- ll_new
      
    }
    
    cpu_time <- proc.time()[["elapsed"]] - time_start
    
    if(!isTRUE(converged) || !is.finite(ll_new)){
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
    entropy <- NA_real_
    
    if(!is.null(estep_final$zhat)){
      
      entropy <- -mean(
        rowSums(
          estep_final$zhat *
            log(estep_final$zhat + 1e-12)
        )
      )
      
    }
    
    
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
      
      cpu_time = cpu_time,
      
      params = params,
      
      zhat = zhat,
      
      initialization = init
      
    )
    
    if(
      isTRUE(converged) &&
      is.finite(ll_new) &&
      ll_new > best_ll
    ){
      
      best_ll  <- ll_new
      best_fit <- current_fit
      
    }
    
  }
  
  if(is.null(best_fit)){
    
    warning("No initialization converged.")
    
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
        initialization = NA_integer_,
        params = NULL,
        zhat = NULL
      )
    )
  }
  
  return(best_fit)
  
}

fit_mfac_once <- function(
    l_mat,
    u_mat,
    X,
    g,
    q,
    n_init = 3,
    max_iter = MAX_ITER,
    tol = TOL
){
    
    n <- nrow(X)
    p <- ncol(l_mat)
    d <- ncol(X)
    
    best_fit <- NULL
    best_ll  <- -Inf
    
    for(init in seq_len(n_init)){
      
      #cat("MFAC initialization:", init, "\n")
      
      params <- rmfac_initialize(
        l_mat,
        u_mat,
        X,
        g = g,
        q = q
      )
      
      ## ----------------------------
      ## Force MFAC
      ## ----------------------------
      
      fixed_nu  <- rep(1e-6, g)
      fixed_eta <- rep(1 + 1e-4, g)
      
      params$nu  <- fixed_nu
      params$eta <- fixed_eta
      
      ll_old <- -Inf
      ll_new <- NA_real_
      
      converged <- FALSE
      
      
      time_start <- proc.time()[["elapsed"]]
      
      for(iter in seq_len(max_iter)){
        
        ## ----------------------------
        ## E-step
        ## ----------------------------
        
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
        
        if(!is.null(fact$Fhat))
          estep$Fhat <- fact$Fhat
        
        ## ----------------------------
        ## Update pi
        ## ----------------------------
        
        pi_new <- colMeans(estep$zhat)
        
        pi_new <- pmax(pi_new, 1e-12)
        
        params$pi <- pi_new / sum(pi_new)
        
        ## ----------------------------
        ## Keep contamination fixed
        ## ----------------------------
        
        params$nu  <- fixed_nu
        params$eta <- fixed_eta
        
        ## ----------------------------
        ## CM updates
        ## ----------------------------
        
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
        
        ## ----------------------------
        ## Enforce MFAC again
        ## ----------------------------
        
        params$nu  <- fixed_nu
        params$eta <- fixed_eta
        
        ## ----------------------------
        ## Parameter check
        ## ----------------------------
        
        rmfac_check_params(params)
        
        
        ## ----------------------------
        ## Log-likelihood
        ## ----------------------------
        
        ll_new <- rmfac_loglikelihood(
          l_mat,
          u_mat,
          X,
          params
        )
        
        
        if(!is.finite(ll_new)){
          
          warning("Non-finite log-likelihood encountered.")
          
          break
        }
        
        ## ----------------------------
        ## Convergence
        ## ----------------------------
        
        delta <- ll_new - ll_old
        
        if(is.finite(ll_old)){
          
          if(delta < -1e-7){
            
            warning(sprintf(
              "Observed log-likelihood decreased by %.6e at iteration %d.",
              delta,
              iter
            ))
            
          }
          
          if(abs(delta)/(abs(ll_old) + 1e-10) < tol){
            
            converged <- TRUE
            
            break
          }
          
        }
        
        ll_old <- ll_new
        
      }
      
      cpu_time <- proc.time()[["elapsed"]] - time_start
      
      if(!isTRUE(converged) || !is.finite(ll_new)){
        next
      }
      
      k <- count_params_mfac(
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
      
      zhat <- get_map_class(estep_final)
      
      entropy <- NA_real_
      
      if(!is.null(estep_final$zhat)){
        
        entropy <- -mean(
          rowSums(
            estep_final$zhat *
              log(estep_final$zhat + 1e-12)
          )
        )
        
      }
      
      
      current_fit <- list(
        
        model = "MFAC",
        
        g = g,
        
        q = q,
        
        logLik = ll_new,
        
        AIC = IC$AIC,
        
        BIC = IC$BIC,
        
        EDC = IC$EDC,
        
        entropy = entropy,
        
        k = k,
        
        iter = iter,
        
        converged = converged,
        
        cpu_time = cpu_time,
        
        initialization = init,
        
        params = params,
        
        zhat = zhat
        
      )
      
      if(
        isTRUE(converged) &&
        is.finite(ll_new) &&
        ll_new > best_ll
      ){
        
        best_ll <- ll_new
        best_fit <- current_fit
        
      }
      
    }
    
    if(is.null(best_fit)){
      
      warning("No initialization converged.")
      
      return(
        list(
          model = "MFAC",
          g = g,
          q = q,
          converged = FALSE,
          logLik = NA_real_,
          AIC = NA_real_,
          BIC = NA_real_,
          EDC = NA_real_,
          k = count_params_mfac(g, p, d, q),
          entropy = NA_real_,
          iter = NA_integer_,
          cpu_time = NA_real_,
          initialization = NA_integer_,
          params = NULL,
          zhat = NULL
        )
      )
    }
    
    return(best_fit)
    
  }
  
  fit_model_once <- function(
    model,
    l_mat,
    u_mat,
    X,
    g,
    q,
    n_init = 3
  ){
    
    if(model == "MCNFAC"){
      
      return(
        
        fit_mcnfac_once(
          l_mat = l_mat,
          u_mat = u_mat,
          X = X,
          g = g,
          q = q,
          n_init = n_init
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
  # 7. Recovery metrics
  ###############################################################
  
  recovery_metrics <- function(fit, scen){
    
    #------------------------------------------------------------
    # Return NA if fit failed
    #------------------------------------------------------------
    if(is.null(fit) || is.null(fit$params)){
      return(
        data.frame(
          beta_RMSE     = NA_real_,
          beta_abs_bias = NA_real_,
          Sigma_RMSE    = NA_real_,
          nu_RMSE       = NA_real_,
          eta_RMSE      = NA_real_,
          nu_abs_bias   = NA_real_,
          eta_abs_bias  = NA_real_
        )
      )
    }
    
    if(!isTRUE(fit$converged)){
      warning("Model did not converge.")
    }
    
    g <- scen$g
    
    perms <- gtools::permutations(g, g)
    
    best_score <- Inf
    best_perm  <- seq_len(g)
    
    #------------------------------------------------------------
    # Solve label switching
    #------------------------------------------------------------
    for(r in seq_len(nrow(perms))){
      
      pp <- perms[r, ]
      
      score <- 0
      valid <- TRUE
      
      for(i in seq_len(g)){
        
        S_hat <- sigma_from_params(
          fit$params,
          pp[i]
        )
        
        S_true <- rmfac_Sigma_i(
          scen$B[[i]],
          scen$D[[i]]
        )
        
        beta_hat <- fit$params$beta[[pp[i]]]
        beta_true <- scen$beta[[i]]
        
        if(
          any(!is.finite(S_hat)) ||
          any(!is.finite(beta_hat))
        ){
          valid <- FALSE
          break
        }
        
        beta_distance <- mean(
          (beta_hat - beta_true)^2,
          na.rm = TRUE
        )
        
        sigma_distance <- mean(
          (S_hat - S_true)^2,
          na.rm = TRUE
        )
        
        score <- score +
          beta_distance +
          sigma_distance
      }
      
      if(valid && score < best_score){
        best_score <- score
        best_perm <- pp
      }
      
    }
    
    #------------------------------------------------------------
    # Initialize accumulators
    #------------------------------------------------------------
    beta_mse  <- 0
    sigma_mse <- 0
    
    beta_abs  <- 0
    
    nu_mse  <- NA_real_
    eta_mse <- NA_real_
    
    nu_abs  <- NA_real_
    eta_abs <- NA_real_
    
    if(fit$model == "MCNFAC"){
      nu_mse  <- 0
      eta_mse <- 0
      nu_abs  <- 0
      eta_abs <- 0
    }
    
    #------------------------------------------------------------
    # Compute recovery measures
    #------------------------------------------------------------
    for(i in seq_len(g)){
      
      h <- best_perm[i]
      
      ## beta
      
      beta_diff <- fit$params$beta[[h]] - scen$beta[[i]]
      
      beta_mse <- beta_mse + mean(beta_diff^2, na.rm = TRUE)
      
      beta_abs <- beta_abs + mean(abs(beta_diff), na.rm = TRUE)
      
      ## covariance
      
      S_hat <- sigma_from_params(
        fit$params,
        h
      )
      
      S_true <- rmfac_Sigma_i(
        scen$B[[i]],
        scen$D[[i]]
      )
      
      sigma_mse <- sigma_mse +
        mean((S_hat - S_true)^2, na.rm = TRUE)
      
      ## contamination parameters
      
      if(fit$model == "MCNFAC"){
        
        nu_diff  <- fit$params$nu[h]  - scen$nu[i]
        eta_diff <- fit$params$eta[h] - scen$eta[i]
        
        nu_mse  <- nu_mse + nu_diff^2
        eta_mse <- eta_mse + eta_diff^2
        
        nu_abs  <- nu_abs + abs(nu_diff)
        eta_abs <- eta_abs + abs(eta_diff)
        
      }
      
    }
    
    #------------------------------------------------------------
    # Return metrics
    #------------------------------------------------------------
    data.frame(
      
      beta_RMSE =
        sqrt(beta_mse / g),
      
      beta_abs_bias =
        beta_abs / g,
      
      Sigma_RMSE =
        sqrt(sigma_mse / g),
      
      nu_RMSE =
        if(fit$model == "MCNFAC")
          sqrt(nu_mse / g)
      else
        NA_real_,
      
      eta_RMSE =
        if(fit$model == "MCNFAC")
          sqrt(eta_mse / g)
      else
        NA_real_,
      
      nu_abs_bias =
        if(fit$model == "MCNFAC")
          nu_abs / g
      else
        NA_real_,
      
      eta_abs_bias =
        if(fit$model == "MCNFAC")
          eta_abs / g
      else
        NA_real_
      
    )
    
  }
  
  ###############################################################
  # 8. Candidate-model helpers
  ###############################################################
  fit_candidate_models <- function(
    dat,
    model,
    fit_g,
    fit_q,
    n_init = 3
  ){
    
    fit_list <- list()
    
    for(gg in fit_g){
      
      for(qq in fit_q){
        
        key <- paste0(
          model,
          "_g", gg,
          "_q", qq
        )
        
        # cat(
        #   "   fitting ",
        #   key,
        #   " (n_init = ",
        #   n_init,
        #   ")\n",
        #   sep = ""
        # )
        
        fit_list[[key]] <- tryCatch(
          
          suppressWarnings(
            
            fit_model_once(
              
              model = model,
              
              l_mat = dat$l_mat,
              
              u_mat = dat$u_mat,
              
              X = dat$X,
              
              g = gg,
              
              q = qq,
              
              n_init = n_init
              
            )
            
          ),
          
          error = function(e){
            
            structure(
              
              list(
                
                message = conditionMessage(e),
                
                model = model,
                
                g = gg,
                
                q = qq,
                
                converged = FALSE
                
              ),
              
              class = "error"
              
            )
            
          }
          
        )
        
      }
      
    }
    
    return(fit_list)
    
  }
  
  ###############################################################
  # Choose best fitted model
  ###############################################################
  
  summarize_best_fit <- function(
    fit_list,
    dat,
    scen,
    criterion="BIC"){
    
    vals <- sapply(fit_list,function(fit){
      
      if(is.null(fit))
        return(Inf)
      
      if(inherits(fit,"error"))
        return(Inf)
      
      if(!isTRUE(fit$converged))
        return(Inf)
      
      if(!is.finite(fit$logLik))
        return(Inf)
      
      x <- fit[[criterion]]
      
      if(is.null(x) || !is.finite(x))
        return(Inf)
      
      x
      
    })
    
    if(all(!is.finite(vals))){
      
      return(
        
        data.frame(
          
          selected_g=NA_integer_,
          selected_q=NA_integer_,
          
          correct_g=NA_integer_,
          correct_q=NA_integer_,
          correct_selection=NA_integer_,
          
          best_AIC=NA_real_,
          best_BIC=NA_real_,
          best_EDC=NA_real_,
          best_logLik=NA_real_,
          
          best_iter=NA_integer_,
          best_converged=FALSE,
          best_cpu_time=NA_real_,
          best_initialization = NA_integer_,
          best_entropy = NA_real_,
          
          ARI=NA_real_,
          misclass=NA_real_,
          CCR=NA_real_
          
        )
        
      )
      
    }
    
    idx <- which.min(vals)
    
    best_fit <- fit_list[[idx]]
    
    ARI <- NA_real_
    MCR <- NA_real_
    CCR <- NA_real_
    
    if(!is.null(best_fit$zhat)){
      
      if(length(best_fit$zhat)==length(dat$z)){
        
        ARI <- tryCatch(
          
          mclust::adjustedRandIndex(
            dat$z,
            best_fit$zhat
          ),
          
          error=function(e) NA_real_
          
        )
        
        MCR <- tryCatch(
          
          misclass_rate(
            dat$z,
            best_fit$zhat
          ),
          
          error=function(e) NA_real_
          
        )
        
        if(is.finite(MCR))
          CCR <- 1-MCR
        
      }
      
    }
    
    data.frame(
      
      selected_g=best_fit$g,
      
      selected_q=best_fit$q,
      
      correct_g =
        as.integer(best_fit$g == scen$g),
      
      correct_q =
        as.integer(best_fit$q == scen$q),
      
      correct_selection =
        as.integer(
          best_fit$g == scen$g &&
            best_fit$q == scen$q
        ),
      
      best_AIC=best_fit$AIC,
      best_BIC=best_fit$BIC,
      best_EDC=best_fit$EDC,
      
      best_logLik=best_fit$logLik,
      
      best_iter=best_fit$iter,
      
      best_converged=best_fit$converged,
      
      best_cpu_time=best_fit$cpu_time,
      
      best_initialization = best_fit$initialization,
      
      best_entropy =
        if(!is.null(best_fit$entropy))
          best_fit$entropy
      else
        NA_real_,
      
      ARI=ARI,
      
      misclass=MCR,
      
      CCR=CCR
      
    )
    
  }
  
  
  ###############################################################
  # Save all candidate-model statistics
  ###############################################################
  
  add_candidate_columns <- function(out, fit_list){
    
    for(nm in names(fit_list)){
      
      fit <- fit_list[[nm]]
      
      if(is.null(fit) || inherits(fit,"error")){
        
        out[[paste0("AIC_",nm)]]     <- NA_real_
        out[[paste0("BIC_",nm)]]     <- NA_real_
        out[[paste0("EDC_",nm)]]     <- NA_real_
        out[[paste0("logLik_",nm)]]  <- NA_real_
        
        out[[paste0("iter_",nm)]]    <- NA_integer_
        out[[paste0("conv_",nm)]]    <- FALSE
        out[[paste0("time_",nm)]]    <- NA_real_
        
        out[[paste0("error_",nm)]] <-
          if(inherits(fit,"error"))
            fit$message
        else
          "NULL fit"
        
      }else{
        
        out[[paste0("AIC_",nm)]]    <- fit$AIC
        out[[paste0("BIC_",nm)]]    <- fit$BIC
        out[[paste0("EDC_",nm)]]    <- fit$EDC
        out[[paste0("logLik_",nm)]] <- fit$logLik
        
        out[[paste0("iter_",nm)]]   <- fit$iter
        out[[paste0("conv_",nm)]]   <- fit$converged
        out[[paste0("time_",nm)]]   <- fit$cpu_time
        
        out[[paste0("error_",nm)]] <- ""
        
      }
      
    }
    
    out
    
  }
  
  ###############################################################
  # 9. Study 1: Parameter recovery under censoring
  ###############################################################
  
  run_study1_task <- function(task){
    set.seed(
      100000 + task$rep +
        task$n*10 +
        round(task$censor*100)
    )
    
    scen <- make_scenario("moderate_g2_q1")
    
    dat <- simulate_mcnfac_data(
      n = task$n,
      scen = scen,
      censor_rate = task$censor
    )
    
    models <- c("MFAC","MCNFAC")
    
    cat(
      "\n[Study1] n =",task$n,
      "| censor =",task$censor,
      "| rep =",task$rep,"\n"
    )
    
    rows <- list()
    counter <- 1
    
    for(model in models){
      
      cat("   model =",model,"\n")
      
      fit <- tryCatch(
        
        suppressWarnings(
          
          fit_model_once(
            model,
            dat$l_mat,
            dat$u_mat,
            dat$X,
            scen$g,
            scen$q
          )
          
        ),
        
        error=function(e){
          
          structure(
            list(message=conditionMessage(e)),
            class="error"
          )
          
        }
        
      )
      
      base <- data.frame(
        
        study="Study1_recovery",
        
        model=model,
        
        scenario=scen$name,
        
        rep=task$rep,
        
        n=task$n,
        
        censor_target=task$censor,
        
        censor_actual=dat$censor_actual,
        
        true_g=scen$g,
        
        true_q=scen$q,
        
        stringsAsFactors=FALSE
        
      )
      
      ###########################################################
      ## fitting failed
      ###########################################################
      
      if(inherits(fit, "error")){
        
        out <- cbind(
          
          base,
          
          data.frame(
            
            converged = FALSE,
            
            iter = NA_integer_,
            
            logLik = NA_real_,
            
            AIC = NA_real_,
            
            BIC = NA_real_,
            
            EDC = NA_real_,
            
            entropy = NA_real_,
            
            cpu_time = NA_real_,
            
            ARI = NA_real_,
            
            misclass = NA_real_,
            
            CCR = NA_real_,
            
            beta_RMSE = NA_real_,
            
            beta_abs_bias = NA_real_,
            
            Sigma_RMSE = NA_real_,
            
            nu_RMSE = NA_real_,
            
            eta_RMSE = NA_real_,
            
            nu_abs_bias = NA_real_,
            
            eta_abs_bias = NA_real_,
            
            error_message = fit$message,
            
            stringsAsFactors = FALSE
            
          )
          
        )
        
      }else{
        
        #########################################################
        ## clustering measures
        #########################################################
        
        ARI <- NA_real_
        MCR <- NA_real_
        CCR <- NA_real_
        
        if(!is.null(fit$zhat)){
          
          if(length(fit$zhat)==length(dat$z)){
            
            ARI <- tryCatch(
              
              mclust::adjustedRandIndex(
                dat$z,
                fit$zhat
              ),
              
              error=function(e) NA_real_
              
            )
            
            MCR <- tryCatch(
              
              misclass_rate(
                dat$z,
                fit$zhat
              ),
              
              error=function(e) NA_real_
              
            )
            
            if(is.finite(MCR))
              CCR <- 1-MCR
            
          }
          
        }
        
        #########################################################
        ## parameter recovery
        #########################################################
        
        if(
          !is.null(fit$params) &&
          isTRUE(fit$converged)
        ){
          
          rec <- tryCatch(
            
            recovery_metrics(
              fit,
              scen
            ),
            
            error=function(e){
              
              data.frame(
                
                beta_RMSE=NA_real_,
                
                beta_abs_bias=NA_real_,
                
                Sigma_RMSE=NA_real_,
                
                nu_RMSE=NA_real_,
                
                eta_RMSE=NA_real_,
                
                nu_abs_bias=NA_real_,
                
                eta_abs_bias=NA_real_
                
              )
              
            }
            
          )
          
        }else{
          
          rec <- data.frame(
            
            beta_RMSE=NA_real_,
            
            beta_abs_bias=NA_real_,
            
            Sigma_RMSE=NA_real_,
            
            nu_RMSE=NA_real_,
            
            eta_RMSE=NA_real_,
            
            nu_abs_bias=NA_real_,
            
            eta_abs_bias=NA_real_
            
          )
          
        }
        
        out <- cbind(
          
          base,
          
          data.frame(
            
            converged = fit$converged,
            
            initialization = fit$initialization,
            
            iter = fit$iter,
            
            logLik = fit$logLik,
            
            AIC = fit$AIC,
            
            BIC = fit$BIC,
            
            EDC = fit$EDC,
            
            entropy = fit$entropy,
            
            cpu_time = fit$cpu_time,
            
            ARI = ARI,
            
            misclass = MCR,
            
            CCR = CCR,
            
            error_message = "",
            
            stringsAsFactors = FALSE
            
          ),
          
          rec
          
        )
        
      }
      
      rows[[counter]] <- out
      
      counter <- counter + 1
      
    }
    
    bind_rows_fill(rows)
    
  }
  
  
  ###############################################################
  # Run Study 1
  ###############################################################
  
  run_study1_recovery <- function(){
    
    n_grid <- c(
      300,
      500,
      1000
    )
    
    censor_grid <- c(
      0.00,
      0.10,
      0.20,
      0.30
    )
    
    tasks <- vector(
      "list",
      length(n_grid)*
        length(censor_grid)*
        R_sim
    )
    
    counter <- 1
    
    for(n_val in n_grid){
      
      for(cens_val in censor_grid){
        
        for(r in seq_len(R_sim)){
          
          tasks[[counter]] <- list(
            
            n=n_val,
            
            censor=cens_val,
            
            rep=r
            
          )
          
          counter <- counter + 1
          
        }
        
      }
      
    }
    
    progress_file <-
      file.path(
        out_dir,
        "Study1_recovery_PROGRESS.csv"
      )
    
    res <- process_tasks(
      
      tasks,
      
      run_study1_task,
      
      progress_file
      
    )
    
    write.csv(
      
      res,
      
      file.path(
        out_dir,
        "Study1_recovery_raw.csv"
      ),
      
      row.names=FALSE
      
    )
    
    cat(
      "\nStudy 1 completed.",
      nrow(res),
      "rows saved.",
      sum(res$converged, na.rm = TRUE),
      "fits converged.\n"
    )
    
    res
    
  }
  
  
  ###############################################################
  # 10. Study 2: Recovery of true number of components g
  ###############################################################
  
  run_study2_task <- function(task){
    
    ## reproducibility
    set.seed(
      100000 +
        task$rep +
        task$n*10 +
        round(task$censor*100)
    )
    
    scen <- make_scenario("moderate_g2_q1")
    
    fit_g <- 1:4
    fit_q <- 1
    
    models <- c("MFAC","MCNFAC")
    
    cat(
      "\n[Study 2] n =",
      task$n,
      "| censor =",
      task$censor,
      "| rep",
      task$rep,
      "\n"
    )
    
    ##----------------------------------------------------------
    ## Data generation
    ##----------------------------------------------------------
    dat <- simulate_mcnfac_data(
      n = task$n,
      scen = scen,
      censor_rate = task$censor
    )
    
    rows <- vector("list", length(models))
    
    for(m in seq_along(models)){
      
      model <- models[m]
      
      cat(" model =", model, "\n")
      
      ##----------------------------------------------------------
      ## Fit all candidate models
      ##----------------------------------------------------------
      fit_list <- fit_candidate_models(
        dat,
        model,
        fit_g,
        fit_q
      )
      
      ##----------------------------------------------------------
      ## Select best model
      ##----------------------------------------------------------
      best <- tryCatch(
        
        summarize_best_fit(
          fit_list,
          dat,
          scen,
          criterion="BIC"
        ),
        
        error=function(e){
          
          data.frame(
            selected_g      = NA_integer_,
            selected_q      = NA_integer_,
            correct_g       = NA_integer_,
            correct_q       = NA_integer_,
            correct_selection = NA_integer_,
            best_AIC        = NA_real_,
            best_BIC        = NA_real_,
            best_EDC        = NA_real_,
            best_logLik     = NA_real_,
            best_iter       = NA_integer_,
            best_converged  = FALSE,
            best_cpu_time   = NA_real_,
            best_initialization = NA_integer_,
            best_entropy    = NA_real_,
            
            ARI             = NA_real_,
            misclass        = NA_real_,
            CCR             = NA_real_
          )
          
        }
        
      )
      
      ##----------------------------------------------------------
      ## Base output
      ##----------------------------------------------------------
      out <- data.frame(
        
        study="Study2_component_selection",
        
        model=model,
        
        scenario=scen$name,
        
        rep=task$rep,
        
        n=task$n,
        
        censor_target=task$censor,
        
        censor_actual=dat$censor_actual,
        
        true_g=scen$g,
        
        true_q=scen$q,
        
        stringsAsFactors=FALSE
        
      )
      
      out <- cbind(out, best)
      
      ##----------------------------------------------------------
      ## Store all candidate-model statistics
      ##----------------------------------------------------------
      out <- tryCatch(
        
        add_candidate_columns(out, fit_list),
        
        error=function(e) out
        
      )
      
      rows[[m]] <- out
      
    }
    
    bind_rows_fill(rows)
    
  }
  
  
  ###############################################################
  # Main function
  ###############################################################
  
  run_study2_component_selection <- function(){
    
    n_grid <- c(300, 1000)
    
    censor_grid <- c(
      0.00,
      0.10,
      0.20,
      0.30
    )
    
    tasks <- list()
    
    counter <- 1
    
    for(n_val in n_grid){
      
      for(cens_val in censor_grid){
        
        for(r in seq_len(R_sim)){
          
          tasks[[counter]] <- list(
            
            n       = n_val,
            
            censor  = cens_val,
            
            rep     = r
            
          )
          
          counter <- counter + 1
          
        }
        
      }
      
    }
    
    res <- process_tasks(
      
      tasks,
      
      run_study2_task,
      
      file.path(
        out_dir,
        "Study2_component_selection_PROGRESS.csv"
      )
      
    )
    
    write.csv(
      
      res,
      
      file.path(
        out_dir,
        "Study2_component_selection_raw.csv"
      ),
      
      row.names=FALSE
      
    )
    
    invisible(res)
    
  }
  
  ###############################################################
  # 11. Study 3: Clustering performance
  ###############################################################
  
  run_study3_task <- function(task){
    
    ##----------------------------------------------------------
    ## Reproducibility
    ##----------------------------------------------------------
    set.seed(
      100000 +
        task$rep +
        task$n*10 +
        round(task$censor*100)
    )
    
    scen <- make_scenario("moderate_g2_q1")
    
    fit_g <- scen$g
    fit_q <- scen$q
    
    models <- c("MFAC","MCNFAC")
    
    cat(
      "\n[Study 3] n =",
      task$n,
      "| censor =",
      task$censor,
      "| rep",
      task$rep,
      "\n"
    )
    
    ##----------------------------------------------------------
    ## Generate data
    ##----------------------------------------------------------
    dat <- simulate_mcnfac_data(
      n = task$n,
      scen = scen,
      censor_rate = task$censor
    )
    
    rows <- vector("list", length(models))
    
    for(m in seq_along(models)){
      
      model <- models[m]
      
      cat(" model =", model, "\n")
      
      ##----------------------------------------------------------
      ## Fit candidate model(s)
      ##----------------------------------------------------------
      fit_list <- fit_candidate_models(
        dat,
        model,
        fit_g,
        fit_q
      )
      
      ##----------------------------------------------------------
      ## Extract fitted model results
      ##----------------------------------------------------------
      best <- tryCatch(
        
        summarize_best_fit(
          fit_list,
          dat,
          scen,
          criterion = "BIC"
        ),
        
        error = function(e){
          
          data.frame(
            
            selected_g = NA_integer_,
            selected_q = NA_integer_,
            
            correct_g = NA_integer_,
            correct_q = NA_integer_,
            correct_selection = NA_integer_,
            
            best_AIC = NA_real_,
            best_BIC = NA_real_,
            best_EDC = NA_real_,
            best_logLik = NA_real_,
            
            best_iter = NA_integer_,
            best_converged = FALSE,
            best_cpu_time = NA_real_,
            best_initialization = NA_integer_,
            best_entropy = NA_real_,
            
            ARI = NA_real_,
            misclass = NA_real_,
            CCR = NA_real_,
            
            stringsAsFactors = FALSE
            
          )
          
        }
        
      )
      
      ##----------------------------------------------------------
      ## Base output
      ##----------------------------------------------------------
      out <- data.frame(
        
        study = "Study3_clustering",
        
        model = model,
        
        scenario = scen$name,
        
        rep = task$rep,
        
        n = task$n,
        
        censor_target = task$censor,
        
        censor_actual = dat$censor_actual,
        
        true_g = scen$g,
        
        true_q = scen$q,
        
        stringsAsFactors = FALSE
        
      )
      
      out <- cbind(out, best)
      
      ##----------------------------------------------------------
      ## Save candidate-model information
      ##----------------------------------------------------------
      out <- tryCatch(
        
        add_candidate_columns(out, fit_list),
        
        error = function(e) out
        
      )
      
      rows[[m]] <- out
      
    }
    
    bind_rows_fill(rows)
    
  }
  
  
  ###############################################################
  # Main Study 3
  ###############################################################
  
  run_study3_clustering <- function(){
    
    n_grid <- c(
      300,
      500,
      1000
    )
    
    censor_grid <- c(
      0.00,
      0.10,
      0.20,
      0.30
    )
    
    tasks <- list()
    
    counter <- 1
    
    for(n_val in n_grid){
      
      for(cens_val in censor_grid){
        
        for(r in seq_len(R_sim)){
          
          tasks[[counter]] <- list(
            
            n = n_val,
            
            censor = cens_val,
            
            rep = r
            
          )
          
          counter <- counter + 1
          
        }
        
      }
      
    }
    
    res <- process_tasks(
      
      tasks,
      
      run_study3_task,
      
      file.path(
        out_dir,
        "Study3_clustering_PROGRESS.csv"
      )
      
    )
    
    write.csv(
      
      res,
      
      file.path(
        out_dir,
        "Study3_clustering_raw.csv"
      ),
      
      row.names = FALSE
      
    )
    
    invisible(res)
    
  }
  
  ###############################################################
  # 12. Study 4: Complex latent structure
  # Recover both g and q
  ###############################################################
  
  run_study4_task <- function(task){
    
    ##----------------------------------------------------------
    ## Reproducibility
    ##----------------------------------------------------------
    set.seed(
      300000 +
        task$rep +
        task$n*10 +
        round(task$censor*100)
    )
    
    scen <- make_scenario("complex_g3_q2")
    
    fit_g <- 1:5
    fit_q <- 1:3
    
    models <- c("MFAC","MCNFAC")
    
    cat(
      "\n[Study 4] n =",
      task$n,
      "| censor =",
      task$censor,
      "| rep",
      task$rep,
      "\n"
    )
    
    ##----------------------------------------------------------
    ## Generate data
    ##----------------------------------------------------------
    dat <- simulate_mcnfac_data(
      n = task$n,
      scen = scen,
      censor_rate = task$censor
    )
    
    rows <- vector("list", length(models))
    
    for(m in seq_along(models)){
      
      model <- models[m]
      
      cat(" model =", model, "\n")
      
      ##----------------------------------------------------------
      ## Fit all candidate models
      ##----------------------------------------------------------
      fit_list <- tryCatch(
        
        fit_candidate_models(
          dat,
          model,
          fit_g,
          fit_q
        ),
        
        error=function(e){
          
          structure(
            list(message=conditionMessage(e)),
            class="error"
          )
          
        }
      )
      
      ##----------------------------------------------------------
      ## Select best model
      ##----------------------------------------------------------
      best <- tryCatch(
        
        summarize_best_fit(
          fit_list,
          dat,
          scen,
          criterion = "BIC"
        ),
        
        error = function(e){
          
          data.frame(
            
            selected_g = NA_integer_,
            selected_q = NA_integer_,
            
            correct_g = NA_integer_,
            correct_q = NA_integer_,
            correct_selection = NA_integer_,
            
            best_AIC = NA_real_,
            best_BIC = NA_real_,
            best_EDC = NA_real_,
            best_logLik = NA_real_,
            
            best_iter = NA_integer_,
            best_converged = FALSE,
            best_cpu_time = NA_real_,
            best_initialization = NA_integer_,
            best_entropy = NA_real_,
            
            ARI = NA_real_,
            misclass = NA_real_,
            CCR = NA_real_,
            
            stringsAsFactors = FALSE
            
          )
          
        }
        
      )
      
      ##----------------------------------------------------------
      ## Recovery for TRUE model only
      ##----------------------------------------------------------
      true_key <- paste0(
        model,
        "_g",
        scen$g,
        "_q",
        scen$q
      )
      
      rec <- tryCatch({
        
        true_fit <- fit_list[[true_key]]
        
        if(!is.null(true_fit) &&
           !inherits(true_fit,"error") &&
           isTRUE(true_fit$converged)){
          
          recovery_metrics(
            true_fit,
            scen
          )
          
        }else{
          
          data.frame(
            
            beta_RMSE = NA_real_,
            beta_abs_bias = NA_real_,
            Sigma_RMSE = NA_real_,
            
            nu_RMSE = NA_real_,
            eta_RMSE = NA_real_,
            
            nu_abs_bias = NA_real_,
            eta_abs_bias = NA_real_,
            
            stringsAsFactors = FALSE
            
          )
          
        }
        
      },
      
      error=function(e){
        
        data.frame(
          
          beta_RMSE = NA_real_,
          beta_abs_bias = NA_real_,
          Sigma_RMSE = NA_real_,
          
          nu_RMSE = NA_real_,
          eta_RMSE = NA_real_,
          
          nu_abs_bias = NA_real_,
          eta_abs_bias = NA_real_,
          
          stringsAsFactors = FALSE
          
        )
        
      })
      
      ##----------------------------------------------------------
      ## Base output
      ##----------------------------------------------------------
      out <- data.frame(
        
        study = "Study4_complex_gq",
        
        model = model,
        
        scenario = scen$name,
        
        rep = task$rep,
        
        n = task$n,
        
        censor_target = task$censor,
        
        censor_actual = dat$censor_actual,
        
        true_g = scen$g,
        
        true_q = scen$q,
        
        stringsAsFactors = FALSE
        
      )
      
      out <- cbind(
        
        out,
        
        best,
        
        rec
        
      )
      
      ##----------------------------------------------------------
      ## Store all candidate-model information
      ##----------------------------------------------------------
      out <- tryCatch(
        
        add_candidate_columns(
          out,
          fit_list
        ),
        
        error=function(e) out
        
      )
      
      rows[[m]] <- out
      
    }
    
    bind_rows_fill(rows)
    
  }
  
  
  ###############################################################
  # Main Study 4
  ###############################################################
  
  run_study4_complex_gq <- function(){
    
    n_grid <- c(
      500
    )
    
    censor_grid <- c(
      0.10,
      0.30
    )
    
    tasks <- list()
    
    counter <- 1
    
    for(n_val in n_grid){
      
      for(cens_val in censor_grid){
        
        for(r in seq_len(R_sim)){
          
          tasks[[counter]] <- list(
            
            n = n_val,
            
            censor = cens_val,
            
            rep = r
            
          )
          
          counter <- counter + 1
          
        }
        
      }
      
    }
    
    res <- process_tasks(
      
      tasks,
      
      run_study4_task,
      
      file.path(
        
        out_dir,
        
        "Study4_complex_gq_PROGRESS.csv"
        
      )
      
    )
    
    write.csv(
      
      res,
      
      file.path(
        
        out_dir,
        
        "Study4_complex_gq_raw.csv"
        
      ),
      
      row.names = FALSE
      
    )
    
    invisible(res)
    
  }
  ###############################################################
  # 12b. Study 5: Heavy-contamination robustness stress test
  ###############################################################
  
  run_study5_task <- function(task){
    
    set.seed(
      400000 +
        task$rep +
        task$n*10 +
        round(task$censor*100)
    )
    scen <- make_scenario("heavy_g2_q1")
    
    fit_g <- 1:5
    fit_q <- 1
    
    models <- c("MFAC","MCNFAC")
    
    cat(
      "\n[Study 5] n =", task$n,
      "| censor =", task$censor,
      "| rep", task$rep, "\n"
    )
    
    dat <- tryCatch(
      simulate_mcnfac_data(task$n, scen, task$censor),
      error = function(e) e
    )
    
    if(inherits(dat,"error")){
      
      out <- data.frame(
        study="Study5_heavy_contamination",
        model=NA_character_,
        scenario=scen$name,
        rep=task$rep,
        n=task$n,
        censor_target=task$censor,
        censor_actual=NA_real_,
        true_g=scen$g,
        true_q=scen$q,
        error=conditionMessage(dat),
        stringsAsFactors=FALSE
      )
      
      return(out)
      
    }
    
    rows <- vector("list", length(models))
    
    for(mm in seq_along(models)){
      
      model <- models[mm]
      
      cat(" model =", model, "\n")
      
      ###########################################################
      ## Fit all candidate models
      ###########################################################
      
      fit_list <- tryCatch(
        
        fit_candidate_models(
          dat,
          model,
          fit_g,
          fit_q
        ),
        
        error=function(e) e
        
      )
      
      base <- data.frame(
        
        study="Study5_heavy_contamination",
        model=model,
        scenario=scen$name,
        rep=task$rep,
        n=task$n,
        censor_target=task$censor,
        censor_actual=dat$censor_actual,
        true_g=scen$g,
        true_q=scen$q,
        stringsAsFactors=FALSE
        
      )
      
      ###########################################################
      ## Entire fitting failed
      ###########################################################
      
      if(inherits(fit_list,"error")){
        
        rows[[mm]] <- cbind(
          
          base,
          
          data.frame(
            
            selected_g = NA_integer_,
            selected_q = NA_integer_,
            
            correct_g = NA_integer_,
            correct_q = NA_integer_,
            correct_selection = NA_integer_,
            
            best_AIC = NA_real_,
            best_BIC = NA_real_,
            best_EDC = NA_real_,
            best_logLik = NA_real_,
            
            best_iter = NA_integer_,
            best_converged = FALSE,
            best_cpu_time = NA_real_,
            best_initialization = NA_integer_,
            best_entropy = NA_real_,
            
            ARI = NA_real_,
            misclass = NA_real_,
            CCR = NA_real_,
            
            error = conditionMessage(fit_list),
            
            stringsAsFactors = FALSE
            
          )
          
        )
        
        next
        
      }
      
      ###########################################################
      ## Best candidate model
      ###########################################################
      
      best <- tryCatch(
        
        summarize_best_fit(
          fit_list,
          dat,
          scen,
          criterion="BIC"
        ),
        
        error=function(e){
          
          data.frame(
            
            selected_g = NA_integer_,
            selected_q = NA_integer_,
            
            correct_g = NA_integer_,
            correct_q = NA_integer_,
            correct_selection = NA_integer_,
            
            best_AIC = NA_real_,
            best_BIC = NA_real_,
            best_EDC = NA_real_,
            best_logLik = NA_real_,
            
            best_iter = NA_integer_,
            best_converged = FALSE,
            best_cpu_time = NA_real_,
            best_initialization = NA_integer_,
            best_entropy = NA_real_,
            
            ARI = NA_real_,
            misclass = NA_real_,
            CCR = NA_real_,
            
            stringsAsFactors = FALSE
            
          )
          
        }
        
      )
      
      ###########################################################
      ## Parameter recovery
      ###########################################################
      
      true_key <- paste0(
        model,
        "_g",
        scen$g,
        "_q",
        scen$q
      )
      
      true_fit <- fit_list[[true_key]]
      
      rec <- tryCatch(
        
        {
          
          if(!is.null(true_fit) &&
             !inherits(true_fit,"error") &&
             isTRUE(true_fit$converged)){
            
            recovery_metrics(true_fit, scen)
            
          }else{
            
            data.frame(
              
              beta_RMSE=NA_real_,
              beta_abs_bias=NA_real_,
              Sigma_RMSE=NA_real_,
              nu_RMSE=NA_real_,
              eta_RMSE=NA_real_,
              nu_abs_bias=NA_real_,
              eta_abs_bias=NA_real_
              
            )
            
          }
          
        },
        
        error=function(e){
          
          data.frame(
            
            beta_RMSE=NA_real_,
            beta_abs_bias=NA_real_,
            Sigma_RMSE=NA_real_,
            nu_RMSE=NA_real_,
            eta_RMSE=NA_real_,
            nu_abs_bias=NA_real_,
            eta_abs_bias=NA_real_
            
          )
          
        }
        
      )
      
      ###########################################################
      ## Final row
      ###########################################################
      
      out <- cbind(
        base,
        best,
        rec
      )
      
      out$error <- ""
      
      out <- add_candidate_columns(
        out,
        fit_list
      )
      
      rows[[mm]] <- out
      
    }
    
    bind_rows_fill(rows)
    
  }
  
  run_study5_heavy_contamination <- function(){
    
    n_grid <- c(500)
    
    censor_grid <- c(
      0.00,
      0.10,
      0.20,
      0.30
    )
    
    tasks <- list()
    
    counter <- 1
    
    for(n_val in n_grid){
      
      for(cens_val in censor_grid){
        
        for(r in seq_len(R_sim)){
          
          tasks[[counter]] <- list(
            
            n=n_val,
            censor=cens_val,
            rep=r
            
          )
          
          counter <- counter + 1
          
        }
        
      }
      
    }
    
    res <- process_tasks(
      
      tasks,
      run_study5_task,
      file.path(
        out_dir,
        "Study5_heavy_contamination_PROGRESS.csv"
      )
      
    )
    
    write.csv(
      
      res,
      
      file.path(
        out_dir,
        "Study5_heavy_contamination_raw.csv"
      ),
      
      row.names=FALSE
      
    )
    
    res
    
  }
  ###############################################################
  # 13. Run selected studies
  ###############################################################
  
  CL <- NULL
  
  ###############################################################
  ## Create parallel cluster
  ###############################################################
  
  if(isTRUE(USE_PARALLEL) && N_CORES > 1){
    
    cat("\nCreating parallel cluster with", N_CORES, "workers ...\n")
    
    CL <- parallel::makeCluster(N_CORES)
    
    parallel::clusterSetRNGStream(CL, 123)
    
    parallel::clusterEvalQ(CL, {
      
      source("MCNFAC_Functions.R")
      
      suppressPackageStartupMessages({
        
        library(mvtnorm)
        library(mclust)
        library(gtools)
        
      })
      
      NULL
      
    })
    
    parallel::clusterExport(
      CL,
      varlist = setdiff(ls(envir = .GlobalEnv), "CL"),
      envir = .GlobalEnv
    )
    
  }
  
  ###############################################################
  
  
  ###############################################################
  ## Helper to run one study safely
  ###############################################################
  
  run_one_study <- function(study_name, study_fun){
    
    if(!(study_name %in% studies_to_run))
      return(NULL)
    
    cat("\n=====================================================\n")
    cat("Running", study_name, "\n")
    cat("=====================================================\n")
    
    out <- tryCatch(
      
      study_fun(),
      
      error = function(e){
        
        warning(
          paste(
            "Study failed:",
            study_name,
            "\n",
            conditionMessage(e)
          )
        )
        
        NULL
        
      }
      
    )
    
    gc()
    
    out
    
  }
  
  ###############################################################
  ## Run requested studies
  ###############################################################
  
  all_results <- list()
  
  all_results[["Study1_recovery"]] <-
    run_one_study(
      "Study1_recovery",
      run_study1_recovery
    )
  
  all_results[["Study2_component_selection"]] <-
    run_one_study(
      "Study2_component_selection",
      run_study2_component_selection
    )
  
  all_results[["Study3_clustering"]] <-
    run_one_study(
      "Study3_clustering",
      run_study3_clustering
    )
  
  all_results[["Study4_complex_gq"]] <-
    run_one_study(
      "Study4_complex_gq",
      run_study4_complex_gq
    )
  
  all_results[["Study5_heavy_contamination"]] <-
    run_one_study(
      "Study5_heavy_contamination",
      run_study5_heavy_contamination
    )
  
  
  ###############################################################
  ## Stop parallel cluster
  ###############################################################
  
  if(!is.null(CL)){
    
    cat("\nStopping parallel cluster ...\n")
    
    parallel::stopCluster(CL)
    
    CL <- NULL
  }
  ###############################################################
  ## Combine all results
  ###############################################################
  
  sim_results <- bind_rows_fill(all_results)
  
  write.csv(
    
    sim_results,
    
    file.path(
      out_dir,
      "MCNFAC_Lin_FINAL_all_raw_results.csv"
    ),
    
    row.names = FALSE
    
  )
  
  cat(
    "\nSimulation finished successfully.\n",
    "Total rows =", nrow(sim_results), "\n"
  )
  
  gc()
  
  ###############################################################
  # 14. Summary tables
  ###############################################################
  
  cat("\n=========================================\n")
  cat("Creating summary tables ...\n")
  cat("=========================================\n")
  
  if(nrow(sim_results)==0){
    
    warning("No simulation results available.")
    
  }else{
    
    #############################################################
    ## convert logical columns
    #############################################################
    
    logical_cols <- intersect(
      
      c("converged","best_converged"),
      
      names(sim_results)
      
    )
    
    for(v in logical_cols)
      sim_results[[v]] <- as.numeric(sim_results[[v]])
    
    #############################################################
    ## helper
    #############################################################
    
    make_summary <- function(vars){
      
      vars <- intersect(vars,names(sim_results))
      
      if(length(vars)==0)
        return(NULL)
      
      aggregate(
        
        sim_results[,vars,drop=FALSE],
        
        by=list(
          
          study=sim_results$study,
          
          model=sim_results$model,
          
          n=sim_results$n,
          
          censor_target=sim_results$censor_target
          
        ),
        
        FUN=safe_mean
        
      )
      
    }
    
    #############################################################
    ## Recovery summary
    #############################################################
    
    cat("\n===== SUMMARY: Parameter recovery =====\n")
    
    recovery_summary <- make_summary(
      
      c(
        
        "beta_RMSE",
        
        "beta_abs_bias",
        
        "Sigma_RMSE",
        
        "nu_RMSE",
        
        "eta_RMSE",
        
        "nu_abs_bias",
        
        "eta_abs_bias",
        
        "entropy",
        
        "cpu_time",
        
        "converged"
        
      )
      
    )
    
    if(!is.null(recovery_summary)){
      
      print(recovery_summary)
      
      write.csv(
        
        recovery_summary,
        
        file.path(
          
          out_dir,
          
          "MCNFAC_Lin_FINAL_recovery_summary.csv"
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    #############################################################
    ## Selection summary
    #############################################################
    
    cat("\n===== SUMMARY: Component selection =====\n")
    
    selection_summary <- make_summary(
      c(
        "correct_g",
        "correct_q",
        "correct_selection",
        "best_converged",
        "best_entropy",
        "best_cpu_time",
        "best_AIC",
        "best_BIC",
        "best_EDC"
      )
    )
    
    if(!is.null(selection_summary)){
      
      print(selection_summary)
      
      write.csv(
        
        selection_summary,
        
        file.path(
          
          out_dir,
          
          "MCNFAC_Lin_FINAL_selection_summary.csv"
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    #############################################################
    ## Clustering summary
    #############################################################
    
    cat("\n===== SUMMARY: Clustering =====\n")
    
    cluster_summary <- make_summary(
      
      c(
        
        "ARI",
        
        "misclass",
        
        "CCR",
        
        "best_entropy",
        
        "best_cpu_time",
        
        "cpu_time"
        
      )
      
    )
    
    if(!is.null(cluster_summary)){
      
      print(cluster_summary)
      
      write.csv(
        
        cluster_summary,
        
        file.path(
          
          out_dir,
          
          "MCNFAC_Lin_FINAL_clustering_summary.csv"
          
        ),
        
        row.names=FALSE
        
      )
      
    }
    
    #############################################################
    ## Selected model frequency
    #############################################################
    
    if(all(
      
      c(
        
        "study",
        
        "model",
        
        "selected_g",
        
        "selected_q"
        
      ) %in% names(sim_results)
      
    )){
      
      cat("\n===== Selected model frequency =====\n")
      
      try(
        
        print(
          
          table(
            
            sim_results$study,
            
            sim_results$model,
            
            sim_results$selected_g,
            
            sim_results$selected_q,
            
            useNA="ifany"
            
          )
          
        ),
        
        silent=TRUE
        
      )
      
    }
    
    #############################################################
    ## Diagnostics
    #############################################################
    
    failed_runs <- 0
    
    if("error" %in% names(sim_results)){
      
      failed_runs <- failed_runs +
        sum(
          !is.na(sim_results$error) &
            sim_results$error != ""
        )
      
    }
    
    if("error_message" %in% names(sim_results)){
      
      failed_runs <- failed_runs +
        sum(
          !is.na(sim_results$error_message) &
            sim_results$error_message != ""
        )
      
    }
    
    cat(
      "\nNumber of failed runs =",
      failed_runs,
      "\n"
    )
    
    conv_vec <- c()
    
    if("converged" %in% names(sim_results)){
      
      conv_vec <- c(
        conv_vec,
        sim_results$converged
      )
      
    }
    
    if("best_converged" %in% names(sim_results)){
      
      conv_vec <- c(
        conv_vec,
        sim_results$best_converged
      )
      
    }
    
    if("converged" %in% names(sim_results)){
      cat(
        "Recovery convergence =",
        round(mean(sim_results$converged, na.rm=TRUE),4),
        "\n"
      )
    }
    
    if("best_converged" %in% names(sim_results)){
      cat(
        "Best-model convergence =",
        round(mean(sim_results$best_converged, na.rm=TRUE),4),
        "\n"
      )
    }
    
  }
  
  
  cat("\nResults saved in:\n",out_dir,"\n")
  cat("Simulation finished.\n")