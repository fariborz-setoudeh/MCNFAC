
# install.packages("mvtnorm")
# install.packages("MomTrunc")  




# ============================================================
# Robust Mixture of Factor-Analytic Regression under
# Detection-limit censoring + contaminated-normal (two-point W)
# AECM implementation 
#-----------------------------------------------------------------
# PART 1: foundations + data interface + TMVN moments + helpers
# ============================================================

rmfac_require_packages <- function() {
  pkgs <- c("mvtnorm", "MomTrunc")
  miss <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(miss) > 0) {
    stop("Missing packages: ", paste(miss, collapse = ", "),
         ". Install them before running.")
  }
  suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE))
  invisible(TRUE)
}

rmfac_require_packages()

# --------------------------
# Numerics / utilities
# --------------------------
rmfac_logsumexp <- function(a) {
  m <- max(a)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(a - m)))
}

rmfac_sym <- function(M) 0.5 * (M + t(M))

rmfac_chol_safe <- function(S, jitter = 1e-8, max_tries = 8) {
  # Cholesky with adaptive jitter to handle near-PD matrices
  S <- rmfac_sym(S)
  for (k in 0:max_tries) {
    eps <- jitter * (10^k)
    out <- try(chol(S + diag(eps, nrow(S))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  stop("Cholesky failed even after jitter. Matrix likely not PD.")
}

rmfac_logdmvnorm <- function(x, mean, sigma) {
  # stable log-density
  mvtnorm::dmvnorm(x = x, mean = mean, sigma = sigma, log = TRUE)
}

# ------------------------------------------------------------
# ------------------------------------------------------------
# TMVN MOMENTS via MomTrunc 
# ------------------------------------------------------------
rmfac_tmvn_moments <- function(lower, upper, mu, Sigma) {
  out <- MomTrunc::meanvarTMD(
    mu = as.numeric(mu),
    Sigma = rmfac_sym(Sigma),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    dist = "normal"
  )
  list(Ey = as.numeric(out$mean), Eyy = out$EYY)
}
  

# --------------------------
# DATA INTERFACE (strict)
# --------------------------
# You must provide:
#   y_obs : numeric matrix n x p, with NA for censored components (optional;
#           can be ignored if you provide l,u fully)
#   l_mat : numeric matrix n x p, lower bounds (can equal upper for exact)
#   u_mat : numeric matrix n x p, upper bounds (can equal lower for exact)
#   X     : numeric matrix n x d (covariates); use matrix(1,n,1) for intercept-only
#
# Exact observation convention:
#   component k of unit j is exact iff l_mat[j,k] == u_mat[j,k]
#
# Fully censored means all l!=u; fully observed means all l==u.

rmfac_build_censor_info <- function(l_row, u_row) {
  # For a single observation j (vectors length p)
  p <- length(l_row)
  if (length(u_row) != p) stop("l_row/u_row length mismatch")
  
  is_exact <- (l_row == u_row)
  idx_o <- which(is_exact)
  idx_c <- which(!is_exact)
  
  # selection matrices are huge; we store indices and use subsetting
  list(
    p = p,
    idx_o = idx_o,
    idx_c = idx_c,
    po = length(idx_o),
    pc = length(idx_c),
    l_o = if (length(idx_o) > 0) l_row[idx_o] else numeric(0),
    u_o = if (length(idx_o) > 0) u_row[idx_o] else numeric(0),
    l_c = if (length(idx_c) > 0) l_row[idx_c] else numeric(0),
    u_c = if (length(idx_c) > 0) u_row[idx_c] else numeric(0)
  )
}

rmfac_build_all_censor_info <- function(l_mat, u_mat) {
  if (!is.matrix(l_mat) || !is.matrix(u_mat)) stop("l_mat/u_mat must be matrices")
  if (!all(dim(l_mat) == dim(u_mat))) stop("l_mat/u_mat dimension mismatch")
  n <- nrow(l_mat)
  lapply(seq_len(n), function(j) rmfac_build_censor_info(l_mat[j, ], u_mat[j, ]))
}

# --------------------------
# Model parameter container
# --------------------------
# params is a list with:
#   g        : number of components
#   p,q,d    : dimensions
#   pi       : length g, sums to 1
#   beta     : list length g, each (p*d) vector (vec(Beta_i))
#   B        : list length g, each p x q
#   D        : list length g, each p x p diagonal matrix (or store diag vector)
#   nu       : length g, contamination prob (P(V=0))
#   eta      : length g, inflation > 1

rmfac_check_params <- function(params) {
  req <- c("g","p","q","d","pi","beta","B","D","nu","eta")
  miss <- setdiff(req, names(params))
  if (length(miss) > 0) stop("params missing: ", paste(miss, collapse = ", "))
  
  g <- params$g; p <- params$p; q <- params$q; d <- params$d
  if (length(params$pi) != g) stop("pi length != g")
  if (abs(sum(params$pi) - 1) > 1e-6) stop("pi must sum to 1")
  if (any(params$pi <= 0)) stop("pi must be > 0")
  
  if (length(params$beta) != g) stop("beta list length != g")
  if (length(params$B)    != g) stop("B list length != g")
  if (length(params$D)    != g) stop("D list length != g")
  if (length(params$nu)   != g) stop("nu length != g")
  if (length(params$eta)  != g) stop("eta length != g")
  
  for (i in 1:g) {
    if (length(params$beta[[i]]) != p*d) stop("beta[[",i,"]] length != p*d")
    Bi <- params$B[[i]]
    if (!all(dim(Bi) == c(p,q))) stop("B[[",i,"]] must be p x q")
    Di <- params$D[[i]]
    if (!all(dim(Di) == c(p,p))) stop("D[[",i,"]] must be p x p")
    if (any(diag(Di) <= 0)) stop("D[[",i,"]] must have positive diagonal")
    if (!isTRUE(all.equal(Di, diag(diag(Di))))) stop("D[[",i,"]] must be diagonal")
    if (!(params$nu[i] > 0 && params$nu[i] < 1)) stop("nu must be in (0,1)")
    if (!(params$eta[i] > 1)) stop("eta must be > 1")
  }
  
  invisible(TRUE)
}

# --------------------------
# Build X_j = I_p ⊗ x_j' and compute mu_ij
# --------------------------
rmfac_Xj_big <- function(xj, p) {
  # Correct X_j = I_p ⊗ x_j^T  →  dimension p × (p d)
  xj <- as.numeric(xj)
  kronecker(diag(p), t(xj))
}


rmfac_mu_ij <- function(xj, beta_i, p) {
  # mu_ij = Xj * beta_i  => p-vector
  Xj <- rmfac_Xj_big(xj, p)
  as.numeric(Xj %*% beta_i)
}

rmfac_Sigma_i <- function(Bi, Di) {
  rmfac_sym(Bi %*% t(Bi) + Di)
}

# --------------------------
# Two-point W from V (
# V=1 typical, V=0 contaminated, W in {1, eta^{-1}}
# --------------------------
rmfac_W_from_V <- function(V, eta) {
  V + (1 - V)/eta
}

# Done: Part 1
cat("RMFAC PART 1 loaded: packages, TMVN moments, censor-info builder, param checks.\n")
# ============================================================
# PART 2: L_ij(t) for 3 censoring scenarios + TMVN moments
#   (i) no censoring
#   (ii) fully censored
#   (iii) mixed censoring
#
# Returns:
#   logL_ij(t)
#   yhat(t)=E(y*|obs,z=1,W=t),  Yhat(t)=E(y*y'|...)
# ============================================================

# ---- Stable MVN log-density (avoid -Inf/NaN) ----
rmfac_logdmvnorm_safe <- function(x, mean, sigma) {
  val <- mvtnorm::dmvnorm(x = x, mean = mean, sigma = sigma, log = TRUE)
  if (!is.finite(val)) val <- -1e300
  val
}

# ---- Conditional Gaussian params for mixed censoring ----
rmfac_cond_gaussian_params <- function(mu_o, mu_c, Sigma_oo, Sigma_oc, Sigma_co, Sigma_cc, y_o) {
  # Conditional of y_c | y_o for joint Gaussian with mean (mu_o, mu_c).
  # Scale factor t^{-1} cancels in conditional mean; conditional covariance scales by t^{-1}.
  if (length(mu_o) == 0) {
    return(list(mu_c_given_o = mu_c, Sigma_cc_given_o = Sigma_cc))
  }
  chol_oo <- rmfac_chol_safe(Sigma_oo)
  diff_o  <- y_o - mu_o
  tmp     <- backsolve(chol_oo, forwardsolve(t(chol_oo), diff_o))  # Sigma_oo^{-1} diff_o
  
  mu_c_go <- as.numeric(mu_c + Sigma_co %*% tmp)
  
  tmp2 <- backsolve(chol_oo, forwardsolve(t(chol_oo), Sigma_oc))   # Sigma_oo^{-1} Sigma_oc
  Sig_cc_go <- rmfac_sym(Sigma_cc - Sigma_co %*% tmp2)
  
  list(mu_c_given_o = mu_c_go, Sigma_cc_given_o = Sig_cc_go)
}

# ---- Stable rectangle probability log P(lower<=Y<=upper) ----
rmfac_prob_rect <- function(lower, upper, mu, Sigma) {
  d <- length(mu)
  if (d == 0) return(list(logprob = 0))
  
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  mu <- as.numeric(mu)
  Sigma <- rmfac_sym(Sigma)
  
  if (d == 1) {
    sd <- sqrt(as.numeric(Sigma))
    pr <- pnorm(upper, mean = mu, sd = sd) -
      pnorm(lower, mean = mu, sd = sd)
    return(list(logprob = log(max(pr, 1e-300))))
  }
  
  pr <- try(
    mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean = mu,
      sigma = Sigma,
      algorithm = mvtnorm::Miwa(steps = 128)
    ),
    silent = TRUE
  )
  
  if (inherits(pr, "try-error")) return(list(logprob = -1e300))
  
  list(logprob = log(max(as.numeric(pr), 1e-300)))
}

# ---- Core: compute L_ij(t) and state-specific TMVN moments ----
rmfac_Lij_t <- function(j, i, t, l_mat, u_mat, X, params, censor_info_list) {
  # Model: y* | z=1, W=t ~ N_p(mu_ij, t^{-1} Sigma_i),   t in {1, eta_i^{-1}}
  # Returns: list(logL, yhat, Yhat)
  
  p <- params$p
  
  # mean and covariance
  xj    <- X[j, ]
  mu    <- rmfac_mu_ij(xj, params$beta[[i]], p)          # p-vector
  Sigma <- rmfac_Sigma_i(params$B[[i]], params$D[[i]])   # p x p
  Sigma_t <- (1 / t) * Sigma
  
  info <- censor_info_list[[j]]
  idx_o <- info$idx_o
  idx_c <- info$idx_c
  po <- info$po
  pc <- info$pc
  
  # Observed exact values are those with l==u
  y_o <- if (po > 0) as.numeric(l_mat[j, idx_o]) else numeric(0)
  
  # (i) no censoring (pc=0): all exact
  if (pc == 0) {
    logL <- rmfac_logdmvnorm_safe(
      x = y_o,
      mean = mu[idx_o],
      sigma = Sigma_t[idx_o, idx_o, drop = FALSE]
    )
    yhat <- as.numeric(l_mat[j, ])  # equals u_mat row
    Yhat <- tcrossprod(yhat)
    if (!is.finite(logL)) logL <- -1e300
    return(list(logL = logL, yhat = yhat, Yhat = Yhat))
  }
  
  # (ii) fully censored (po=0)
  if (po == 0) {
    lower <- as.numeric(l_mat[j, ])
    upper <- as.numeric(u_mat[j, ])
    
    logL <- rmfac_prob_rect(lower, upper, mu, Sigma_t)$logprob
    
    mom <- try(rmfac_tmvn_moments(lower = lower, upper = upper, mu = mu, Sigma = Sigma_t),
               silent = TRUE)
    if (inherits(mom, "try-error")) {
      cat("TMVN MOMENT FAILURE: fully censored case; j=", j, " i=", i, " t=", t, "\n")
      
      if (!is.finite(logL)) logL <- -1e300
      
      return(list(
        logL = logL,
        yhat = mu,
        Yhat = Sigma_t + tcrossprod(mu)
      ))
    }
    
    if (!is.finite(logL)) logL <- -1e300
    return(list(logL = logL, yhat = mom$Ey, Yhat = mom$Eyy))
  }
  
  # (iii) mixed censoring (po>0, pc>0)
  mu_o <- mu[idx_o]
  mu_c <- mu[idx_c]
  
  Sigma_oo <- Sigma[idx_o, idx_o, drop = FALSE]
  Sigma_oc <- Sigma[idx_o, idx_c, drop = FALSE]
  Sigma_co <- Sigma[idx_c, idx_o, drop = FALSE]
  Sigma_cc <- Sigma[idx_c, idx_c, drop = FALSE]
  
  # Density of observed exact part under scaled covariance
  log_dens_o <- rmfac_logdmvnorm_safe(
    x = y_o,
    mean = mu_o,
    sigma = (1 / t) * Sigma_oo
  )
  
  # Conditional params (unscaled); then scale conditional covariance by t^{-1}
  cond <- rmfac_cond_gaussian_params(mu_o, mu_c, Sigma_oo, Sigma_oc, Sigma_co, Sigma_cc, y_o)
  mu_c_go <- cond$mu_c_given_o
  Sig_cc_go_t <- (1 / t) * cond$Sigma_cc_given_o
  
  lower_c <- as.numeric(l_mat[j, idx_c])
  upper_c <- as.numeric(u_mat[j, idx_c])
  
  log_prob_c <- rmfac_prob_rect(lower_c, upper_c, mu_c_go, Sig_cc_go_t)$logprob
  logL <- log_dens_o + log_prob_c
  
  # Moments for censored block: y_c | y_o, W=t, y_c in [l_c,u_c]
  mom_c <- try(rmfac_tmvn_moments(lower = lower_c, upper = upper_c, mu = mu_c_go, Sigma = Sig_cc_go_t),
               silent = TRUE)
  if (inherits(mom_c, "try-error")) {
    cat("TMVN MOMENT FAILURE: mixed case; j=", j, " i=", i, " t=", t, "\n")
    
    if (!is.finite(logL)) logL <- -1e300
    
    return(list(
      logL = logL,
      yhat = mu,
      Yhat = Sigma_t + tcrossprod(mu)
    ))
  }
  
  Ey_c  <- mom_c$Ey
  Eyy_c <- mom_c$Eyy
  
  # Reconstruct full-vector moments
  yhat <- numeric(p)
  yhat[idx_o] <- y_o
  yhat[idx_c] <- Ey_c
  
  Yhat <- matrix(0, p, p)
  Yhat[idx_o, idx_o] <- tcrossprod(y_o)
  Yhat[idx_o, idx_c] <- y_o %*% t(Ey_c)
  Yhat[idx_c, idx_o] <- Ey_c %*% t(y_o)
  Yhat[idx_c, idx_c] <- Eyy_c
  
  if (!is.finite(logL)) logL <- -1e300
  list(logL = logL, yhat = yhat, Yhat = Yhat)
}

# ---- Convenience: compute both states t=1 and t=eta^{-1} ----
rmfac_Lij_states <- function(j, i, l_mat, u_mat, X, params, censor_info_list) {
  t1 <- 1
  t0 <- 1 / params$eta[i]  # eta^{-1}
  
  out1 <- rmfac_Lij_t(j, i, t1, l_mat, u_mat, X, params, censor_info_list)
  out0 <- rmfac_Lij_t(j, i, t0, l_mat, u_mat, X, params, censor_info_list)
  
  list(
    t = c(t1, t0),
    logL = c(out1$logL, out0$logL),
    yhat = list(`1` = out1$yhat, `eta_inv` = out0$yhat),
    Yhat = list(`1` = out1$Yhat, `eta_inv` = out0$Yhat)
  )
}

cat("RMFAC PART 2 loaded: L_ij(t) + TMVN moments for 3 censoring scenarios.\n")
# ============================================================

# PART 3 — FULL E-STEP 
# ============================================================

rmfac_E_step <- function(l_mat, u_mat, X, params) {
  
  rmfac_check_params(params)
  
  n <- nrow(l_mat)
  g <- params$g
  p <- params$p
  
  censor_info_list <- rmfac_build_all_censor_info(l_mat, u_mat)
  
  # storage
  logL1  <- matrix(0, n, g)   # log L_ij(1)
  logL0  <- matrix(0, n, g)   # log L_ij(eta^{-1})
  
  yhat1  <- array(0, c(n,p,g))  # E[y* | W=1]
  yhat0  <- array(0, c(n,p,g))  # E[y* | W=eta^{-1}]
  Yhat1  <- array(0, c(p,p,n,g))
  Yhat0  <- array(0, c(p,p,n,g))
  
  # --------------------------------------------------
  # Step 1: L_ij(t) + truncated moments
  # --------------------------------------------------
  for (i in 1:g) {
    for (j in 1:n) {
      
      out <- rmfac_Lij_states(j, i, l_mat, u_mat, X, params, censor_info_list)
      
      logL1[j,i] <- out$logL[1]
      logL0[j,i] <- out$logL[2]
      
      yhat1[j,,i] <- out$yhat[[1]]
      yhat0[j,,i] <- out$yhat[[2]]
      
      Yhat1[,,j,i] <- out$Yhat[[1]]
      Yhat0[,,j,i] <- out$Yhat[[2]]
    }
  }
  
  # --------------------------------------------------
  # Step 2: zhat
  # --------------------------------------------------
  zhat <- matrix(0,n,g)
  
  for (j in 1:n) {
    log_num <- rep(0,g)
    for (i in 1:g) {
      a <- log(params$pi[i]) + log(1-params$nu[i]) + logL1[j,i]
      b <- log(params$pi[i]) + log(params$nu[i])   + logL0[j,i]
      log_num[i] <- rmfac_logsumexp(c(a,b))
    }
    log_den <- rmfac_logsumexp(log_num)
    zhat[j,] <- exp(log_num - log_den)
  }
  
  # --------------------------------------------------
  # Step 3: vhat (typical-state posterior)
  # --------------------------------------------------
  vhat <- matrix(0,n,g)
  
  for (i in 1:g) {
    for (j in 1:n) {
      a <- log(1-params$nu[i]) + logL1[j,i]
      b <- log(params$nu[i])   + logL0[j,i]
      denom <- rmfac_logsumexp(c(a,b))
      vhat[j,i] <- exp(a - denom)
    }
  }
  
  # --------------------------------------------------
  # Step 4: W moments
  # --------------------------------------------------
  omegahat <- vhat
  ellW     <- vhat * 0
  
  for(i in 1:g){
    omegahat[,i] <- vhat[,i] + (1 - vhat[,i]) / params$eta[i]
    ellW[,i]     <- -(1 - vhat[,i]) * log(params$eta[i])
  }
  
  # --------------------------------------------------
  # Step 5: W-weighted conditional response moments
  # Wyhat = E(W y* | obs, z_ij = 1)
  # WYhat = E(W y* y*' | obs, z_ij = 1)
  # --------------------------------------------------
  Wyhat <- array(0, c(n,p,g))
  WYhat <- array(0, c(p,p,n,g))
  
  for (i in 1:g) {
    t0 <- 1/params$eta[i]
    for (j in 1:n) {
      
      tau1 <- vhat[j,i]
      tau0 <- 1 - tau1
      
      Wyhat[j,,i]  <- tau1*yhat1[j,,i] + tau0*t0*yhat0[j,,i]
      WYhat[,,j,i] <- tau1*Yhat1[,,j,i] + tau0*t0*Yhat0[,,j,i]
    }
  }
  
  list(
    zhat=zhat,
    vhat=vhat,
    omegahat=omegahat,
    ellW=ellW,
    yhat1=yhat1,
    yhat0=yhat0,
    Yhat1=Yhat1,
    Yhat0=Yhat0,
    Wyhat=Wyhat,
    WYhat=WYhat,
    logL1=logL1,
    logL0=logL0
  )
}

cat("RMFAC PART 3  loaded.\n")
# ============================================================
# PART 4 — Cycle 1 CM updates: pi, nu, eta
#   Uses:
#     zhat, vhat from E-step
#     contaminated-state moments (yhat0, Yhat0) for delta^(0)
# ============================================================

rmfac_CM1_update <- function(params, estep, X, eta_eps = 1e-6) {
  
  g <- params$g
  p <- params$p
  n <- nrow(estep$zhat)
  
  zhat <- estep$zhat
  vhat <- estep$vhat
  
  # --------------------------
  # Update pi
  # --------------------------
  pi_new <- colMeans(zhat)
  pi_new <- pmax(pi_new, 1e-12)
  pi_new <- pi_new / sum(pi_new)
  
  # --------------------------
  # Update nu (contamination prob = P(V=0))
  # nu_i = sum_j zhat_ij * (1 - vhat_ij) / sum_j zhat_ij
  # --------------------------
  n_i <- colSums(zhat)
  a_i <- colSums(zhat * (1 - vhat))  # expected # contaminated in comp i
  
  nu_new <- a_i / pmax(n_i, 1e-300)
  nu_new <- pmin(pmax(nu_new, 1e-8), 1 - 1e-8)
  
  # --------------------------
  # Update eta
  # eta_i = max(1+eps, b_i / (p * a_i))
  # where b_i = sum_j zhat_ij * (1-vhat_ij) * delta_ij^(0)
  # and delta_ij^(0) = tr( Sigma_i^{-1} R_ij^(0) )
  # with R_ij^(0) built from state t = eta^{-1} moments (yhat0, Yhat0)
  # --------------------------
  eta_new <- params$eta
  
  for (i in 1:g) {
    
    if (a_i[i] <= 0) {
      # no contaminated mass -> eta not identified this iter
      eta_new[i] <- params$eta[i]
      next
    }
    
    Sigma_i <- rmfac_Sigma_i(params$B[[i]], params$D[[i]])
    Sigma_inv <- solve(Sigma_i)
    
    b_sum <- 0
    
    for (j in 1:n) {
      mu_ij <- rmfac_mu_ij(X[j,], params$beta[[i]], p)
      
      y0 <- estep$yhat0[j,,i]
      Y0 <- estep$Yhat0[,,j,i]
      
      # R_ij^(0) = E( (y*-mu)(y*-mu)' | W=eta^{-1}, obs )
      R0 <- Y0 - y0 %*% t(mu_ij) - mu_ij %*% t(y0) + tcrossprod(mu_ij)
      
      delta0_ij <- sum(diag(Sigma_inv %*% R0))
      
      b_sum <- b_sum + zhat[j,i] * (1 - vhat[j,i]) * delta0_ij
    }
    
    eta_star <- b_sum / (p * a_i[i])
    
    # ---- HARD CONSTRAINT: eta must be strictly > 1 ----
    eta_tmp <- max(1 + eta_eps, eta_star)
    
    
    # protect against NaN / Inf / collapse to 1
    if(!is.finite(eta_tmp) || eta_tmp <= 1){
      eta_tmp <- 1 + eta_eps
    }
    
    eta_new[i] <- eta_tmp
    
  }
  
  params_new <- params
  params_new$pi  <- pi_new
  params_new$nu  <- nu_new
  params_new$eta <- eta_new
  
  rmfac_check_params(params_new)
  params_new
}

cat("RMFAC PART 4 loaded: Cycle-1 CM updates (pi, nu, eta).\n")
# ============================================================
# PART 5 — INITIALIZATION
# Creates the initial "params" object required by the EM/AECM
# ============================================================


# ------------------------------------------------------------
# 5.1 initialization imputation via TMVN expectation
# 
# ------------------------------------------------------------
rmfac_init_impute_TMNV <- function(l_mat, u_mat) {
  
  n <- nrow(l_mat)
  p <- ncol(l_mat)
  
  Yimp <- matrix(0, n, p)
  
  mu0 <- rep(0, p)
  Sigma0 <- diag(p)
  
  for(j in 1:n){
    
    lower <- l_mat[j,]
    upper <- u_mat[j,]
    
    # fully observed → keep value
    if(all(lower == upper)){
      Yimp[j,] <- lower
      next
    }
    
    # censored → use TMVN mean under N(0,I)
    mom <- rmfac_tmvn_moments(lower, upper, mu0, Sigma0)
    Yimp[j,] <- mom$Ey
  }
  
  Yimp
}

# ------------------------------------------------------------
# 5.1b Finite initialization values for censored observations
#
# This function is used only to initialize k-means, beta, B and D.
# It does not alter the censoring bounds used in the likelihood.
# ------------------------------------------------------------
rmfac_init_impute_bounds <- function(l_mat, u_mat) {
  
  if(!is.matrix(l_mat) || !is.matrix(u_mat)){
    stop("l_mat and u_mat must be matrices.")
  }
  
  if(!all(dim(l_mat) == dim(u_mat))){
    stop("l_mat and u_mat must have identical dimensions.")
  }
  
  Yimp <- matrix(
    NA_real_,
    nrow = nrow(l_mat),
    ncol = ncol(l_mat)
  )
  
  # Exact observations
  exact <- is.finite(l_mat) &
    is.finite(u_mat) &
    (l_mat == u_mat)
  
  # Left-censored observations: (-Inf, upper]
  left_censored <- is.infinite(l_mat) &
    (l_mat < 0) &
    is.finite(u_mat)
  
  # Right-censored observations: [lower, Inf)
  right_censored <- is.finite(l_mat) &
    is.infinite(u_mat) &
    (u_mat > 0)
  
  # Finite interval-censored observations
  interval_censored <- is.finite(l_mat) &
    is.finite(u_mat) &
    (l_mat < u_mat)
  
  Yimp[exact] <- l_mat[exact]
  
  # For initialization only, use the censoring threshold.
  Yimp[left_censored] <- u_mat[left_censored]
  Yimp[right_censored] <- l_mat[right_censored]
  
  # Use interval midpoint for finite interval censoring.
  Yimp[interval_censored] <-
    0.5 * (
      l_mat[interval_censored] +
        u_mat[interval_censored]
    )
  
  if(any(!is.finite(Yimp))){
    stop(
      "Finite initialization failed: Yimp contains NA, NaN or Inf."
    )
  }
  
  Yimp
}
# ------------------------------------------------------------
# 5.2 initialize beta by least squares inside each cluster
# ------------------------------------------------------------
rmfac_init_beta_cluster <- function(Y, X) {
  
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)
  
  # Build stacked design matrix X*
  Xbig <- do.call(rbind,
                  lapply(1:n, function(j) rmfac_Xj_big(X[j,], p)))
  # Xbig : (n p) × (p d)
  
  # Correct vec(Y) stacking (column-wise stacking by variable)
  # vec(Y) = (y1^T, y2^T, ..., yp^T)^T
  yvec <- as.vector(t(Y))
  
  
  beta <- solve(t(Xbig) %*% Xbig + diag(1e-6, ncol(Xbig))) %*%
    t(Xbig) %*% yvec
  
  as.numeric(beta)
}


# ------------------------------------------------------------
# 5.3 initialize factor-analytic covariance via PCA
# ------------------------------------------------------------
rmfac_init_FA <- function(Y, q) {
  
  p <- ncol(Y)
  n <- nrow(Y)
  
  # ---- robust covariance (works even for small clusters)
  if(n <= p){
    S <- diag(apply(Y,2,var) + 1e-3)   # fallback diagonal covariance
  } else {
    S <- cov(Y)
  }
  
  # numerical safety
  S <- S + diag(1e-6, p)
  
  # eigen decomposition
  eig <- eigen(S, symmetric = TRUE)
  
  # ensure non-negative eigenvalues
  eigvals <- pmax(eig$values, 1e-8)
  
  # IMPORTANT: q cannot exceed rank
  q_eff <- min(q, sum(eigvals > 1e-6))
  
  if(q_eff == 0){
    # pure diagonal model fallback
    B <- matrix(0, p, q)
    D <- diag(diag(S))
    return(list(B=B, D=D))
  }
  
  # factor loadings
  Btemp <- eig$vectors[,1:q_eff,drop=FALSE] %*% 
    diag(sqrt(eigvals[1:q_eff]), q_eff)
  
  # pad with zeros if q > q_eff
  if(q_eff < q){
    B <- matrix(0,p,q)
    B[,1:q_eff] <- Btemp
  } else {
    B <- Btemp
  }
  
  # uniquenesses
  D <- diag(pmax(diag(S - B %*% t(B)), 1e-6))
  
  list(B=B, D=D)
}


# ------------------------------------------------------------
# 5.4 MAIN INITIALIZER — builds params list
# ------------------------------------------------------------
rmfac_initialize <- function(l_mat, u_mat, X, g, q) {
  
  Ymid <- rmfac_init_impute_bounds(
    l_mat = l_mat,
    u_mat = u_mat
  )
  
  if(any(!is.finite(Ymid))){
    stop(
      "Initialization matrix contains NA, NaN or Inf before k-means."
    )
  }
  
  n <- nrow(Ymid)
  p <- ncol(Ymid)
  d <- ncol(X)
  
  # K-means clustering
  cl <- kmeans(Ymid, centers = g, nstart = 20)$cluster
  
  pi <- table(cl) / n
  pi <- as.numeric(pi)
  pi <- pmax(pi,1e-6)
  pi <- pi / sum(pi)
  
  beta_list <- vector("list", g)
  B_list    <- vector("list", g)
  D_list    <- vector("list", g)
  
  for (i in 1:g) {
    
    idx <- which(cl == i)
    if (length(idx) < 2) idx <- sample(1:n, 2)
    
    Yi <- Ymid[idx,,drop=FALSE]
    Xi <- X[idx,,drop=FALSE]
    
    beta_list[[i]] <- rmfac_init_beta_cluster(Yi, Xi)
    
    FA <- rmfac_init_FA(Yi, q)
    B_list[[i]] <- FA$B
    D_list[[i]] <- FA$D
  }
  
  params <- list(
    g=g, p=p, q=q, d=d,
    pi=pi,
    beta=beta_list,
    B=B_list,
    D=D_list,
    nu=rep(0.05,g),
    eta=rep(10,g)
  )
  
  rmfac_check_params(params)
  params
}

cat("RMFAC PART 5 loaded: initialization ready.\n")
#---------------------------------------------------------
#---------------------------------------------------------

# ============================================================
# Cycle-2 helper
#
# All conditional moments below are computed from the E-step performed at
# iteration h.  The factor moments therefore remain based on mu^(h), B^(h),
# and D^(h), even after beta has been updated in Cycle 1.  Only the response
# residual moments used in the B and D updates are re-centered at mu^(h+1).
# ============================================================

rmfac_cycle2_stats <- function(params_old, params_new, estep, Xrow, i, j)
{
  # old parameters (h)
  B_old <- params_old$B[[i]]
  D_old <- params_old$D[[i]]
  Sigma_old <- rmfac_Sigma_i(B_old, D_old)
  Sigma_old_inv <- solve(Sigma_old)
  
  Gamma <- t(B_old) %*% Sigma_old_inv
  Omega <- diag(params_old$q) - Gamma %*% B_old
  
  # regression means at iterations h and h+1
  mu_old <- rmfac_mu_ij(Xrow, params_old$beta[[i]], params_old$p)
  mu_new <- rmfac_mu_ij(Xrow, params_new$beta[[i]], params_new$p)
  
  omega_hat <- estep$omegahat[j,i]
  ty <- estep$Wyhat[j,,i]
  tY <- estep$WYhat[,,j,i]

  # WR^(h): E_h[W (y* - mu^h)(y* - mu^h)' | obs, z_ij = 1]
  WR_old <- tY -
    ty %*% t(mu_old) -
    mu_old %*% t(ty) +
    omega_hat * tcrossprod(mu_old)

  # Wf^(h): E_h[W f | obs, z_ij = 1]
  Wf_old <- as.numeric(
    Gamma %*% (ty - omega_hat * mu_old)
  )

  # WF^(h): E_h[W f f' | obs, z_ij = 1]
  WF_old <- Omega + Gamma %*% WR_old %*% t(Gamma)

  # WR^(h,+): the same E-step response moments re-centered at mu^(h+1)
  WR_plus <- tY -
    ty %*% t(mu_new) -
    mu_new %*% t(ty) +
    omega_hat * tcrossprod(mu_new)

  # WS^(h,+): E_h[W (y* - mu^(h+1)) f' | obs, z_ij = 1]
  WS_plus <- WR_old %*% t(Gamma) +
    (mu_old - mu_new) %*% t(Wf_old)

  list(
    WR_old  = WR_old,
    WR_plus = WR_plus,
    Wf_old  = Wf_old,
    WF_old  = WF_old,
    WS_plus = WS_plus
  )
}
# ============================================================
#new
# ============================================================
# ============================================================
# CM update for beta
# Marginal AECM regression update:
# factors are integrated out and Sigma_i = B_i B_i' + D_i
# ============================================================

rmfac_CM1_update_beta <- function(params, estep, X)
{
  g <- params$g
  n <- nrow(estep$zhat)
  p <- params$p
  d <- params$d
  
  beta_new <- vector("list", g)
  
  for(i in seq_len(g))
  {
    Sigma_i <- rmfac_Sigma_i(
      params$B[[i]],
      params$D[[i]]
    )
    
    Sigma_inv <- solve(Sigma_i)
    
    A <- matrix(0, p * d, p * d)
    b <- numeric(p * d)
    
    for(j in seq_len(n))
    {
      z_ij     <- estep$zhat[j, i]
      omega_ij <- estep$omegahat[j, i]
      
      if(z_ij < 1e-12) next
      
      Xj <- rmfac_Xj_big(X[j, ], p)
      
      # A = sum z_ij omega_ij X_j' Sigma_i^{-1} X_j
      A <- A +
        z_ij * omega_ij *
        t(Xj) %*% Sigma_inv %*% Xj
      
      # Wyhat is already E(W_ij y*_j | obs, z_ij = 1).
      # Therefore omega_ij must NOT be multiplied again here.
      b <- b +
        z_ij *
        t(Xj) %*% Sigma_inv %*%
        estep$Wyhat[j, , i]
    }
    
    beta_new[[i]] <- as.numeric(
      solve(
        A + diag(1e-8, p * d),
        b
      )
    )
  }
  
  params$beta <- beta_new
  params
}

# ============================================================
# CM update for B 
# ============================================================

rmfac_CM2_update_B <- function(params_old, params_new, estep, X)
{
  g <- params_old$g
  n <- nrow(estep$zhat)
  p <- params_old$p
  q <- params_old$q
  
  B_new <- vector("list", g)
  
  for(i in 1:g)
  {
    S_sum <- matrix(0,p,q)
    F_sum <- matrix(0,q,q)
    
    for(j in 1:n)
    {
      z_ij <- estep$zhat[j,i]
      if(z_ij < 1e-12) next
      
      stats <- rmfac_cycle2_stats(params_old, params_new, estep, X[j,], i, j)
      
      # WS_plus and WF_old already contain the factor W_ij.
      S_sum <- S_sum + z_ij * stats$WS_plus
      F_sum <- F_sum + z_ij * stats$WF_old
    }
    
    B_new[[i]] <- S_sum %*% solve(F_sum + diag(1e-8,q))
  }
  
  params_new$B <- B_new
  params_new
}
# ============================================================
# CM update for D 
# ============================================================

rmfac_CM2_update_D <- function(params_old, params_new, estep, X, d_floor = 1e-6)
{
  g <- params_old$g
  n <- nrow(estep$zhat)
  p <- params_old$p
  
  D_new <- vector("list", g)
  
  for(i in 1:g)
  {
    n_i <- sum(estep$zhat[,i])
    if(n_i < 1e-12){
      D_new[[i]] <- params_old$D[[i]]
      next
    }
    
    A_sum <- matrix(0,p,p)
    
    for(j in 1:n)
    {
      z_ij <- estep$zhat[j,i]
      if(z_ij < 1e-12) next
      
      stats <- rmfac_cycle2_stats(params_old, params_new, estep, X[j,], i, j)
      
      B_new_i <- params_new$B[[i]]
      
      # E_ij(B^{h+1}, μ^{h+1})  (Eq E_mu_new)
      Eij <- stats$WR_plus -
        B_new_i %*% t(stats$WS_plus) -
        stats$WS_plus %*% t(B_new_i) +
        B_new_i %*% stats$WF_old %*% t(B_new_i)
      
      # Every matrix in Eij is already W_ij-weighted.
      A_sum <- A_sum + z_ij * Eij
    }
    
    Dvec <- diag(A_sum / n_i)
    Dvec <- pmax(Dvec, d_floor)
    D_new[[i]] <- diag(Dvec, p)
  }
  
  params_new$D <- D_new
  params_new
}




# ============================================================
# PART 10 — Observed log-likelihood of the model
# ============================================================

rmfac_loglikelihood <- function(l_mat, u_mat, X, params){
  
  rmfac_check_params(params)
  
  n <- nrow(l_mat)
  g <- params$g
  
  censor_info_list <- rmfac_build_all_censor_info(l_mat, u_mat)
  
  loglik <- 0
  
  for(j in 1:n){
    
    log_terms <- rep(0, g)
    
    for(i in 1:g){
      
      out <- rmfac_Lij_states(j, i, l_mat, u_mat, X, params, censor_info_list)
      
      logL1 <- out$logL[1]   # W=1
      logL0 <- out$logL[2]   # W=eta^{-1}
      
      # log( (1-ν)L1 + νL0 )
      a <- log(1 - params$nu[i]) + logL1
      b <- log(params$nu[i])     + logL0
      
      log_mixW <- rmfac_logsumexp(c(a, b))
      
      # log π_i + log mixture over W
      log_terms[i] <- log(params$pi[i]) + log_mixW
    }
    
    # log sum over components
    loglik <- loglik + rmfac_logsumexp(log_terms)
  }
  
  as.numeric(loglik)
}

cat("RMFAC PART 10 loaded: log-likelihood function.\n")








#####################################



################
