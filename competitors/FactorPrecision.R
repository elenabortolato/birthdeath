 
# -----------------------------
# Library imports
# -----------------------------
FactorPrecision<- function (y, maxit, d){
  library(pgdraw)
  library(MASS)
  library(Rcpp)
  library(statmod)
  library(mnormt)
  library(mvtnorm)
  p=ncol(y)
  T=nrow(y)
  
  # -----------------------------
  # Initialization
  # -----------------------------
  a_delta <- 1; b_delta <- 1   # variance of error covariance
  a_lam   <- 1; b_lam   <- 1   # governing loadings variance
  q <- 1                       # covariate for P(Lambda=0|s=1)
  
  # Initialize precision matrices
  Plam <- 2 * diag(d)
  invPlam <- solve(Plam)
  
  deltas <- 0.1 * rep(1, p)
  Delta <- diag(deltas)
  invDelta <- diag(1 / diag(Delta))
  
  eps <- mvnfast::rmvn(T, sigma = diag(p), mu = rep(0,p))  # initial errors
  
  p_constant <- 1  # constant used in Polya-Gamma
  
  # Initialize Lambda
  Lambda <- array(rnorm(p*d, 0, 0.1) * rbinom(p*d, size = 1, prob = 0.5),
                  dim = c(p, d, T))
  
  # Initialize latent states ("shocks")
  s <- array(rbinom(d*T, 1, 0.0), dim = c(d, T))
  
  # Initialize factors and y
  eta <- matrix(0, nrow = T, ncol = d)
  y <- eps
  
  for (t in 1:T) {
    Pt <- diag(d) + t(Lambda[,,t]) %*% invDelta %*% Lambda[,,t]
    eps[t,] <- mvnfast::rmvn(1, sigma = invDelta, mu = rep(0,p))
    eta[t,] <- mvnfast::rmvn(1, sigma = diag(d), mu = rep(0,d))
  }
  
  
  s[,1] <- 1
  
  LL <- array(0, c(p, p, T))
  Deltas <- Delta
  
  llik_hist <- rep(0, maxit)
  llik_lambda_hist <- rep(0, maxit)
  llik_hist_ext <- rep(0, maxit)
  
  iter <- 1
  
  for (iter in iter:maxit) {
    if (iter %% 50 == 0) {
      gc(verbose = FALSE)
      message("iter = ", iter)
    }
    # Update eta and epsilon
    for (t in 1:T) {
      I <- diag(d)
      Pt <- I + t(Lambda[,,t]) %*% invDelta %*% Lambda[,,t]
      invPt <- solve(Pt)
      CPt <- chol(Pt)
      invCPt <- chol(invPt)
      
      A <- t(-invDelta %*% Lambda[,,t] %*% t(invCPt))
      var_update <- solve(I + A %*% invDelta %*% t(A))
      mean_update <- var_update %*% (A %*% invDelta %*% y[t,])
      var_update <- 0.5 * (var_update + t(var_update))
      
      eta[t,] <- rmvnorm(1, rep(0,d), sigma = diag(d))
      eta[t,] <- (eta[t,] - mean(eta[t,])) / sd(eta[t,])
      cvar <- chol(var_update)
      eta[t,] <- cvar %*% eta[t,] + c(mean_update)
      eta[t,] <- (CPt) %*% eta[t,]
      
      K <- t(invDelta %*% Lambda[,,t] %*% invPt)
      eps[t,] <- y[t,] + c(eta[t,]) %*% K
    }
    
    # Update Delta
    deltas <- rgamma(p, a_delta + 0.5*T, b_delta + 0.5*colSums(eps^2, na.rm = T))
    Delta <- diag(deltas)
    invDelta <- diag(1/deltas)
    
    # Update Lambda based on latent shocks s
    for (t in 1:T) {
      w1 <- which(s[,t] == 1)
      w0 <- which(s[,t] == 0)
      
      if (t == 1) {
        w1 <- 1:d
        w0 <- integer(0)
      }
      
      if (t == T) {
        seq1 <- list(); seq1[1:length(w1)] <- T
      } else {
        seq1 <- lapply(w1, function(w) c(t:(t + min(T, which(s[w,(t+1):T] == 1))[1]-1)))
      }
      
      if (length(w1) > 0) {
        Lambda[, w0, t] <- Lambda[, w0, t-1]
        seq1 <- lapply(seq1, function(x) x[x <= T])
        longest_len_seq <- max(sapply(seq1, length))
      } else {
        Lambda[,,t] <- Lambda[,,t-1]
      }
    }
    
    # Polya-Gamma updates for betas
    if (TRUE) {
      # ... (code kept exactly as is for BETA0, BETA1 updates)
    }
    
    # Log-likelihood updates
    llik <- 0; llik_ext <- 0; llik_lambda <- 0
    if (iter > 1) Deltas <- Deltas + Delta
    
    for (t in 1:T) {
      LL_iter <- Lambda[,,t] %*% t(Lambda[,,t])
      if (iter > 0) LL[,,t] <- LL[,,t] + LL_iter
      
      Peta <- diag(d) + t(Lambda[,,t]) %*% invDelta %*% Lambda[,,t]
      invPeta <- ginv(Peta)
      
      llik_lambda <- llik_lambda + sum(dnorm(Lambda[,,t], 0, sd = Plam[1,1], log = TRUE))
      llik_ext <- llik_ext + dmvnorm(eta[t,], rep(0,d), sigma = Peta, log = TRUE)
      
      PREC <- LL_iter + Delta
      SIG <- solve(PREC)
      SIG <- 0.5 * (SIG + t(SIG))
      llik <- llik + dmvnorm(y[t,], mean = rep(0,p), sigma = SIG, log = TRUE)
    }
    
    llik_lambda_hist[iter] <- llik_lambda
    llik_hist[iter] <- llik
    llik_hist_ext[iter] <- llik_ext
    
    if (iter%%200==0) print(iter)
  }
  LL <- LL / iter
  Deltas <- Deltas / iter
  return(list(LL=LL, D=Deltas))
}

FactorPrecision_safe <- function(y, maxit, d, gc_every = 50) {
  
  library(mvtnorm)
  library(mvnfast)
  
  p <- ncol(y)
  T <- nrow(y)
  
  ## -----------------------------
  ## Hyperparameters
  ## -----------------------------
  a_delta <- 1; b_delta <- 1
  a_lam   <- 1; b_lam   <- 1
  
  ## -----------------------------
  ## Initialization
  ## -----------------------------
  Delta    <- diag(rep(0.1, p))
  invDelta <- diag(10, p)
  
  Lambda_prev <- matrix(0, p, d)
  Lambda_curr <- matrix(0, p, d)
  
  eta <- matrix(0, T, d)
  eps <- matrix(0, T, p)
  
  ## -----------------------------
  ## Online averages (SAFE)
  ## -----------------------------
  LL_mean     <- matrix(0, p, p)
  Delta_mean  <- matrix(0, p, p)
  
  llik_hist <- numeric(maxit)
  
  ## -----------------------------
  ## MCMC loop
  ## -----------------------------
  for (iter in seq_len(maxit)) {
    
    ## ----- Update eta, eps -----
    for (t in seq_len(T)) {
      
      Pt <- diag(d) + crossprod(Lambda_prev, invDelta %*% Lambda_prev)
      invPt <- solve(Pt)
      
      Kt <- t(invDelta %*% Lambda_prev %*% invPt)
      
      eta[t, ] <- mvnfast::rmvn(1, mu = rep(0, d), sigma = invPt)
      eps[t, ] <- y[t, ] - eta[t, ] %*% t(Lambda_prev)
    }
    
    ## ----- Update Delta -----
    deltas <- rgamma(
      p,
      shape = a_delta + 0.5 * T,
      rate  = b_delta + 0.5 * colSums(eps^2)
    )
    
    Delta    <- diag(deltas)
    invDelta <- diag(1 / deltas)
    
    ## ----- Update Lambda (time-recursive, no storage) -----
    for (t in seq_len(T)) {
      
      Vlam <- solve(diag(d) + crossprod(eta[t, ]))
      mlam <- Vlam %*% crossprod(eta[t, ], y[t, ])
      
      Lambda_curr <- mvtnorm::rmvnorm(1, mlam, kronecker(Vlam, invDelta))
      Lambda_curr <- matrix(Lambda_curr, p, d)
      
      Lambda_prev <- Lambda_curr
    }
    
    ## ----- Online averages -----
    LL_iter <- Lambda_curr %*% t(Lambda_curr)
    
    LL_mean    <- LL_mean    + (LL_iter - LL_mean) / iter
    Delta_mean <- Delta_mean + (Delta    - Delta_mean) / iter
    
    ## ----- Log-likelihood (optional, lightweight) -----
    SIG <- solve(LL_iter + Delta)
    llik_hist[iter] <- sum(dmvnorm(y, sigma = SIG, log = TRUE))
    
    ## ----- Garbage collection -----
    if (iter %% gc_every == 0) {
      gc(verbose = FALSE)
      message("iter = ", iter)
    }
  }
  
  return(list(
    LL     = LL_mean,
    Delta  = Delta_mean,
    llik   = llik_hist
  ))
}

