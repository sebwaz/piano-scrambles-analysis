###
### Fit bilinear model using SVD, test models with additional components ---------------------------
###

bilinear_svd <- function(dp) {
  # Use svd instead of mcmc to get the point estimates for the bilinear model
  dp_decomp <- svd(dp)
  
  # Each row of D %*% t(V)), normalized, is a set of F estimates
  F_svd <- diag(dp_decomp$d) %*% t(dp_decomp$v)
  scale_factor <- rowMeans(F_svd)
  for (r in 1:6) {
    F_svd[r, ] <- F_svd[r, ] / scale_factor[r]
  }
  
  # Each column of U, normalized is a set of R estimates
  R_svd <- dp_decomp$u
  for (k in 1:6) {
    R_svd[, k] <- R_svd[, k] * scale_factor[k]
  }
  
  # Perform the F-tests for each model relative to the previous model
  grand_mean <- mean(as.vector(dp))
  sst        <- sum(as.vector(dp - grand_mean) ^ 2) # Total squared error in the data
  
  dfx   <- c()   # 'Extra' df, i.e., df current model minus df previous model
  dfres <- c()   # Total df minus sum of df's of previous models
  ssr   <- c()   # Sum of squares explained by current model
  ssx   <- c()   # 'Extra' sum of squares explained by current model relative to previous model
  ssres <- c()
  msx   <- c()
  msres <- c()
  fstat <- c()
  pvals <- c()
  
  for (k in 1:6) {
    # degrees of freedom calculations
    dfx[k]   <- length(R_svd[, k]) + length(F_svd[k, ]) - (2 * k - 1)
    dfres[k] <- length(R_svd[, k]) * length(F_svd[k, ]) - sum(dfx)
    
    # prediction error calculations
    #   Note: it seems F_svd[k, ] gives the k^th row's data, organized as a *column*
    #   so the transpose in t(F_svd[1:1, ]) is used to account for R behavior, not linear algebra
    #   However, if you index multiple rows/columns, orientation is preserved :upside-down smiley face:
    if (k == 1) {
      err <- R_svd[, k] %*% t(F_svd[k, ]) - dp
    } else {
      err <- R_svd[, 1:k] %*% F_svd[1:k, ] - dp
    }
    
    # sum-of-squares calculations
    ssres[k] <- sum(as.vector(err) ^ 2)
    ssr[k]   <- sst - ssres[k]
    if (k > 1) {
      ssx[k] <- ssr[k] - ssr[k - 1]
    } else {
      ssx[k] <- ssr[k]
    }
    
    # mean-squared-error calculations
    msx[k]   <- ssx[k] / dfx[k]
    msres[k] <- ssres[k] / dfres[k]
    
    # F-statistic calculation
    fstat[k] <- msx[k] / msres[k]
    
    # p-values
    pvals[k] <- 1 - pf(fstat[k], dfx[k], dfres[k])
  }
  
  # R-squared of each model
  Rsq <- ssr / sst
  
  # Done
  return(list(
    "R"      = R_svd,
    "F"      = F_svd,
    "pvals"  = pvals,
    "Rsq"    = Rsq,
    "fstat"  = fstat,
    "ssr"    = ssr,
    "ssx"    = ssx,
    "ssres"  = ssres,
    "dfx"    = dfx,
    "dfres"  = dfres
  ))
}