# Update the sigmas used to generate candidates
update_sigmas <- function(sigmas, samples_block, like_ratios, sigma_scalar) {
  if (all(apply(samples_block, 2, sd) < 1e-12)) {
    sigmas <- 0.9 * sigmas
  } else {
    med <- median(like_ratios, na.rm = TRUE)
    if (med < 0.5) {
      sigma_scalar <- 0.9 * sigma_scalar
    } else {
      sigma_scalar <- 1.1 * sigma_scalar
    }
    sigmas <- sigma_scalar * apply(samples_block, 2, sd)
  }
  return(list("sigmas" = sigmas, "sigma_scalar" = sigma_scalar))
}


# Get next candidate
get_candidate <- function(last_sample, mins, maxes, sigmas, dp) {
  while (TRUE) {
    # Generate a candidate
    tmp <- last_sample + sigmas * rnorm(length(sigmas))

    # Constrain sum(F_t) = 7
    idx <- (nrow(dp) + 1):(nrow(dp) + ncol(dp))
    tmp[idx] <- 7 * tmp[idx] / sum(tmp[idx])

    # Return candidate only if it falls within bounds
    if (all(tmp > mins) & all(tmp < maxes)) {
      return(tmp)
    }
  }
}


# Get log likelihood of given parameter estimates under bilinear model
get_log_likelihood <- function(params, dp) {
  # d' = R_sF_t + err
  R_s <- params[1:nrow(dp)]
  F_t <- params[(nrow(dp) + 1):(nrow(dp) + ncol(dp))]
  stdev <- tail(params, n = 1)
  dp_hat <- t(t(R_s)) %*% t(F_t)

  # Log likelihood
  log_p <- log(dnorm(dp, mean = dp_hat, sd = stdev))
  return(sum(log_p))
}


# Function to estimate the parameters using MCMC
mcmc <- function(d_prime, n_samples, n_thin, last_sample, last_sigmas, last_sigma_scalar) {
  # Use a shorter variable name
  dp <- d_prime
  
  # Count the number of parameters
  # Number subjects (R_s) + Number tasks (F_t) + one variance parameter
  n_params <- nrow(dp) + ncol(dp) + 1
  
  # Block size for updating the MCMC sigmas
  update_block_sz <- 5000
  
  # Prior on all parameters ~ Unif(-100, 100) except sigma ~ Unif(0, 100)
  mins <- c(rep(-100, n_params - 1), 0)
  maxes <- rep(100, n_params)
  
  # Default for last_sample
  if (missing(last_sample)) {
    last_sample <- c(rep(0, nrow(dp)), rep(1, ncol(dp)), 1)
  }
  
  # Default for last_sigmas
  if (missing(last_sigmas)) {
    sigmas <- rep(0.1, n_params)
  } else {
    sigmas <- last_sigmas
  }
  
  # Default for last_sigma_scalar
  if (missing(last_sigma_scalar)) {
    sigma_scalar <- 1
  } else {
    sigma_scalar <- last_sigma_scalar
  }
  
  
  # MCMC setup
  # Create arrays to store the *kept* parameter estimates and log likelihoods
  samples     <- matrix(, nrow = floor(n_samples / n_thin), ncol = n_params)
  candids     <- matrix(, nrow = floor(n_samples / n_thin), ncol = n_params)
  loglike_s   <- matrix(, nrow = floor(n_samples / n_thin), ncol = 1)
  loglike_c   <- matrix(, nrow = floor(n_samples / n_thin), ncol = 1)
  
  
  # Create arrays to temporarily store above + thinned samples
  tmp_samples     <- matrix(, nrow = update_block_sz, ncol = n_params)
  tmp_candids     <- matrix(, nrow = update_block_sz, ncol = n_params)
  tmp_loglike_s   <- matrix(, nrow = update_block_sz, ncol = 1)
  tmp_loglike_c   <- matrix(, nrow = update_block_sz, ncol = 1)
  tmp_like_ratios <- matrix(, nrow = update_block_sz, ncol = 1)
  
  # Initialize the chain
  tmp_samples[1, ] <- last_sample
  tmp_loglike_s[1] <- get_log_likelihood(tmp_samples[1, ], dp)
  
  
  # Run the MCMC
  for (k in 2:n_samples) {
    
    # To index temp arrays
    k_this <- ((k - 1) %% update_block_sz + 1)
    k_last <- if (k_this - 1 == 0) update_block_sz else k_this - 1 
    
    # Get a candidate
    candidate <- get_candidate(tmp_samples[k_last, ], mins, maxes, sigmas, dp)
    
    # Use likelihood ratio to determine the probability of accepting the new candidate
    tmp_candids[k_this, ] <- candidate
    tmp_loglike_c[k_this] <- get_log_likelihood(candidate, dp)
    tmp_like_ratios[k_this] <- exp(tmp_loglike_c[k_this] - tmp_loglike_s[k_last])
    if (runif(1) < tmp_like_ratios[k_this]) {
      # Accept the candidate
      tmp_samples[k_this, ] <- candidate
      tmp_loglike_s[k_this] <- tmp_loglike_c[k_this]
    } else {
      # Reject the candidate
      tmp_samples[k_this, ] <- tmp_samples[k_last, ]
      tmp_loglike_s[k_this] <- tmp_loglike_s[k_last]
    }
    
    # Update the sigmas at after each update block
    if (k_this == update_block_sz) {
      # Announce progress
      cat("\f")
      print(paste("Iteration", k, "of", n_samples))
      print(paste("Median LRatio:", median(tmp_like_ratios, na.rm = TRUE)))
      
      # Do the update
      update           <- update_sigmas(sigmas, tmp_samples, tmp_like_ratios, sigma_scalar)
      sigmas           <- update$sigmas
      sigma_scalar     <- update$sigma_scalar
    }
    
    # Keep every n_thin^th sample
    if (k %% n_thin == 0) {
      samples[k / n_thin, ] <- tmp_samples[k_this, ]
      candids[k / n_thin, ] <- tmp_candids[k_this, ]
      loglike_s[k / n_thin] <- tmp_loglike_s[k_this]
      loglike_c[k / n_thin] <- tmp_loglike_c[k_this]
    }
  }
  
  # Done
  return(list(
    "samples"         = samples,
    "candidates"      = candids,
    "loglike_samples" = loglike_s,
    "loglike_candids" = loglike_c,
    "sigmas"          = sigmas,
    "sigma_scalar"    = sigma_scalar
  ))
}