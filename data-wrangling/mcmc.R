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
mcmc <- function(d_prime, n_samples, last_sample, last_sigmas, last_sigma_scalar) {
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
  # Create arrays to store the parameter estimates and log likelihoods
  samples     <- matrix(, nrow = n_samples, ncol = n_params)
  candids     <- matrix(, nrow = n_samples, ncol = n_params)
  loglike_s   <- matrix(, nrow = n_samples, ncol = 1)
  loglike_c   <- matrix(, nrow = n_samples, ncol = 1)
  like_ratios <- matrix(, nrow = update_block_sz, ncol = 1)

  # Initialize the chain
  samples[1, ] <- last_sample
  loglike_s[1] <- get_log_likelihood(samples[1, ], dp)
  update_block_ctr <- 0


  # Run the MCMC
  for (k in 2:n_samples) {
    # Get a candidate
    update_block_ctr <- update_block_ctr + 1
    candidate <- get_candidate(samples[k - 1, ], mins, maxes, sigmas, dp)

    # Use likelihood ratio to determine the probability of accepting the new candidate
    candids[k, ] <- candidate
    loglike_c[k] <- get_log_likelihood(candidate, dp)
    like_ratios[update_block_ctr] <- exp(loglike_c[k] - loglike_s[k - 1])
    if (runif(1) < like_ratios[update_block_ctr]) {
      # Accept the candidate
      samples[k, ] <- candidate
      loglike_s[k] <- loglike_c[k]
    } else {
      # Reject the candidate
      samples[k, ] <- samples[k - 1, ]
      loglike_s[k] <- loglike_s[k - 1]
    }

    # Update the sigmas at after each update block
    if (update_block_ctr == update_block_sz) {
      # Announce progress
      cat("\f")
      print(paste("Iteration", k, "of", n_samples))
      print(paste("Median LRatio:", median(like_ratios, na.rm = TRUE)))

      # Do the update
      update           <- update_sigmas(sigmas, samples[(k - update_block_sz + 1):k, ], like_ratios, sigma_scalar)
      sigmas           <- update$sigmas
      sigma_scalar     <- update$sigma_scalar
      update_block_ctr <- 0
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