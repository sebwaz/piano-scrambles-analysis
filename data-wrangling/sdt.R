sdt <- function(answer, response, nhit = NaN, nmiss = NaN, ncr = NaN, nfa = NaN) {
  if (is.nan(nhit) | is.nan(nmiss) | is.nan(ncr) | is.nan(nfa)) {
    # Labels for the trial answers
    noise  <- 1
    signal <- 2
  
    # Number total of trials, number of signal trials, and number of noise trials
    n_trials <- length(answer)
    n_noise  <- sum(answer == noise)
    n_signal <- sum(answer == signal)
  
    # Get the number of hits, misses, correct rejections, and false alarms
    hit  <- sum((answer == response) * (answer == signal))
    miss <- sum((answer != response) * (answer == signal))
    cr   <- sum((answer == response) * (answer == noise))
    fa   <- sum((answer != response) * (answer == noise))
  } else {
    # Number total of trials, number of signal trials, and number of noise trials
    n_trials <- sum(nhit, nmiss, ncr, nfa)
    n_noise  <- sum(ncr, nfa)
    n_signal <- sum(nhit, nmiss)

    # Get the number of hits, misses, correct rejections, and false alarms
    hit  <- nhit
    miss <- nmiss
    cr   <- ncr
    fa   <- nfa
  }

  # Calculate p(hit) and p(false alarm)
  p_hit <- hit / n_signal
  p_fa <- fa / n_noise
  
  # Adjustments based on Macmillan & Creelman, 2nd edition (2005), page 8
  if (p_hit == 1) {
    p_hit <- (n_signal - 0.5) / n_signal
  }

  if (p_fa == 1) {
    p_fa <- (n_noise - 0.5) / n_noise
  }

  if (p_hit == 0) {
    p_hit <- 0.5 / n_signal
  }

  if (p_fa == 0) {
    p_fa <- 0.5 / n_noise
  }

  # Calculate d-prime
  d_prime <- qnorm(p_hit) - qnorm(p_fa)

  # Calculate bias
  bias <- qnorm(1 - p_fa) - d_prime / 2

  # Calculate the overall accuracy
  accuracy <- (hit + cr) / n_trials

  # Done
  return(list(
    "dprime"   = d_prime,
    "bias"     = bias,
    "accuracy" = accuracy,
    "p_hit"    = p_hit,
    "p_fa"     = p_fa,
    "hit"      = hit,
    "miss"     = miss,
    "cr"       = cr,
    "fa"       = fa
  ))
}