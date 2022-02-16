library(lubridate)
library(tidyverse)
library(ggplot2)
library(audio)
library(seewave)
library(reshape2)

data_path <- here::here("data", "2022_01_18")
source(here::here("data-wrangling", "sdt.R"))

# Get the list of CSV files
csv_list <- list.files(path = data_path, pattern = "*.csv")

# Initialize a data.frame to contain the summary for subjects
all_data = list()

for (k in 1:length(csv_list)) {
  
  ##### Load each subject's data, one-by-one #####
  subj_data <- read.csv(paste(data_path, "/", csv_list[k], sep = ""))
  
  # Start time is encoded in the file name (not quite start time; initial entry time)
  time_start <- as.numeric(strsplit(strsplit(csv_list[k], "_")[[1]][2],  ".", fixed = TRUE)[[1]][1])
  
  # Omit the training trials
  subj_data <- subj_data[subj_data$training_bool == "false", ]

  # Also separate out the headphone trials
  subj_headphonecheck <- subj_data[subj_data$block == -1, ]
  subj_data           <- subj_data[subj_data$block != -1, ]
  
  
  
  ##### Validate the data #####
  # Check that all minor trials and all major trials have the same note histogram
  for (r in 1:nrow(subj_data)) {
    if (!is.nan(subj_data$scramble_type[r]) & !is.nan(subj_data$cur_trial[r])) {
      
      # Check signal trials
      if (subj_data$scramble_type[r] == 1) {
        
        # Note that conditions 2 and 5 are fast and therefore have more pips
        if (subj_data$condition[r] == 2 || subj_data$condition[r] == 5) {
          if (any(tabulate(as.numeric(unlist(strsplit(subj_data$seq[r], split = ","))) + 1) != c(8, 0, 0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 8)))
            stop(paste("Bad stimulus found for subject", k, "row", r))
        } else {
          if (any(tabulate(as.numeric(unlist(strsplit(subj_data$seq[r], split = ","))) + 1) != c(3, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 3)))
            stop(paste("Bad stimulus found for subject", k, "row", r))
        }
        
      } else if (subj_data$scramble_type[r] == 2) {
        
         # Note that conditions 2 and 5 are fast and therefore have more pips
         if (subj_data$condition[r] == 2 || subj_data$condition[r] == 5) {
           if (any(tabulate(as.numeric(unlist(strsplit(subj_data$seq[r], split = ","))) + 1) != c(8, 0, 0, 0, 8, 0, 0, 8, 0, 0, 0, 0, 8)))
             stop(paste("Bad stimulus found for subject", k, "row", r))
         } else {
           if (any(tabulate(as.numeric(unlist(strsplit(subj_data$seq[r], split = ","))) + 1) != c(3, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 3)))
             stop(paste("Bad stimulus found for subject", k, "row", r))
         }
      }
    }
  }
  
  # Check that there's an even split of S/N trials in each block for first 10
  mask <- (subj_data$cur_trial <= 10) & !is.nan(subj_data$cur_trial)
  for (block in 0:11) {
    if (mean(subj_data$scramble_type[mask & (subj_data$block == block)] - 1) != 0.5) {
      stop(paste("Uneven number of S/N trials detected in training segment of block", block))
    }
  }
  
  # Check that there's an even split of S/N trials in each block for last 30
  mask <- (subj_data$cur_trial > 10) & (subj_data$cur_trial < 999) & !is.nan(subj_data$cur_trial)
  for (block in 0:11) {
    if (mean(subj_data$scramble_type[mask & (subj_data$block == block)] - 1) != 0.5) {
      stop(paste("Uneven number of S/N trials detected in testing segment of block", block))
    }
  }
  
  
  
  ##### Analyze the data #####
  
  # Headphone check
  n_attempts <- nrow(subj_headphonecheck) / 6
  first_score <- sum(subj_headphonecheck$scramble_type[1:6]  == subj_headphonecheck$response[1:6])
  final_score <- if (n_attempts == 2) sum(subj_headphonecheck$scramble_type[7:12] == subj_headphonecheck$response[7:12]) else sum(subj_headphonecheck$scramble_type[1:6]  == subj_headphonecheck$response[1:6])
  passed <- final_score >= 5
  
  # Get dprime by block
  #  0 : slow, piano, G5
  #  1 : slow, pure,  G5
  #  2 : fast, pure,  G5
  #  3 : slow, piano, C4
  #  4 : slow, pure,  C4
  #  5 : fast, pure,  C4
  nHit_nMiss_nCR_nFA <- matrix(, nrow = 12, ncol = 4)
  dp <- matrix(, nrow = 1, ncol = 12)
  bias <- matrix(, nrow = 1, ncol = 12)
  for (block in 0:11) {
    res <- sdt(subj_data$scramble_type[mask & (subj_data$block == block)], subj_data$response[mask & (subj_data$block == block)])
    nHit_nMiss_nCR_nFA[block + 1, ] <- c(res$hit, res$miss, res$cr, res$fa)
    dp[block + 1] <- res$dprime
    bias[block + 1] <- res$bias
  }


  # Get the condition order
  cond_order <- as.numeric(unlist(strsplit(subj_data$cond_order[1], split=",")))

  # Also get dprime, bias, and % correct, marginalizing over blocks of the same condition
  dp_cond <- matrix(, nrow = 1, ncol = 6)
  bias_cond <- matrix(, nrow = 1, ncol = 6)
  p_cond <- matrix(, nrow = 1, ncol = 6)
  for (cond in 1:6) {
    bool_cond <- which(cond == cond_order)
    stats <- sdt(nhit = sum(nHit_nMiss_nCR_nFA[bool_cond, 1]), nmiss = sum(nHit_nMiss_nCR_nFA[bool_cond, 2]),
                 ncr  = sum(nHit_nMiss_nCR_nFA[bool_cond, 3]), nfa   = sum(nHit_nMiss_nCR_nFA[bool_cond, 4]))
    dp_cond[cond] <- stats$dprime
    bias_cond[cond] <- stats$bias
    p_cond[cond] <- stats$accuracy
  }

  # Get the finishing time
  time_finish <- tail(subj_data, 1)$resp_time

  # Put all the data in a list
  all_data[[k]] <- list(
    "inits" = subj_data$subj_id[1],
    "n_attempts_heaphone" = n_attempts,
    "first_score_heaphone" = first_score,
    "final_score_headphone" = final_score,
    "passed" = passed,
    "cond_order" = cond_order,
    "nHit_nMiss_nCR_nFA" = nHit_nMiss_nCR_nFA,
    "dprime" = dp,
    "bias" = bias,
    "dprime_cond" = dp_cond,
    "bias_cond" = bias_cond,
    "p_cond" = p_cond,
    "native_lang" = if (is.na(subj_data$lang_other[1])) subj_data$lang[1] else subj_data$lang_other[1],
    "years_train" = if (subj_data$years_train[1] == "none") 0 else subj_data$years_train[1],
    "latin_row" = subj_data$latin_row[1],
    "device_samp_hz" =  subj_data$device_samp_hz[1],
    "time_start" = time_start,
    "time_finish" = time_finish,
    "mins_to_finish" = (time_finish - time_start)/60000
  )
}



###
### Rearrange the data into a tibble ------------------------------
###

# Get N
n  <- length(all_data)

# Matrices to store the data we are interested in
dp <- matrix(, nrow = n, ncol = 6)  # dprimes
bi <- matrix(, nrow = n, ncol = 6)  # biases
pc <- matrix(, nrow = n, ncol = 6)  # proportion correct
tc <- matrix(, nrow = n, ncol = 1)  # time to complete
fs <- matrix(, nrow = n, ncol = 1)  # system sampling rate
nl <- character(n)                  # native language
yt <- matrix(, nrow = n, ncol = 1)  # years of music training
ph <- matrix(, nrow = n, ncol = 1)  # passed headphone check?

# Fill the matrices 
for (k in 1:n) {
  dp[k, ] <- all_data[[k]]$dprime_cond
  bi[k, ] <- all_data[[k]]$bias_cond
  pc[k, ] <- all_data[[k]]$p_cond
  tc[k]   <- all_data[[k]]$mins_to_finish
  fs[k]   <- all_data[[k]]$device_samp_hz
  nl[k]   <- all_data[[k]]$native_lang
  
  if (all_data[[k]]$years_train == "None") {
    yt[k] <- 0
  } else {
    yt[k] <- all_data[[k]]$years_train
  }
  
  if (all_data[[k]]$passed) {
    ph[k] <- 1
  } else {
    ph[k] <- 0
  }
}

# Label the data in the 6-column matrices (one column per condition)
labs <- c("0", "1", "2", "3", "4", "5")
labs <- paste("val", labs, sep = "_")
colnames(dp) <- paste("dprime", labs, sep = "_")
colnames(bi) <- paste("bias", labs, sep = "_")
colnames(pc) <- paste("pcorrect", labs, sep = "_")



# Create a tibble X of the data
tmp <- cbind(1:n, nl)
colnames(tmp)[1] <- "subj_num"
colnames(tmp)[2] <- "native_lang"
tmp <- as_tibble(tmp) %>%
  mutate_at(vars("subj_num"), as.numeric)

X <- cbind(1:n, ph, tc, fs, yt, dp, bi, pc)
colnames(X)[1] <- "subj_num"
colnames(X)[2] <- "passed"
colnames(X)[3] <- "mins_complete"
colnames(X)[4] <- "device_samp_hz"
colnames(X)[5] <- "years_train"
X <- as_tibble(X) %>%
  left_join(tmp) %>%
  pivot_longer(cols = contains("val"),
               names_to = c("column_name", ".value", "condition"),
               names_sep = "_",
               values_to = "column_value") %>%
  pivot_wider(names_from = "column_name",
              values_from = "val")



###
### Some preliminary analyses -----------------
###

# Look at N, N passed, and % of participants who passed the headphone check
n_passed <- sum(ph)
n_passed
n
n_passed / n


# Fit a linear model to measure the effect of pass/no pass, and whether it modulates any of the condition-specific effects
#  A: it does not affect results
summary(lm(dprime ~ passed + factor(condition) + factor(condition) * passed - 1, data = X))


# Fit a linear model to measure the effect of the system's sampling rate
#  A: It has no effect
summary(lm(dprime ~ factor(condition) + factor(condition) * factor(device_samp_hz) - 1, data = X))


# Plot the d's in each condition, using color to separate pass/no pass
ggplot(X, aes(x = condition, y = dprime, color = passed)) +
  geom_point()


# Spaghetti plots of dprimes over blocks
plot(NaN, NaN, xlim = c(1, 12), ylim = c(-1, 5), xlab = "Block number", ylab = "dprime")
for (k in 1:length(all_data)) {
  lines(1:12, all_data[[k]]$dprime)
}


# Spaghetti plot of dprime over conditions
plot(NaN, NaN, xlim = c(0, 5), ylim = c(-1, 5), xlab = "Condition", ylab = "dprime")
for (k in 1:length(all_data)) {
  lines(0:5, all_data[[k]]$dprime_cond)
}


# Look at how long people are taking
mins_to_finish <- matrix(, ncol = 1, nrow = length(all_data))
for (k in 1:length(all_data)) {
  mins_to_finish[k] <- all_data[[k]]$mins_to_finish
}
summary(mins_to_finish)


# Get the Latin rows and check numbers in each ordering
lr <- c()
for (k in 1:length(all_data)) {
  lr <- c(lr, all_data[[k]]$latin_row)
}
sort(lr)


# Plot the means and standard deviations
X %>%
  group_by(condition) %>%
  summarize(mean = mean(dprime), ci = 1.96 * sd(dprime)/sqrt(n)) %>%
  ggplot(aes(x = condition, y = mean)) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_point()


# Paired-samples t-test comparing performance on
#  the standard tone-scramble task vs.
#  the task in which notes are drawn from C4 - C5, and have piano timbre
#  0 : slow, piano, G5
#  1 : slow, pure,  G5
#  2 : fast, pure,  G5
#  3 : slow, piano, C4
#  4 : slow, pure,  C4
#  5 : fast, pure,  C4
# Remember that R is 1 indexed
t.test(dp[, 3], dp[, 4], paired = TRUE, alternative = "two.sided")

# Plot comparing d' on conditions 2 and 3
X %>%
  filter(condition == 2 | condition == 3) %>%
  pivot_wider(id_cols = subj_num, names_from = condition, values_from = dprime) %>%
  ggplot(aes(x = `2`, y = `3`)) +
  geom_point() +
  geom_abline() +
  labs(x = "d' (pure tones, 32 pips, G5)",
       y = "d' (piano tones, 12 pips, C4)") +
  xlim(c(-1, 5)) +
  ylim(c(-1, 5)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))



###
### More plots ---------------------
###

# Condition labeller
#  0 : slow, piano, G5
#  1 : slow, pure,  G5
#  2 : fast, pure,  G5
#  3 : slow, piano, C4
#  4 : slow, pure,  C4
#  5 : fast, pure,  C4
cond_names <- list(
  '0'="slow, piano, G5",
  '1'="slow, pure,  G5",
  '2'="fast, pure,  G5",
  '3'="slow, piano, C4",
  '4'="slow, pure,  C4",
  '5'="fast, pure,  C4"
)
cond_labeller <- function(variable, value){
  return(cond_names[value])
}

# Plot counter
ctr <- 0

## Reproduce basic finding (in d's)
ggplot(data = X, aes(dprime)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue3", size = 1) +
  facet_grid(. ~ condition, labeller = cond_labeller) +
  scale_x_continuous(breaks = seq(-1, 4, 1)) +
  labs(x = expression("Sensitivity (d')"),
       title = expression(~bold("Bimodality in sensitivity?"))
  ) +
  ylim(c(0, 50)) +
  theme(text = element_text(size = 16))
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



## Reproduce basic finding (in % correct)
ggplot(data = X, aes(pcorrect)) +
  geom_histogram(binwidth = 0.1, color = "black", fill = "steelblue3", size = 1) +
  facet_grid(. ~ condition, labeller = cond_labeller) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  labs(x = expression("% correct"),
       title = expression(~bold("Bimodality in percent correct"))
  ) +
  ylim(c(0, 50)) +
  theme(text = element_text(size = 16))
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



## Plot % correct overall against time to complete
X %>%
  group_by(subj_num) %>%
  summarize(p_corr_overall = mean(pcorrect), mins_complete = max(mins_complete), .groups = 'drop') %>%
  ggplot() +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", size = 1) +
  geom_point(aes(y = p_corr_overall, x = mins_complete), color = "springgreen4", size = 3) +
  scale_x_continuous(breaks = seq(20, 200, 10), limits = c(20, 200)) +
  labs(x = "Minutes to complete",
       y = "% correct overall",
       title = expression(~bold("At least two participants appear to have rushed (one extremely long observation omitted)"))
  ) +
  theme(text = element_text(size = 12))
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



# Use ggpairs to plot the dprimes pairwise
# Thanks to this stackoverflow page for how to quickly add the x=y line to each subplot:
#  https://stackoverflow.com/questions/26780332/how-to-add-geom-abline-to-ggpairs
library(GGally)
hdl <- X %>%
  pivot_wider(id_cols = subj_num,
              names_from = condition,
              names_prefix = "cond",
              values_from = dprime) %>%
  ggpairs(columns = 2:7)
for (i in 2:hdl$nrow) {
  for (j in 1:(i-1)) {
    hdl[i,j] <- hdl[i,j] + geom_abline(intercept=0,slope=1)
  }
}
hdl
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



# Plot differences in condition 3 vs. condition 2 against musical training
X %>%
  filter(condition == 2 | condition == 3) %>%
  pivot_wider(id_cols = c("subj_num", "years_train", "native_lang"), names_from = condition, values_from = dprime) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1) +
  geom_point(aes(y = `3` - `2`, x = years_train, color = native_lang), size = 2) +
  # scale_color_manual(values = c("springgreen3", "violetred3", "violetred3", "midnightblue", "springgreen3", "springgreen3", "springgreen3", "springgreen3", "springgreen3")) +
  labs(x = "Years of musical training",
       y = expression("d'"[slowPianoC4]*" - d'"[fastPureG5]),
       title = expression(~bold("No clear relationship between training and effect of piano (impairs highly trained Mandarin speakers?)")),
       color = "Native language"
  ) +
  theme(text = element_text(size = 12))
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



# Plot biases
X %>%
  filter(condition == 2 | condition == 3) %>%
  pivot_wider(id_cols = subj_num, names_from = condition, values_from = bias) %>%
  ggplot() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1) +
  geom_point(aes(y = `3`, x = `2`), size = 2) +
  labs(x =  expression("bias"[fastPureG5]),
       y =  expression("bias"[slowPianoC4]),
       title = expression(~bold("No systematic biases, one bad performer"))
  ) +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  theme(text = element_text(size = 16))
ctr <- ctr + 1
ggsave(here::here("data-wrangling", paste(today(), "-",ctr, ".png", sep = "")), plot = last_plot(), width = 10, height = 6, units = "in")



###
### Fit the bilinear model ---------------------
###

# Run MCMC
source(here::here("data-wrangling", "mcmc.R"))
if (!file.exists(here::here("data", "posterior_samples.rds"))) {
  burnin    <- mcmc(dp, 10000000, 1000)
  posterior <- mcmc(dp, 10000000, 1000, tail(burnin$samples, 1), burnin$sigmas, burnin$sigma_scalar)

  # Save the samples just taken
  save(burnin, posterior,
       file = here::here("data", "posterior_samples.rds"))
} else {
  load(here::here("data", "posterior_samples.rds"))
}



# Process samples a little bit to get credible intervals
tbl <- as_tibble(posterior$samples)

posterior_F <- tbl %>%
  pivot_longer(cols = everything(), names_to = "parameter") %>%
  filter(parameter %in% paste("F", 1:6, sep = "_"))

posterior_F_CI <- posterior_F %>%
  group_by(parameter) %>%
  summarize(Median = median(value),
            CI_L = quantile(value, 0.025),
            CI_U = quantile(value, 0.975))
  
# Plot the posterior estimates as violin plot
ggplot(data = posterior_F,
       aes(x = parameter, y = value)) +
geom_violin(fill = "gray80", color = "gray70") +
geom_errorbar(data = posterior_F_CI,
              mapping = aes(x = parameter, ymin = CI_L, ymax = CI_U),
              width = 0.2,
              inherit.aes = FALSE) +

# Add median point
geom_point(data = posterior_F_CI,
           mapping = aes(x = parameter, y = Median),
           inherit.aes = FALSE) +

# Add labels
labs(x = expression(paste("Task ", italic(t))),
     y = expression("F"[t])
) +

# Set y-limits, background color, text size
ylim(0, 1.4) +
theme_bw() +
theme(text = element_text(size = 12))



# Look at the shape of the distributions
tbl %>%
  select(contains("F")) %>%
  ggpairs(1:6)



# Look at the difference between fast-pure-G5 and slow-piano-C4
tbl %>%
  select(contains("F")) %>%
  mutate(diff = F_4 - F_3) %>%
  summarize(Median = median(diff),
            CI_L = quantile(diff, 0.025),
            CI_U = quantile(diff, 0.975))


  
# Look at the shape of the posterior difference
tbl %>%
  select(contains("F")) %>%
  mutate(diff = F_4 - F_3) %>%  
  ggplot(aes(diff)) +
  geom_histogram(binwidth = 0.02, color = "black", fill = "steelblue3", size = 1)



# Look at the distribution of R
tbl %>%
  select(contains("R")) %>%
  pivot_longer(cols = everything(), names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(Median = median(value),
            CI_L = quantile(value, 0.025),
            CI_U = quantile(value, 0.975)) %>%
  ggplot(aes(Median)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "steelblue3", size = 1) +
  labs(x = expression("R"[s]),
       title = expression(~bold("Distribution of median posterior R estimates"))
  )
  


###
### Calculate R^2 of the bilinear model ---------------------
###

# Take posterior medians as parameter point estimate
median_posterior_samples <- tbl %>%
  pivot_longer(cols = everything(), names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(Median = median(value))

# F estimates (take extra pains to ensure they come out in the right order)
F_hat <- median_posterior_samples %>%
  filter(grepl("F",parameter)) %>%
  mutate(cond_num = as.numeric(str_replace(parameter, "F_", "")) - 1) %>%
  arrange(cond_num) %>%
  select(Median)

# R estimates (take extra pains to ensure they come out in the right order)
R_hat <- median_posterior_samples %>%
  filter(grepl("R",parameter)) %>%
  mutate(subj_num = as.numeric(str_replace(parameter, "R_", ""))) %>%
  arrange(subj_num) %>%
  select(Median)

# Predicted d's
dp_pred <- as.matrix(R_hat) %*% as.matrix(t(F_hat))

# Calculate R^2
grand_mean <- mean(as.vector(dp))
ss_tot     <- sum(as.vector(dp - grand_mean) ^ 2)
ss_res     <- sum(as.vector(dp -    dp_pred) ^ 2)
Rsq <- (ss_tot - ss_res) / ss_tot
Rsq

# Create a matrix with predicted and observed d's to plot
colnames(dp_pred) <- paste("predicted", labs, sep = "_")
dp_tbl <- cbind(1:n, dp, dp_pred)
colnames(dp_tbl)[1] <- "subj_num"
dp_tbl <- as_tibble(dp_tbl) %>%
  pivot_longer(cols = contains("val"),
               names_to = c("column_name", ".value", "condition"),
               names_sep = "_",
               values_to = "column_value") %>%
  pivot_wider(names_from = "column_name",
              values_from = "val")

# Plot observed vs. predicted d's
dp_tbl %>%
  ggplot(aes(x = dprime, y = predicted, color)) +
  geom_point(aes(color = condition), size = 2) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(-1, 5)) +
  ylim(c(-1, 5)) +
  labs(x = "observed d'",
       y = "predicted d'",
       title = "bilinear model fit"
  )



###
### Do F-tests of the bilinear model and models with additional components ---------------------------
###

source(here::here("data-wrangling", "bilinear_svd.R"))
mdl_svd <- bilinear_svd(dp)

# Verify that what we get from SVD is reasonably close to the point estimates from MCMC (yes)
summary(mdl_svd$F[1, ] - F_hat)
summary(mdl_svd$R[, 1] - R_hat)



###
### Try to fit a mixture of Gaussians to the d's --------------------------
###

# Try to categorize listeners into high- and low-performers using a Gaussian mixture
library(mclust)
fit  <- list()
crit <- list()
par(mfrow = c(1, 6))
for (k in 1:6) {
  # Fit the mixture of two Gaussians, add to the plot
  fit[[k]] <- Mclust(dp[, k], G = 2)
  plot.Mclust(fit[[k]], what = "density")
  #print(fit[[k]]$parameters$mean)
  
  # Find the approximate intersection of the two Gaussian densities, i.e., the point where the MLE of category changes
  # (Just use the categorization MClust provides, no genuine root finding)
  # (Note that Gaussians will technically have two intersections if unequal variance, so very low individuals may end up in the high group w/o this step)
  lo <- max(dp[fit[[k]]$classification == 1 & dp[, k] < fit[[k]]$parameters$mean[2], k])
  hi <- min(dp[fit[[k]]$classification == 2 & dp[, k] > fit[[k]]$parameters$mean[1], k])
  cr <- mean(c(lo, hi))
  
  crit[[k]] <- list(
    "lo" = lo,
    "hi" = hi,
    "crit" = cr
  )
}



