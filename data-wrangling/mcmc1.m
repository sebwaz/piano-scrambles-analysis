% Load the d's
dp = readmatrix("dprimes.csv");

% Do the burn-in
[samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(dp, 10000000, 1000, [], [], []);

% Do the final posterior samples
[samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(dp, 10000000, 1000, samples(1, :), sigmas, sigma_scalar);

% Save results
save("mcmc1_samples.mat");
writematrix(samples, "mcmc1_samples.csv");

% Update the sigmas used to generate candidates
function [sigmas, sigma_scalar] = update_sigmas(sigmas, samples_block, like_ratios, sigma_scalar)
if (all(std(samples_block) < 1e-12))
    sigmas = 0.9 * sigmas;
else
    med = median(like_ratios, 'omitnan');
    if (med < 0.5)
        sigma_scalar = 0.9 * sigma_scalar;
    else
        sigma_scalar = 1.1 * sigma_scalar;
    end
    sigmas = sigma_scalar * std(samples_block);
end
end


% Get next candidate
function tmp = get_candidate(last_sample, mins, maxes, sigmas, dp)
while true
    % Generate a candidate
    tmp = last_sample + sigmas .* randn(size(sigmas));
    
    % Constrain sum(F) = length(F)
    idx = (size(dp, 1) + 1):(size(dp, 1) + size(dp, 2));
    tmp(idx) = size(dp, 2) * tmp(idx) / sum(tmp(idx));
    
    % Return candidate only if it falls within bounds
    if (all(tmp > mins) && all(tmp < maxes))
        break
    end
end
end


% Get log likelihood of given parameter estimates under bilinear model
function sum_log_p = get_log_likelihood(params, dp)
% d' = R_sF_t + err
R_s = params(1:size(dp, 1));
F_t = params((size(dp, 1) + 1):(size(dp, 1) + size(dp, 2)));
stdev = params(end);
dp_hat = R_s' * F_t;

% Log likelihood
log_p = log(normpdf(dp, dp_hat, stdev));
sum_log_p = sum(log_p(:));
end


% Function to estimate the parameters using MCMC
function [samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(d_prime, n_samples, n_thin, last_sample, last_sigmas, last_sigma_scalar)
% Use a shorter variable name
dp = d_prime;

% Count the number of parameters
% Number subjects (R_s) + Number tasks (F_t) + one variance parameter
n_params = size(dp, 1) + size(dp, 2) + 1;

% Block size for updating the MCMC sigmas
update_block_sz = 5000;

% Prior on all parameters ~ Unif(-100, 100) except sigma ~ Unif(0, 100)
mins  = [-100 * ones(1, n_params - 1), 0];
maxes = 100 * ones(1, n_params);

% Default for last_sample
if isempty(last_sample)
    last_sample = [zeros(1, size(dp, 1)), ones(1, size(dp, 2)), 1];
end

% Default for last_sigmas
if isempty(last_sigmas)
    sigmas = 0.1 * ones(1, n_params);
else
    sigmas = last_sigmas;
end

% Default for last_sigma_scalar
if isempty(last_sigma_scalar)
    sigma_scalar = 1;
else
    sigma_scalar = last_sigma_scalar;
end


% MCMC setup
% Create arrays to store the *kept* parameter estimates and log likelihoods
samples     = nan(floor(n_samples / n_thin), n_params);
candids     = nan(floor(n_samples / n_thin), n_params);
loglike_s   = nan(floor(n_samples / n_thin), 1);
loglike_c   = nan(floor(n_samples / n_thin), 1);


% Create arrays to temporarily store above + thinned samples
tmp_samples     = nan(update_block_sz, n_params);
tmp_candids     = nan(update_block_sz, n_params);
tmp_loglike_s   = nan(update_block_sz, 1);
tmp_loglike_c   = nan(update_block_sz, 1);
tmp_like_ratios = nan(update_block_sz, 1);

% Initialize the chain
tmp_samples(1, :) = last_sample;
tmp_loglike_s(1)  = get_log_likelihood(tmp_samples(1, :), dp);


% Run the MCMC
for k = 2:n_samples
    
    % To index temp arrays
    k_this = mod(k - 1, update_block_sz) + 1;
    if (k_this - 1 == 0)
        k_last = update_block_sz;
    else
        k_last = k_this - 1;
    end
    
    % Get a candidate
    candidate = get_candidate(tmp_samples(k_last, :), mins, maxes, sigmas, dp);
    
    % Use likelihood ratio to determine the probability of accepting the new candidate
    tmp_candids(k_this, :)  = candidate;
    tmp_loglike_c(k_this)   = get_log_likelihood(candidate, dp);
    tmp_like_ratios(k_this) = exp(tmp_loglike_c(k_this) - tmp_loglike_s(k_last));
    if (rand < tmp_like_ratios(k_this))
        % Accept the candidate
        tmp_samples(k_this, :) = candidate;
        tmp_loglike_s(k_this)  = tmp_loglike_c(k_this);
    else
        % Reject the candidate
        tmp_samples(k_this, :) = tmp_samples(k_last, :);
        tmp_loglike_s(k_this)  = tmp_loglike_s(k_last);
    end
    
    % Update the sigmas at after each update block
    if (k_this == update_block_sz)
        % Announce progress
        clc;
        fprintf("Iteration %d of %d\n", k, n_samples);
        fprintf("Median LRatio: %10.9f\n", median(tmp_like_ratios, 'omitnan'));
        
        % Do the update
        [sigmas, sigma_scalar] = update_sigmas(sigmas, tmp_samples, tmp_like_ratios, sigma_scalar);
    end
    
    % Keep every n_thin^th sample
    if (mod(k, n_thin) == 0)
        samples(k / n_thin, :) = tmp_samples(k_this, :);
        candids(k / n_thin, :) = tmp_candids(k_this, :);
        loglike_s(k / n_thin)  = tmp_loglike_s(k_this);
        loglike_c(k / n_thin)  = tmp_loglike_c(k_this);
    end
end
end