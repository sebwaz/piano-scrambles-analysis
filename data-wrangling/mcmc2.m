close all

% Load the d's
dp = readmatrix("dprimes.csv");
R  = readmatrix("bilinear_fit_R.csv");
F  = readmatrix("bilinear_fit_F.csv");

% last_sample = [R(:, 1)', F(1, :), R(:, 2)', F(2, :), 1, 1];
% [samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(dp, 500000, 1000, last_sample, [], []);

%% Simulation parameters
% N = 80;
% T = 6;
% 
% % Randomly generate Rs (SVD tends to fit better when one eigenvalue smaller than the other)
% R1 = rand(1, N) * 4;
% R2 = rand(1, N) * 0.5;
% 
% % Randomly generate Fs
% F1 = rand(1, T);
% F2 = rand(1, T);
% 
% % Orthogonalize Fs using Gram-Schmidt
% Q = gram_schmidt([F1', F2']);
% F1 = Q(:, 1)';
% F2 = Q(:, 2)';
% 
% % Scale Fs so that sum(F) == length(F)
% F1 = T * F1 / sum(F1);
% F2 = T * F2 / sum(F2);
% 
% % Randomly genereate SDs
% s1 = gamrnd(1, 1) * 0.1;
% s2 = gamrnd(1, 1) * 0.1;
% 
% % Generate the data according to the two-step model
% dp1 = R1' * F1 + randn(length(R1), length(F1)) * s1;
% dp2 = R2' * F2 + randn(length(R2), length(F2)) * s2;
% dp = dp1 + dp2;
% 
% % Organize the true parameters like they will be organized in mcmc
% true_params = [R1, F1, R2, F2, s1, s2];
% 
% % Get the bilinear model fit to the data using SVD
% [Rb, Fb] = fit_bilinear_model(dp);
% bfit = [Rb(:, 1)', Fb(1,:), Rb(:, 2)', Fb(2, :), nan, nan];
% 
% % Compare how well SVD recovers the parameters
% % R2 estimates are in blue (these tend to be consistently under estimated)
% % F2 estimates are in red (they tend to diverge more than other params)
% figure;
% scatter([R1, F1], bfit(1:end-2-N-T), 'ok'); hold on;
% scatter(R2, bfit(N+T+1:N+T+N), 'ob');
% scatter(F2, bfit(end-2-T+1:end-2), 'or');



%% MCMC
% Do the burn-in
[samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(dp, 1000000, 1000, [], [], []);

% Do the final posterior samples
[samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(dp, 1000000, 1000, samples(end, :), sigmas, sigma_scalar);

% Save results
% save("mcmc2_samples_simulated_dataX.mat");
% writematrix(samples, "mcmc2_samples.csv");

% Plot the chains
figure;
plot(samples);

% % Compare the median posterior to the SVD  estimates
% figure;
% plot(median(samples)); hold on;
% plot(true_params); hold on;
% plot(bfit);

% % Same as above, but using scatter plot instead of two functions
% figure;
% scatter(median(samples), true_params, 'ok'); hold on;
% xlim([-30, 30]);
% ylim([-30, 30]);
% plot([-30, 30], [-30, 30]);
% axis square

% Same as above, but using scatter plot instead of two functions
figure;
N = size(R, 1);
T = size(F, 2);
x = median(samples(:, 1:end-2));
y = [R(:, 1)', F(1,:), R(:, 2)', F(2, :)];
subplot(121);
scatter(x(1:N),         y(1:N),         'xb'); hold on;
scatter(x(N+1:N+T),     y(N+1:N+T),     'xr'); hold on;
xlim([-1, 5]);
ylim([-1, 5]);
xlabel("median posterior estimate");
ylabel("SVD estimate");
title("Component 1");
plot([-1, 5], [-1, 5], ':k');
legend({'R_1', 'F_1', 'x = y'});
axis square

subplot(122);
scatter(x(N+T+1:N+T+N), y(N+T+1:N+T+N), 'xb'); hold on;
scatter(x(2*N+T+1:end), y(2*N+T+1:end), 'xr'); hold on;
xlim([-30, 30]);
ylim([-30, 30]);
xlabel("median posterior estimate");
ylabel("SVD estimate");
title("Component 2");
plot([-30, 30], [-30, 30], ':k');
legend({'R_2', 'F_2', 'x = y'});
axis square

% Autocorrelation plots, one by one
figure;
hold off;
for k = 1:size(samples, 2)
    plot(autocorr(samples(:, k)), 'ok');
    ylim([-1, 1]);
    input('');
end

% % Sanity checks
% % (sum(F) == length(F) and F_1 * F_2' == 0)
% sum(samples(:, (size(R, 1) + 1):(size(R, 1) + 6)), 2)
% sum(samples(:, ((size(R, 1) + 1):(size(R, 1) + 6)) + sum(size(dp))), 2)
% for k = 1:size(samples, 1)
%     samples(k, (size(R, 1) + 1):(size(R, 1) + 6)) * samples(k, ((size(R, 1) + 1):(size(R, 1) + 6)) + sum(size(dp)))'
% end



%% Functions
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
% Remember: proposal distribution is allowed to be anything. Don't fuss
% with derivations, as long as constraints on candidates work.
% Cite: https://stats.stackexchange.com/questions/207496/proposal-distribution-metropolis-hastings-mcmc
function tmp = get_candidate(last_sample, mins, maxes, sigmas, dp)
while true
    % Generate a candidate
    tmp = last_sample + sigmas .* randn(size(sigmas));
    
    % Constrain F_1 and F_2 to be orthonormal by putting them into
    % Gram-Schmidt with F_1 first
    idx1 = (size(dp, 1) + 1):(size(dp, 1) + size(dp, 2));
    idx2 = ((size(dp, 1) + 1):(size(dp, 1) + size(dp, 2))) + size(dp, 2) + size(dp, 1);
    Q = gram_schmidt([tmp(idx1)', tmp(idx2)']);
    tmp(idx1) = Q(:, 1);
    tmp(idx2) = Q(:, 2);
    
    % Constrain sum(F) = length(F)
    tmp(idx1) = size(dp, 2) * tmp(idx1) / sum(tmp(idx1));
    tmp(idx2) = size(dp, 2) * tmp(idx2) / sum(tmp(idx2));
    
%     % If F_2 is more similar (cosine) to the mean d', swap F_1 and F_2 candidates
%     if pdist([mean(dp); tmp(idx1)], "cosine") > pdist([mean(dp); tmp(idx2)], "cosine")
%         buf = tmp(idx1);
%         tmp(idx1) = tmp(idx2);
%         tmp(idx2) = buf;
%     end
    
    % Return candidate only if it falls within bounds
    if (all(tmp > mins) && all(tmp < maxes))
        break
    end
end
end


% Get log likelihood of given parameter estimates under bilinear model
function sum_log_p = get_log_likelihood(params, dp)
% d' = R_s_1 * F_t_1 + R_s_2 * F_t_2 + err
R_s_1 = params( 1:size(dp, 1) );
R_s_2 = params((1:size(dp, 1)) + size(dp, 2) + size(dp, 1));
F_t_1 = params( (size(dp, 1) + 1):(size(dp, 1) + size(dp, 2)) );
F_t_2 = params(((size(dp, 1) + 1):(size(dp, 1) + size(dp, 2))) + size(dp, 2) + size(dp, 1));
stdev1 = params(end-1);
stdev2 = params(end);
dp_hat1 = R_s_1' * F_t_1;
dp_hat2 = R_s_2' * F_t_2;

% Log likelihood
log_p1 = log(normpdf(dp,           dp_hat1, stdev1));
log_p2 = log(normpdf(dp - dp_hat1, dp_hat2, stdev2));
log_p = [log_p1; log_p2];
% log_p = log(normpdf(dp, dp_hat1 + dp_hat2, stdev1 + stdev2));
sum_log_p = sum(log_p(:));
end


% Function to estimate the parameters using MCMC
function [samples, candids, loglike_s, loglike_c, sigmas, sigma_scalar] = mcmc(d_prime, n_samples, n_thin, last_sample, last_sigmas, last_sigma_scalar)
% Use a shorter variable name
dp = d_prime;

% Count the number of parameters
% 2* Number subjects (R_s) + 2 * Number tasks (F_t) + one variance parameter
n_params = 2 * size(dp, 1) + 2 * size(dp, 2) + 2;

% Block size for updating the MCMC sigmas
update_block_sz = 5000;

% Prior on all parameters ~ Unif(-100, 100) except sigma ~ Unif(0, 100)
mins  = [-100 * ones(1, n_params - 2), 0, 0];
maxes = 100 * ones(1, n_params);

% Default for last_sample
if isempty(last_sample)
    last_sample = [zeros(1, size(dp, 1)), ...
                   ones(1,  size(dp, 2)), ...
                   zeros(1, size(dp, 1)), ...
                   ones(1,  size(dp, 2)), 1, 1];
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


% Function that performs Gram-Schmidt process to get orthonormal basis
% vectors, stored in Q
function Q = gram_schmidt(A)
    % source: https://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end