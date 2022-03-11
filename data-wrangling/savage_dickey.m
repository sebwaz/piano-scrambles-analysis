% Load the d's (just so we know the shape of things)
dp = readmatrix("dprimes.csv");
[N, T] = size(dp);

% Load the samples from mcmc2
samples = readmatrix("mcmc2_samples.csv");

% Get the R2 posterior samples specifically
R2 = samples(:, (N + T + 1):(N + T + N));

% Check for covariance among the R2 posterior samples
% (This is important because we're looking at the joint distribution of R,
% we're but applying a smoother on the marginals)
r = tril(corrcoef(R2), -1);
histogram(r(r ~= 0));

%% Use kernel density estimator to approximate the posterior density of R2

% Go through each R2_k and plot its marginal posterior density estimate
figure;
for k = 1:N
    [f, xi] = ksdensity(R2(:, k));
    plot(xi, f); hold on;
    xlim([-.2, .2]);
end

% Matlab has a multivariate kernel density estimator!
% https://www.mathworks.com/help/stats/mvksdensity.html#bu62r12-bw
% Bandwidth must be specified by hand, but the doc gives a rule of thumb
[n, d] = size(R2);
b = std(R2) * (4 / ((d + 2) * n)) ^ (1 / (d + 4));
[f, xi] = mvksdensity(R2, zeros(1, N), 'Bandwidth', b);

% According to Savage-Dickey method, posterior density at null values should
% be divided by the prior density at the same point:
bf_01 = f / (1 / 200)^d

% Plot each parameter's posterior density, conditioned on the zero value of the others, one at a time
figure;
for k = 1:d
    xslice = linspace(-.2, .2, 1000);
    [f, xi] = mvksdensity(R2, [zeros(1000, k-1), xslice', zeros(1000, N-k)], 'Bandwidth', b);
    plot(xslice, f);
    input('');
end

