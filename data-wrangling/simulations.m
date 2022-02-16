
clc

%% 1D data
% Tuning parameters of the simulation
R = [zeros(43, 1); linspace(0, 5, 43)'];
F = ones(1, 6);
nTrialsPerCond = 60;
nIters = 10000;

% To save results
data1_R     = nan(length(R), length(F), nIters);
data1_F     = nan(length(F), length(F), nIters);
data1_rsq   = nan(nIters, length(F));
data1_fstat = nan(nIters, length(F));
data1_pval  = nan(nIters, length(F));

data1_ssr   = nan(nIters, length(F));
data1_ssx   = nan(nIters, length(F));
data1_ssres = nan(nIters, length(F));
data1_dfx   = nan(nIters, length(F));
data1_dfres = nan(nIters, length(F));

% Iteration loop
disp('Simulating 1D data');
parfor i = 1:nIters
    
    % Simulate data from experiment
    nHit_nMiss_nCR_nFA = simulate_data(R, F, nTrialsPerCond);
    
    % Re-estimate d's based on simulated data
    [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    data1_R(:, :, i)  = R_est;
    data1_F(:, :, i)  = F_est;
    data1_rsq(i, :)   = Rsq;
    data1_fstat(i, :) = Fstat;
    data1_pval(i, :)  = P;
    
    data1_ssr(i, :)   = ssr;
    data1_ssx(i, :)   = ssx;
    data1_ssres(i, :) = ssres;
    data1_dfx(i, :)   = dfx;
    data1_dfres(i, :) = dfres;
end


%% 2D data
% Tuning parameters of the simulation
R2 = 5 * ones(86, 1);
F2 = [0, 0, 0, 0, 1, 1];

% To save results
data2_R     = nan(length(R), length(F), nIters);
data2_F     = nan(length(F), length(F), nIters);
data2_rsq   = nan(nIters, length(F));
data2_fstat = nan(nIters, length(F));
data2_pval  = nan(nIters, length(F));

data2_ssr   = nan(nIters, length(F));
data2_ssx   = nan(nIters, length(F));
data2_ssres = nan(nIters, length(F));
data2_dfx   = nan(nIters, length(F));
data2_dfres = nan(nIters, length(F));

% Iteration loop
disp('Simulating 2D data');
parfor i = 1:nIters
    
    % Simulate data from experiment
    nHit_nMiss_nCR_nFA = simulate_data2(R, F, R2, F2, nTrialsPerCond);
    
    % Re-estimate d's based on simulated data
    [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    data2_R(:, :, i)  = R_est;
    data2_F(:, :, i)  = F_est;
    data2_rsq(i, :)   = Rsq;
    data2_fstat(i, :) = Fstat;
    data2_pval(i, :)  = P;
    
    data2_ssr(i, :)   = ssr;
    data2_ssx(i, :)   = ssx;
    data2_ssres(i, :) = ssres;
    data2_dfx(i, :)   = dfx;
    data2_dfres(i, :) = dfres;
end


%% Simulate data without SDT simulation, just use a fixed random normal
sigma = 3;

% To save results
data3_R     = nan(length(R), length(F), nIters);
data3_F     = nan(length(F), length(F), nIters);
data3_rsq   = nan(nIters, length(F));
data3_fstat = nan(nIters, length(F));
data3_pval  = nan(nIters, length(F));

data3_ssr   = nan(nIters, length(F));
data3_ssx   = nan(nIters, length(F));
data3_ssres = nan(nIters, length(F));
data3_dfx   = nan(nIters, length(F));
data3_dfres = nan(nIters, length(F));

disp('Simulating 1D data, just Gaussian');
parfor i = 1:nIters
    % Generate d's directly from model, no SDT stimulation
    dprimes = normrnd(R * F, sigma);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    data3_R(:, :, i)  = R_est;
    data3_F(:, :, i)  = F_est;
    data3_rsq(i, :)   = Rsq;
    data3_fstat(i, :) = Fstat;
    data3_pval(i, :)  = P;
    
    data3_ssr(i, :)   = ssr;
    data3_ssx(i, :)   = ssx;
    data3_ssres(i, :) = ssres;
    data3_dfx(i, :)   = dfx;
    data3_dfres(i, :) = dfres;
end


%% Save simulation results
% save('simulations.mat');
% load('simulations.mat');


%% Plot the results
close all
figure;
boxplot([data1_rsq(:, 1), data2_rsq(:, 1)]);
title('Goodness-of-fit of the bilinear model under different true models');
xticklabels({'1D model', '2D model'});
xlabel("True model");
ylabel("R-squared");
ylim([0.6, 1]);


%% Look at observed vs. predicted d' plots
% 1D data
nHit_nMiss_nCR_nFA = simulate_data(R, F, nTrialsPerCond);
[dprimes, ~, ~] = estimate_dprime(nHit_nMiss_nCR_nFA);
[R_est, F_est] = fit_bilinear_model(dprimes);

figure;
x = dprimes;
y = (R_est(:, 1) * F_est(1, :));
c = repmat(1:6, [length(R), 1]);
for k = 1:6
    scatter(x(:, k), y(:, k), [], c(:, k), 'filled'); hold on;
end
xlim([-2, 6]);
ylim([-2, 6]);
xlabel("Observed d' (2D true model)");
ylabel("Predicted d' (bilinear model)");
line([-2, 6], [-2, 6], 'LineStyle', '--');
legend({'1', '2', '3', '4', '5', '6'});
title('Example fit to 1D data');


% 2D data
nHit_nMiss_nCR_nFA = simulate_data2(R, F, R2, F2, nTrialsPerCond);
[dprimes, ~, ~] = estimate_dprime(nHit_nMiss_nCR_nFA);
[R_est, F_est] = fit_bilinear_model(dprimes);

figure;
x = dprimes;
y = (R_est(:, 1) * F_est(1, :));
c = repmat(1:6, [length(R), 1]);
for k = 1:6
    scatter(x(:, k), y(:, k), [], c(:, k), 'filled'); hold on;
end
xlim([-2, 6]);
ylim([-2, 6]);
xlabel("Observed d' (2D true model)");
ylabel("Predicted d' (bilinear model)");
line([-2, 6], [-2, 6], 'LineStyle', '--');
legend({'1', '2', '3', '4', '5', '6'});
title('Example fit to 2D data');


% 1D data, no SDT, just Gaussian
dprimes = normrnd(R * F, sigma);
[R_est, F_est] = fit_bilinear_model(dprimes);

figure;
x = dprimes;
y = (R_est(:, 1) * F_est(1, :));
c = repmat(1:6, [length(R), 1]);
for k = 1:6
    scatter(x(:, k), y(:, k), [], c(:, k), 'filled'); hold on;
end
xlim([-2, 6]);
ylim([-2, 6]);
xlabel("Observed d' (2D true model)");
ylabel("Predicted d' (bilinear model)");
line([-2, 6], [-2, 6], 'LineStyle', '--');
legend({'1', '2', '3', '4', '5', '6'});
title('Example fit to 1D data, no SDT, just Gaussian');


% Histograms of 1D F estimates
figure;
for k = 1:6
    subplot(1, 6, k);
    histogram(data1_F(1, k, :)); hold on;
    line([1, 1], [0, 800], 'LineStyle', '--');
    title(sprintf('Cond %d, mean = %5.3f', k, mean(data1_F(1, k, :))));
%     xlim([0.85, 1.25]);
%     ylim([0, 800]);
end
sgtitle('Simulated F estimates (1D data)');


% Histograms of 1D F2 estimates
figure;
for k = 1:6
    subplot(1, 6, k);
    histogram(data1_F(2, k, :)); hold on;
    line([1, 1], [0, 800], 'LineStyle', '--');
    title(sprintf('Cond %d, mean = %5.3f', k, mean(data1_F(2, k, :))));
%     xlim([0.85, 1.25]);
%     ylim([0, 800]);
end
sgtitle('Simulated F2 estimates (1D data)');


% Histogram of 2D F estimates
figure;
for k = 1:6
    subplot(1, 6, k);
    histogram(data2_F(1, k, :));
    line([1, 1], [0, 800], 'LineStyle', '--');
    title(sprintf('Cond %d, mean = %5.3f', k, mean(data2_F(1, k, :))));
%     xlim([0.85, 1.25]);
%     ylim([0, 800]);
end
sgtitle('Simulated F estimates (2D data)');


% Histogram of 2D F2 estimates
figure;
for k = 1:6
    subplot(1, 6, k);
    histogram(data2_F(2, k, :));
    line([1, 1], [0, 800], 'LineStyle', '--');
    title(sprintf('Cond %d, mean = %5.3f', k, mean(data2_F(2, k, :))));
%     xlim([0.85, 1.25]);
%     ylim([0, 800]);
end
sgtitle('Simulated F2 estimates (2D data)');


%% F-test
% Histogram of p-values
figure;
subplot(311);
H = histogram(data1_pval(:, 2));
title('1D data');
xlim([0, 1]);
subplot(312);
histogram(data2_pval(:, 2), 'BinEdges', H.BinEdges);
title('2D data');
xlim([0, 1]);
subplot(313);
histogram(data3_pval(:, 2), 'BinEdges', H.BinEdges);
title('1D data, no SDT, just Gaussian');
xlim([0, 1]);

% The above doesn't look right. Check that numerator and denominator for 1D
% data are distributed according to chi2 (within a scaling factor)

% Numerator distribution
ssx = data1_ssx(:, 2);
[nx, edgesx] = histcounts(ssx, 'Normalization', 'probability');

% Approximate mode of numerator distribution
midsx = (edgesx(1:end - 1) + edgesx(2:end)) / 2;
[~, argmaxx] = max(nx);
ssx_mode = midsx(argmaxx);

% Denominator distribution
ssres = data1_ssres(:, 2);
[nres, edgesres] = histcounts(ssres, 'Normalization', 'probability');

% Approximate mode of denominator distribution
midsres = (edgesres(1:end - 1) + edgesres(2:end)) / 2;
[~, argmaxres] = max(nres);
ssres_mode = midsres(argmaxres);

% The intended modes (according to chi2 distribution) should be dfx-2 and
% dfres-2, respectively. Scale the msres to match, and scale msr by the
% same factor.
null_ratio = (data1_dfx(1, 2) - 2) / (data1_dfres(1, 2) - 2);
obs_ratio  = ssx_mode / ssres_mode;

% Plot the num and den distributions
figure;
subplot(211);
histogram(ssx);
xlim([0, max([ssx; ssres])]);
title('SS_X');
subplot(212);
histogram(ssres);
xlim([0, max([ssx; ssres])]);
title('SS_r_e_s');



%% Simulations based on 2022 Feb 16 meeting
R_fit     = readmatrix("bilinear_fit_R.csv");
F_fit     = readmatrix("bilinear_fit_F.csv");
SSres_fit = readmatrix("bilinear_fit_SSres.csv");
R_fit     = R_fit(:, 1);
F_fit     = F_fit(1, :);
SSres_fit = SSres_fit(1);

nTrialsPerCond = 60;
nIters = 10000;

% To save results
sim_R     = nan(length(R_fit), length(F_fit), nIters);
sim_F     = nan(length(F_fit), length(F_fit), nIters);
sim_rsq   = nan(nIters, length(F_fit));
sim_fstat = nan(nIters, length(F_fit));
sim_pval  = nan(nIters, length(F_fit));

sim_ssr   = nan(nIters, length(F_fit));
sim_ssx   = nan(nIters, length(F_fit));
sim_ssres = nan(nIters, length(F_fit));
sim_dfx   = nan(nIters, length(F_fit));
sim_dfres = nan(nIters, length(F_fit));

% Iteration loop
disp('Simulating data from bilinear fit');
parfor i = 1:nIters
    
    % Simulate data from experiment
    nHit_nMiss_nCR_nFA = simulate_data(R_fit, F_fit, nTrialsPerCond);
    
    % Re-estimate d's based on simulated data
    [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    sim_R(:, :, i)  = R_est;
    sim_F(:, :, i)  = F_est;
    sim_rsq(i, :)   = Rsq;
    sim_fstat(i, :) = Fstat;
    sim_pval(i, :)  = P;
    
    sim_ssr(i, :)   = ssr;
    sim_ssx(i, :)   = ssx;
    sim_ssres(i, :) = ssres;
    sim_dfx(i, :)   = dfx;
    sim_dfres(i, :) = dfres;
end

% Compare the observed SSres to the simulated distribution of SSres
figure;
histogram(sim_ssres(:, 1)); hold on;
line([SSres_fit, SSres_fit], [0, 500], 'LineStyle', ':', 'Color', 'k');
title('Histogram of SS_R_e_s');
xlabel('SS_R_e_s');

% Wow, we're way off! Verify this calculation
dp = readmatrix("dprimes.csv");
err = (R_fit * F_fit) - dp;
SSres_ver1 = sum(err(:).^2);

[~, ~, ~, ~, ~, SSres_ver2, ~, ~] = do_tests(readmatrix("bilinear_fit_R.csv"), readmatrix("bilinear_fit_F.csv"), dp);

% Look at the distribution of F estimates when you generate the data this
% way
figure;
for k = 1:6
    subplot(1, 6, k);
    histogram(sim_F(1, k, :)); hold on;
    line([F_fit(k), F_fit(k)], [0, 800], 'LineStyle', '--');
%     title(sprintf('Cond %d, mean = %5.3f', k, mean(sim_F(1, k, :))));
end
sgtitle('Simulated F estimates (data generated from model fit to true data)');