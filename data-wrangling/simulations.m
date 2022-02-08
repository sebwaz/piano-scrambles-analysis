
%% Model 1
% Tuning parameters of the simulation
R = [zeros(43, 1); linspace(0, 5, 43)'];
F = ones(1, 6);
nTrialsPerCond = 60;
nIters = 10000;

% To save results
model1_R     = nan(nIters, length(R));
model1_F     = nan(nIters, length(F));
model1_rsq   = nan(nIters, length(F));
model1_fstat = nan(nIters, length(F));
model1_pval  = nan(nIters, length(F));

% Iteration loop
for i = 1:nIters
    
    % Announce progress
    if mod(i, 100) == 0
        clc;
        fprintf('Iteration %d of %d', i, nIters)
    end
    
    % Simulate data from experiment
    nHit_nMiss_nCR_nFA = simulate_data(R, F, nTrialsPerCond);
    
    % Re-estimate d's based on simulated data
    [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    model1_R(i, :)     = R_est(:, 1);
    model1_F(i, :)     = F_est(1, :);
    model1_rsq(i, :)   = Rsq;
    model1_fstat(i, :) = Fstat;
    model1_pval(i, :)  = P;
end


%% Model 2
% Tuning parameters of the simulation
R2 = ones(86, 1);
F2 = [0, 0, 0, 0, 1, 1];

% To save results
model2_R     = nan(nIters, length(R));
model2_F     = nan(nIters, length(F));
model2_rsq   = nan(nIters, length(F));
model2_fstat = nan(nIters, length(F));
model2_pval  = nan(nIters, length(F));

% Iteration loop
for i = 1:nIters
    
    % Announce progress
    if mod(i, 100) == 0
        clc;
        fprintf('Iteration %d of %d', i, nIters)
    end
    
    % Simulate data from experiment
    nHit_nMiss_nCR_nFA = simulate_data2(R, F, R2, F2, nTrialsPerCond);
    
    % Re-estimate d's based on simulated data
    [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA);
    
    % Fit the bilinear model to the data 
    [R_est, F_est] = fit_bilinear_model(dprimes);
    
    % Do fTests, also get goodness-of-fit measures
    [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R_est, F_est, dprimes);
    
    % Save results
    model2_R(i, :)     = R_est(:, 1);
    model2_F(i, :)     = F_est(1, :);
    model2_rsq(i, :)   = Rsq;
    model2_fstat(i, :) = Fstat;
    model2_pval(i, :)  = P;
end


%% Save simulation results
save('simulations.mat');


%% Plot the results
boxplot([model1_rsq(:, 1), model2_rsq(:, 1)]);
title('Goodness-of-fit of the bilinear model under different true models');
xticklabels({'1D model', '2D model'});
xlabel("True model");
ylabel("R-squared");
ylim([0.85, 1]);