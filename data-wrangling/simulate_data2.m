function nHit_nMiss_nCR_nFA = simulate_data2(R1, F1, R2, F2, nTrials)
%SIMULATE_DATA Simulates data based on bilinear model (and SDT model)
%   R : nsubj x 1 vector
%   F : 1 x ncond vector 
%   nTrials : an even scalar (i.e., nTrials = 2 x # signal trials = 2 x # noise trials)
%   Resulting tensor is length(R) x length(F) x 4


% Initialize the output tensor
nHit_nMiss_nCR_nFA = nan(length(R1), length(F1), 4);

% Calculate the predicted d's
dp = R1 * F1 + R2 * F2;

% Assume no bias, equal variance SDT model
for r = 1:length(R1)
    for c = 1:length(F1)
        
        % Generate proximal stimulus level
        signalTrials = randn(1, nTrials / 2) + dp(r, c);
        noiseTrials  = randn(1, nTrials / 2);
        
        % Optimal criterion under this model is d'/2
        nHit  = sum(signalTrials > dp(r, c) / 2);
        nMiss = sum(signalTrials < dp(r, c) / 2);
        nCR   = sum(noiseTrials  < dp(r, c) / 2);
        nFA   = sum(noiseTrials  > dp(r, c) / 2);
        
        % Save the results
        nHit_nMiss_nCR_nFA(r, c, :) = [nHit, nMiss, nCR, nFA];
    end
end

end

