function [P, Rsq, Fstat, ssr, ssx, ssres, dfx, dfres] = do_tests(R, F, dprimes)


% Total sum of squares (unexplained variance)
sst = sum((dprimes(:) - mean(dprimes(:))).^2);

%% F-test for each bilinear model relative to previous model
dfx   = []; %'Extra' df, i.e., df current model minus df previous model
dfres = []; % Total df minus sum of df's of previous models
ssr   = []; % Sum of squares explained by current model
ssx   = []; %'Extra' sum of squares explained by current model relative to previous model
ssres = [];
msx   = [];
msres = [];
Fstat = [];
P = [];

for k = 1:size(dprimes,2)
    
    % Degrees of freedom calculations
    dfx(k)   = size(R(:, k), 1) + size(F(k, :), 2) - (2 * k - 1); % The extra is likely due to the constraint that each component must be orthogonal to previous components
    dfres(k) = size(R(:, k), 1) * size(F(k, :), 2) - sum(dfx);
    
    % SS calculations
    err = R(:, 1:k) * F(1:k, :) - dprimes;
    ssres(k) = sum(err(:).^2);
    ssr(k)   = sst - ssres(k);
    if k > 1
        ssx(k) = ssr(k) - ssr(k-1);
    else
        ssx(k) = ssr(k);
    end
    
    % MSE calculations
    msx(k)   = ssx(k) / dfx(k);
    msres(k) = ssres(k) / dfres(k);
    
    % Fstat
    Fstat(k) = msx(k) / msres(k);
    
    % p-values
    P(k) = 1 - fcdf(Fstat(k), dfx(k), dfres(k));
end

%Rsq of each model
Rsq = ssr / sst;

end