function [R, F] = fit_bilinear_model(dprimes)

% The number of conditions
nCond = size(dprimes, 2);

% Do svd (bilinear model is first principal component and first vector of loadings)
[C, E, L] = svd(dprimes, 'econ');
Lt = L';

% Negate C and L if Lt(1,:) not positive
if all(Lt(1,:)<0)
    C  = -C;
    L  = -L;
    Lt = -Lt;
end

% Get R and F
R = [];
F = [];
for k = 1:nCond
    % Get R and F estimates from SVD output
    R(:,k) = C(:,k) * E(k,k);
    F(k,:) = Lt(k,:);

    % Scale so that sum of Fs = nCond
    x = nCond / sum(F(k,:));
    R(:,k) = R(:,k) / x;
    F(k,:) = x * F(k,:);
end

end