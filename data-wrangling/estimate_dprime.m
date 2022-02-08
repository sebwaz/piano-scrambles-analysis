function [dprimes, bias, accuracy] = estimate_dprime(nHit_nMiss_nCR_nFA)
%ESTIMATE_DPRIME Calculates d's from a tensor of counts of hits, misses and
%correct rejections.


% Initialize output matrices
inputSize = size(nHit_nMiss_nCR_nFA);
dprimes  = nan(inputSize(1), inputSize(2));
bias     = nan(inputSize(1), inputSize(2));
accuracy = nan(inputSize(1), inputSize(2));

for r = 1:size(nHit_nMiss_nCR_nFA, 1)
    for c = 1:size(nHit_nMiss_nCR_nFA, 2)
        
        % Calculate # signal and noise trials
        nSignal = sum(nHit_nMiss_nCR_nFA(r, c, 1:2));
        nNoise  = sum(nHit_nMiss_nCR_nFA(r, c, 3:4));
        
        % Calculate p(hit) and p(false alarm)
        pHit = nHit_nMiss_nCR_nFA(r, c, 1) / nSignal;
        pFA  = nHit_nMiss_nCR_nFA(r, c, 4) / nNoise;
        
        % Adjustments based on Macmillan & Creelman, 2nd edition (2005), page 8
        if pHit == 1
            pHit = (nSignal - 0.5) / nSignal;
        end
        
        if pFA == 1
            pFA = (nNoise - 0.5) / nNoise;
        end
        
        if pHit == 0
            pHit = 0.5 / nSignal;
        end
        
        if pFA == 0
            pFA = 0.5 / nNoise;
        end
        
        % Calculate d-prime
        dprimes(r, c) = norminv(pHit) - norminv(pFA);
        
        % Calculate bias
        bias(r, c) = norminv(1 - pFA) - dprimes(r, c) / 2;
        
        % Calculate the overall accuracy
        accuracy(r, c) = (nHit_nMiss_nCR_nFA(r, c, 1) + nHit_nMiss_nCR_nFA(r, c, 3)) / (nSignal + nNoise);
    end
end

end

