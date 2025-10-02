function model_vols = addLogVols(optionPrices, strikes, fwdPrices, disc, ttm, s)
% Calculate model-implied volatilities using additive logistic model
%
% INPUTS:
% optionPrices - Matrix of market option prices
% strikes - Vector of strike prices
% fwdPrices - Vector of forward prices for each maturity
% disc - Vector of discount factors
% ttm - Vector of times to maturity
% s - Vector of additive logistic model parameters
%
% OUTPUT:
% model_vols - Matrix of implied volatilities from additive logistic model (maturities x strikes)

% Initialize output
nTTM = length(ttm);
model_vols = NaN(nTTM, length(strikes));

% Main calibration loop
for i = 1:nTTM
    for j = 1:length(strikes)
        % Process only valid option prices
        if ~isnan(optionPrices(2*i-1, j))
            % Objective function: Bachelier price - additive logistic model price
            objective = @(sigma) callPriceBac(disc(i+1), ttm(i), strikes(j), fwdPrices(i), sigma) - ...
                                 callPriceAddLog(s(i), fwdPrices(i), strikes(j));
            
            % Solve for implied volatility that matches additive logistic model price
            model_vols(i, j) = fzero(objective, 20);
        end
    end
end
end