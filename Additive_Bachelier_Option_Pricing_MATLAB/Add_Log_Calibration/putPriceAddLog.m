function putPrice = putPriceAddLog(s, S0, K)
% Put option pricing using additive logistic model
%
% Inputs:
%   s  - model parameter
%   S0 - Current spot price
%   K  - Strike price
%
% Output:
%   putPrice - Pall option price

    putPrice = s*log(1+exp((K-S0)/s));
    
end