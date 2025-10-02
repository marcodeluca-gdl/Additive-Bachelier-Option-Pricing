function callPrice = callPriceAddLog(s, S0, K)
% Call option pricing using additive logistic model
%
% Inputs:
%   s  - model parameter
%   S0 - Current spot price
%   K  - Strike price
%
% Output:
%   callPrice - Call option price

    callPrice = s*log(1+exp((S0-K)/s));
    
end