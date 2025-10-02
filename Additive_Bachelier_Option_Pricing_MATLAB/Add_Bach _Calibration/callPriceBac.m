function callPrice = callPriceBac(B, ttm, K, F0, sigma)
% callPriceBac - Calculate call option price using the Bachelier model
%
% Inputs:
%   B     - Discount factor for the option maturity
%   ttm   - Time to maturity in years
%   K     - Strike price
%   F0    - Forward price of the underlying asset
%   sigma - Volatility parameter (normal volatility)
%
% Output:
%   callPrice - Call option price under the Bachelier model

 % Calculate normalized moneyness
 x = K - F0;
 y = x / (sigma * sqrt(ttm));
 
 % Calculate the option value using the Bachelier formula
 % This uses the standard normal CDF and PDF functions
 cb = -y .* normcdf(-y) + normpdf(-y);
 
 % Scale by discount factor and volatility term
 callPrice = B * sigma * sqrt(ttm) * cb;
end