function sigmaATM = ATMvols(optionPrices, strikes, fwdPrices, disc, ttm)
% ATMvols - Calculate At-The-Money implied volatilities from option prices
%
% Inputs:
%   optionPrices - Matrix of option prices, odd rows for calls, even rows for puts
%   strikes      - Vector of strike prices
%   fwdPrices    - Vector of forward prices corresponding to each maturity
%   disc         - Vector of discount factors (disc(i+1) corresponds to ttm(i))
%   ttm          - Vector of time-to-maturity values in years
%
% Outputs:
%   sigmaATM     - Vector of ATM implied volatilities for each maturity

% Get number of maturities
nTTM = length(ttm);

% Initialize output array
sigmaATM = zeros(1, nTTM);

% Loop through each maturity
for i = 1:nTTM
    % Interpolate call price at the ATM forward level
    callMkt = interp1(strikes, optionPrices(2*i-1,:), fwdPrices(i), 'spline');
    
    % Define objective function for root-finding:
    % Find volatility that makes model price match market price
    minimizer = @(sigma) callPriceBac(disc(i+1), ttm(i), fwdPrices(i), fwdPrices(i), sigma) - callMkt;
    
    % Solve for ATM volatility using fzero with initial guess of 20 
    sigmaATM(i) = fzero(minimizer, 20);
end