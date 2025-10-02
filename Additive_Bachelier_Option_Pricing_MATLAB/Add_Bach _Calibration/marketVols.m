function implied_vols = marketVols(optionPrices, strikes, fwdPrices, disc, ttm)
% marketVols - Calculate implied volatilities across all strikes and maturities
%
% Inputs:
%   optionPrices - Matrix of option prices, odd rows for calls, even rows for puts
%   strikes      - Vector of strike prices
%   fwdPrices    - Vector of forward prices corresponding to each maturity
%   disc         - Vector of discount factors (disc(i+1) corresponds to ttm(i))
%   ttm          - Vector of time-to-maturity values in years
%
% Outputs:
%   implied_vols - Matrix of implied volatilities (rows=maturities, columns=strikes)

% Get dimensions
nTTM = length(ttm);
nStrikes = length(strikes);

% Initialize output matrix with NaN values
implied_vols = NaN(nTTM, nStrikes); % vol matrix: row=ttm, column=strike

% Double loop through all maturities and strikes
for i = 1:nTTM
    % Loop through all strikes to build the volatility smile for this maturity
    for j = 1:nStrikes
        % Only process if we have a valid call option price
        if ~isnan(optionPrices(2*i-1, j))
            % Get the call price for this maturity and strike
            priceCall = optionPrices(2*i-1, j);
            
            % Define objective function for root-finding
            minimizer2 = @(sigma) callPriceBac(disc(i+1), ttm(i), strikes(j), fwdPrices(i), sigma) - priceCall;
            
            % Solve for implied volatility using fzero with initial guess of 0.2 (20%)
            implied_vols(i, j) = fzero(minimizer2, 0.2);
        end
    end
end