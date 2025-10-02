function callInfo = callVolatilityBump(idxTop, fwdPrices, strikes, ...
    optionPricesNew, sigmat1, ttm, disc, ...
    deltaSigmaTop)
% CALLVOLATILITYBUMP - Constructs hedging call options for volatility exposure
% This function calculates the required call options to hedge volatility risk
% by determining optimal strikes, notionals, and maturities for hedging instruments
%
% INPUTS:
%   idxTop           - Indices of top volatility risk contributors (from Pareto analysis)
%   fwdPrices        - Forward prices vector
%   strikes          - Available strike prices for hedging options
%   optionPricesNew  - Market option prices matrix
%   sigmat1          - Volatility term structure
%   ttm              - Time to maturity vector
%   disc             - Discount factors
%   deltaSigmaTop    - Volatility sensitivities to be hedged
%
% OUTPUT:
%   callInfo (5 × nCall matrix):
%   Row 1 → Call option notional amounts
%   Row 2 → Associated strike prices (closest to forward price)
%   Row 3 → Time-to-maturity of each call option
%   Row 4 → Time-to-maturity indexes of each call option
%   Row 5 → Number of calls

% Initialize based on number of hedging instruments needed
nCall = numel(idxTop);
bp=1e-4;

% Pre-allocate output vectors
notionalCall = zeros(1, nCall);  % Row 1: notional amounts
strikeCall = zeros(1, nCall);    % Row 2: strike prices
ttmCall = zeros(1, nCall);       % Row 3: time to maturity
ttmIdx = zeros(1, nCall);        % Row 4: ttm indexes
numCall = zeros(1,nCall);        % Row 5: number of calls

% Loop through each volatility risk factor to construct hedging calls
for i = 1:nCall
    
    % ---------------------- STRIKE SELECTION ----------------------
    % Get forward price for current risk factor
    fwdPrice = fwdPrices(idxTop(i));
    
    % Calculate distance from each available strike to forward price
    % Objective: find at-the-money (ATM) strike for maximum vega exposure
    dist = abs(strikes - fwdPrice);
    
    % Exclude strikes with missing market prices (NaN values)
    % Set distance to infinity for unavailable options
    dist(isnan(optionPricesNew(2*idxTop(i)-1, :))) = Inf;
    
    % Select strike closest to forward price (most liquid ATM option)
    [~, idx] = min(dist);
    closestStrike = strikes(idx);
    
    % ----------------------------- VEGA CALCULATION -------------------------------
    % Calculate vega (sensitivity to volatility) of the selected call option
    % VegaOption returns vega per 1% volatility change
    % Multiply by 100*bp to match the shock size used in sensitivity calculation
    vegaCall = VegaOption(closestStrike, fwdPrice, ...
                         sigmat1(idxTop(i)), ttm(idxTop(i)), ...
                         disc(idxTop(i)+1)) * 100 * bp;
    
    % --------------------------- HEDGING CALCULATION ----------------------------
    % Calculate required notional to hedge the volatility exposure
    % Negative sign: if exotic has positive vega exposure, sell calls (negative notional)
    % Formula: Required Notional = -Risk_to_Hedge / Vega_of_Hedging_Instrument
    notionalCall(i) = -deltaSigmaTop(i) / vegaCall;
    
    % Store strike and maturity information for the hedging call
    strikeCall(i) = closestStrike;
    ttmCall(i) = ttm(idxTop(i));
    ttmIdx(i) = idxTop(i);
    numCall(i) = notionalCall(i) / optionPricesNew(ttmIdx(i)*2-1, idx);

    
end

% Construct output matrix: 3 rows × nCall columns
% Row 1: Notional amounts for each hedging call
% Row 2: Strike prices for each hedging call  
% Row 3: Time to maturity for each hedging call
% Row 4: ttm indexes
% Row 5: Number of each call
callInfo = [notionalCall;
            strikeCall;
            ttmCall;
            ttmIdx;
            ceil(numCall)];

end
