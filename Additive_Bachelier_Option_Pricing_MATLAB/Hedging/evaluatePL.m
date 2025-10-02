function [PLNoHedging, PLHedging, PLHedgingCost] = evaluatePL(callDir, putDir, targets, datesPL, alpha, ...
    a, eta0, k0, Nsim, increment, Notional, callVolatility,...
    futuresInfoVol, bullInfo, numberFutureEta, strangleInfo,...
    numberFutureKappa, ptfValue, CostHedging, priceBac)
% Evaluate profit and loss (P&L) for portfolio with and without hedging
%
% INPUTS:
% callDir - Directory path for call option data
% putDir - Directory path for put option data  
% targets - Vector of target dates for P&L evaluation
% datesPL - Vector of corresponding P&L evaluation dates
% alpha - Model parameter (volatility scaling)
% a - Model parameter (damping factor)
% eta0 - Initial eta parameter (skewness)
% k0 - Initial kappa parameter (jump intensity)
% Nsim - Number of Monte Carlo simulations
% increment - Time increment for simulation
% Notional - Portfolio notional amount
% callVolatility - Call option volatility data
% futuresInfoVol - Futures information for volatility hedging
% bullInfo - Bull spread information for eta hedging
% numberFutureEta - Number of futures contracts for eta hedging
% strangleInfo - Strangle information for kappa hedging
% numberFutureKappa - Number of futures contracts for kappa hedging
% ptfValue - Initial portfolio value
% CostHedging - Total cost of hedging strategy
% PriceBac - Exotic price
%
% OUTPUTS:
% PLNoHedging - P&L without hedging (portfolio value change only)
% PLHedging - P&L with hedging (including hedge performance)
% PLHedgingCost - P&L with hedging after deducting hedging costs

% Initialize output vectors
PLNoHedging = zeros(1, length(targets));
PLHedging = zeros(1, length(targets));
priceBacNew = zeros(1, length(targets));

% Loop through each target evaluation date
for i = 1:length(targets)
    % Calculate hedged P&L and new portfolio price for current target
    [PLHedging(i), priceBacNew(i)] = hedgingQuality(callDir, putDir, targets(i), datesPL(i), alpha, a, eta0, k0, Nsim, increment,...
        Notional, callVolatility, futuresInfoVol, bullInfo, numberFutureEta, strangleInfo, numberFutureKappa, ptfValue);
    
    % Calculate unhedged P&L as difference between new and initial portfolio price
    PLNoHedging(i) = priceBacNew(i) - priceBac;
end

% Calculate net P&L after hedging costs
PLHedgingCost = PLHedging - CostHedging;

end