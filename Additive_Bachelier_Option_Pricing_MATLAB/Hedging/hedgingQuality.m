function [P_L, priceNew] = hedgingQuality(callDir, putDir, target, t0, alpha, a, eta0, k0, Nsim, increment,...
    Notional, callVolatility, futuresInfoVol, bullInfo, numberFutureEta, strangleInfo, nFutureKappa, ptfValue)

% HEDGING QUALITY FUNCTION
% This function calculates the profit and loss (P&L) of a hedged portfolio
% by comparing the portfolio value at a target date with the original portfolio value
% 
% INPUTS:
% - callDir, putDir: directories containing call and put option data
% - target: target date for P&L calculation
% - t0: initial date
% - alpha, a, eta0, k0: model parameters for the additive Bachelier model
% - Nsim: number of Monte Carlo simulations
% - increment: simulation time increment
% - Notional: notional amount of the exotic option
% - callVolatility: matrix containing call option volatility information
% - futuresInfoVol: futures contract information for volatility hedging
% - bullInfo: bull spread option information
% - numberFutureEta: number of futures contracts for eta hedging
% - strangleInfo: strangle option information
% - nFutureKappa: number of futures contracts for kappa hedging
% - ptfValue: original portfolio value
%
% OUTPUTS:
% - P_L: profit and loss of the hedged portfolio
% - priceNew: price of the exotic option using the calibrated model

%% STEP 1: DEFINE CHARACTERISTIC FUNCTION COMPONENTS

% Define the cumulant generating function (psi_fun) for the additive Bachelier model
% This function captures the non-Gaussian features of the underlying asset returns
psi_fun = @(u,kappa) (1./kappa) .* (1-alpha)./alpha .* ...
    ( 1 - (1 + (u.*kappa)./(1-alpha)).^alpha );

% Define the characteristic function (phi_fun) used for option pricing
% This incorporates both the drift (eta) and volatility (sigma) components
phi_fun = @(u,sigma,eta,kappa,t) ...
    exp( psi_fun( 1i*u.*eta.*sigma.*sqrt(t) + ...
                 0.5*u.^2 .* sigma.^2 .* t , ...
                 kappa) ...
         + 1i*u.*eta.*sigma.*sqrt(t) );

%% STEP 2: SET FFT PARAMETERS FOR NUMERICAL INTEGRATION

M = 15;      % Grid size parameter for FFT (determines accuracy)
dz = 0.01;   % Grid spacing for FFT (determines resolution)

%% STEP 3: DEFINE LEWIS FORMULA FOR OPTION PRICING

% Lewis formula implementation using FFT for efficient option pricing
% This allows for fast calculation of option prices under the additive Bachelier model
C_price_Lewis = @(x,t,eta,kappa,sigma_t,B,Ra) ...
    B.* (Ra + (exp(-a*x)/(2*pi)) .* ...
    FFTBachelier(@(u) phi_fun(u, sigma_t, eta, kappa, t),M,dz,x,a));

%% STEP 4: BUILD MARKET DATA AND CALIBRATE MODEL

% Extract option prices and strikes from market data directories
[optionPrices, strikes] = buildOptionPrices(callDir, putDir, target);

% Bootstrap discount factors and forward prices from option market data
[disc, fwdPrices] = bootstrap(optionPrices, strikes);

% Set up time structure using ACT/365 day count convention
ACT_365 = 3;
dates = getDates(t0);
ttm = yearfrac(dates(1), dates(2:end), ACT_365);  % time to maturity for each date

% Apply put-call parity to ensure consistency in option prices
optionPricesNew = put_call_parity(optionPrices(3:end,:), strikes, disc, fwdPrices);

% Calculate at-the-money implied volatilities from market option prices
sigmaATM = ATMvols(optionPricesNew, strikes, fwdPrices, disc, ttm);

% Calibrate the additive Bachelier model parameters (eta, kappa) to market data
% This step fits the model to observed option prices
[eta, kappa, I0] = calibrateAddBach(optionPricesNew, sigmaATM, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, 2, length(ttm));

% Scale the volatility by the calibration factor
sigmat = sigmaATM/I0;

%% STEP 5: PRICE THE EXOTIC OPTION

% Set up time grid including initial time
ttmlong = [0; ttm];
strike = fwdPrices(end);  % Use final forward price as strike

% Price the exotic option using Monte Carlo simulation under the calibrated model
[priceBac, ~] = pricingExotic(strike, disc, ttmlong, fwdPrices, Nsim, increment, eta, kappa, sigmat, alpha, 0, 0, Notional, 'Bac');

%% STEP 6: VALUE HEDGING INSTRUMENTS

% Value call options used for volatility hedging
priceCallVol = zeros(1, length(callVolatility(1,:)));
for i = 1:length(callVolatility(1,:))
    t = callVolatility(4,i);  % Time index for this call option
    % Price call option using Lewis formula with calibrated parameters
    priceCallVol(i) = C_price_Lewis((callVolatility(2,i) - fwdPrices(t)')/(sigmaATM(t)*sqrt(ttm(t))), ...
        1, eta, kappa, 1, 1, 0) * sigmaATM(t) * sqrt(ttm(t)) * disc(t);
end

% Value futures contracts used for volatility hedging
% Futures value = forward price * number of contracts
valueFuturesVol = fwdPrices(futuresInfoVol(3,:))' .* futuresInfoVol(4,:);

% Value bull spread (long lower strike call, short higher strike call)
bullValue = C_price_Lewis((bullInfo(2) - fwdPrices(end))/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, 1, 1, 0) * sigmaATM(end) * sqrt(ttm(end)) * disc(end) - ...
            C_price_Lewis((bullInfo(3) - fwdPrices(end))/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, 1, 1, 0) * sigmaATM(end) * sqrt(ttm(end)) * disc(end);

% Value strangle (long call + long put with different strikes)
% The second term adjusts for the put option using put-call parity
strangleValue = C_price_Lewis((strangleInfo(2) - fwdPrices(end))/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, 1, 1, 0) * sigmaATM(end) * sqrt(ttm(end)) * disc(end) - ...
                (C_price_Lewis((strangleInfo(3) - fwdPrices(end))/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, 1, 1, 0) * (sigmaATM(end)*sqrt(ttm(end))*disc(end)) - disc(end)*(fwdPrices(end) - strangleInfo(3)));

%% STEP 7: CALCULATE TOTAL PROFIT AND LOSS

% Sum all portfolio components to get total P&L
P_L = priceBac + ...                                    % Exotic option value
      sum(priceCallVol .* callVolatility(5,:)) + ...    % Call options weighted by positions
      sum(valueFuturesVol) + ...                        % Futures for volatility hedging
      bullValue * bullInfo(4) + ...                     % Bull spread weighted by position
      numberFutureEta * fwdPrices(end) + ...           % Futures for eta hedging
      strangleValue + ...                               % Strangle option value
      nFutureKappa * fwdPrices(end) - ...              % Futures for kappa hedging
      ptfValue;                                         % Subtract original portfolio value

% Return the new price of the exotic option
priceNew = priceBac;

end