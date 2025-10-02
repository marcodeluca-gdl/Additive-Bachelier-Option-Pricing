function [deltaPriceSigma, sigmatnew] = deltaPriceVolatility(sigmaATM,optionPricesNew, strikes, strike, fwdPrices, ttm, ttmlong, disc, alpha, a, eta0, k0, Nsim, increment, priceBac, Notional)
% DELTAPRICESIGMA - Calculates price sensitivity to volatility changes
% This function computes the delta (sensitivity) of an exotic option price
% to changes in each volatility parameter using a bump-and-revalue approach
%
% INPUTS:
%   sigmaATM        - At-the-money volatility vector
%   optionPricesNew - Market option prices for model calibration
%   strikes         - Strike prices for market options
%   strike          - Strike price of the exotic option
%   fwdPrices       - Forward prices
%   ttm             - Time to maturity
%   ttmlong         - Time to maturity
%   disc            - Discount factors
%   alpha, a        - Model parameters
%   eta0, k0        - Initial calibration parameters
%   Nsim            - Number of Monte Carlo simulations
%   increment       - Step size for simulation
%   priceBac        - Base price of exotic option (unshocked)
%   Notional        - Notional amount
%
% OUTPUT:
%   deltaPriceSigma - Vector of price sensitivities to each volatility parameter
%   sigmatnew - Model parameter sigmat after bump

% Define shock size: 1 basis point = 0.0001 = 0.01%
bp = 1e-4;

% Initialize output vector - stores sensitivity for each volatility parameter
deltaPriceSigma = zeros(1, length(sigmaATM));

sigmatnew = zeros(1, length(sigmaATM));

% Loop through each volatility parameter to calculate individual sensitivities
for i=1:length(sigmaATM)
    
    % Create a copy of the original volatility vector
    sigmaATMShocked = sigmaATM;
    
    % Apply shock to the i-th volatility parameter
    % Shock size: 100 * bp = 100 * 0.0001 = 0.01 = 1%
    sigmaATMShocked(i) = sigmaATM(i) + 100*bp;
    
    % Recalibrate the model with the shocked volatility
    % This ensures model consistency after the volatility shock
    [eta1, kappa1, I0_1] = calibrateAddBach(optionPricesNew, sigmaATMShocked, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, 2, length(ttm));
    
    % Adjust volatility term structure using calibration factor
    sigmat1 = sigmaATMShocked/I0_1;
    sigmatnew(i) = sigmat1(i);
    
    % Price the exotic option with shocked parameters
    % Using Bachelier model ('Bac' flag) and recalibrated parameters
    [temp, ~] = pricingExotic(strike, disc, ttmlong, fwdPrices, Nsim, increment, eta1, kappa1, sigmat1, alpha, 0, 0, Notional, 'Bac');
    
    % Calculate the price sensitivity (delta)
    % Difference between shocked price and base price
    deltaPriceSigma(i) = temp - priceBac;
    
end % End of volatility parameter loop

end % End of function