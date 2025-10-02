function model_vols = addBachVols(optionPrices, strikes, fwdPrices, disc, ttm, sigmaATM, eta, kappa,I0, a, alpha)
% Calculate model-implied volatilities using characteristic function and FFT
%
% INPUTS:
%   optionPrices - Matrix of market option prices
%   strikes      - Vector of strike prices
%   fwdPrices    - Vector of forward prices for each maturity
%   disc         - Vector of discount factors
%   ttm          - Vector of times to maturity
%   sigmaATM     - Vector of at-the-money volatilities
%   eta          - Model parameter (skewness)
%   kappa        - Model parameter (jump intensity)
%   a            - Damping parameter for FFT
%   alpha        - Model parameter
%
% OUTPUT:
%   model_vols   - Matrix of implied volatilities (maturities x strikes)

% Cumulant generating function
psi_fun = @(u,kappa) (1./kappa) .* (1-alpha)./alpha .* ...
    ( 1 - (1 + (u.*kappa)./(1-alpha)).^alpha );

% Characteristic function
phi_fun = @(u,sigma,eta,kappa,t) ...
    exp( psi_fun( 1i*u.*eta.*sigma.*sqrt(t) + ...
                 0.5*u.^2 .* sigma.^2 .* t , ...
                 kappa) ...
         + 1i*u.*eta.*sigma.*sqrt(t) );

% FFT parameters
M = 15;     % Grid size
dz = 0.01;  % Grid spacing

% Lewis formula for option pricing
C_price_Lewis = @(x,t,eta,kappa,sigma_t,B,Ra) ...
    B.* (Ra + (exp(-a*x)/(2*pi)) .* ...
         FFTBachelier(@(u) phi_fun(u, sigma_t, eta, kappa, t),M,dz,x,a));

% Initialize output
nTTM = length(ttm);
model_vols = NaN(nTTM, length(strikes));

% Main calibration loop
for i = 1:nTTM
    for j = 1:length(strikes)
        % Process only valid option prices
        if ~isnan(optionPrices(2*i-1, j))
            % Normalized moneyness
            mon_norm = (strikes(j) - fwdPrices(i))/(sqrt(ttm(i))*sigmaATM(i));
            
            % Objective function: market price - model price
            objective = @(sigma) callPriceBac(disc(i+1), ttm(i), strikes(j), fwdPrices(i), sigma) - ...
                C_price_Lewis(mon_norm,1,eta,kappa,1/I0,1,0) * disc(i+1)*sigmaATM(i)*sqrt(ttm(i));
            
            % Solve for implied volatility
            model_vols(i, j) = fzero(objective, 20);
        end
    end
end

end