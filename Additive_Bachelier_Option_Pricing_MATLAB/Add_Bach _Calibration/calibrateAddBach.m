function [eta, kappa, I0, MSE] = calibrateAddBach(optionPrices, sigmaATM, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, start, last)
% calibrateAddBach - Calibrate parameters for an additive Bachelier model with smile
%
% Inputs:
%   optionPrices - Matrix of option prices, odd rows for calls, even rows for puts
%   sigmaATM     - Vector of ATM implied volatilities for each maturity
%   strikes      - Vector of strike prices
%   fwdPrices    - Vector of forward prices corresponding to each maturity
%   ttm          - Vector of time-to-maturity values in years
%   disc         - Vector of discount factors
%   alpha        - Model parameter
%   a            - Parameter in Lewis Formula
%   eta0, k0     - Initial guesses
%   start, last  - Time interval of the calibration
%
% Outputs:
%   eta          - Fitted volatility-of-volatility parameter
%   kappa        - Fitted mean-reversion or skew parameter

% Define characteristic function components
% psi_fun represents the cumulant generating function
psi_fun = @(u,kappa) (1./kappa) .* (1-alpha)./alpha .* ...
    ( 1 - (1 + (u.*kappa)./(1-alpha)).^alpha );
% phi_fun is the characteristic function
phi_fun = @(u,sigma,eta,kappa,t) ...
    exp( psi_fun( 1i*u.*eta.*sigma.*sqrt(t) + ...
        0.5*u.^2 .* sigma.^2 .* t , ...
        kappa) ...
        + 1i*u.*eta.*sigma.*sqrt(t) );

% FFT parameters
M = 15;    % Grid size parameter
dz = 0.01; % Grid spacing

% Lewis formula for option pricing using FFT
C_price_Lewis = @(x,t,eta,kappa,sigma_t,B,Ra) ...
    B.* (Ra + (exp(-a*x)/(2*pi)) .* ...
    FFTBachelier(@(u) phi_fun(u, sigma_t, eta, kappa, t),M,dz,x,a));

I0 = @(eta,kappa) sqrt(2*pi) * C_price_Lewis(0, 1, eta, kappa,1, 1, 0);

pricesModel = {}; % Cell array to store model price functions
pricesMkt = [];   % Vector to store normalized market prices

% Loop through maturities (skipping the first one)
for i = start:last

    mon = strikes - fwdPrices(i);

    % Select only valid call and put options that are out-of-the-money
    maskCall = ~isnan(optionPrices(2*i-1,:)) & mon>0 & mon<30;
    maskPut = ~isnan(optionPrices(2*i,:)) & mon<0 & mon>-30;
    
    % Calculate normalized log-moneyness (standardized by ATM vol and sqrt time)
    xCall = (strikes(maskCall) - fwdPrices(i)) ./ (sigmaATM(i)*sqrt(ttm(i)));
    xPut = (strikes(maskPut) - fwdPrices(i)) ./ (sigmaATM(i)*sqrt(ttm(i)));

    disc_i = disc(i+1);  % Discount factor for this maturity
    
    % Store pricing function for this maturity 
    % Note: Many alternative pricing formulations are commented out
    pricesModel{end+1} = @(eta,kappa) C_price_Lewis( ...
        xCall,1, eta, kappa, ...
        1/I0(eta,kappa), 1, 0);
    pricesModel{end+1} = @(eta,kappa) C_price_Lewis( ...
        xPut, 1, eta, kappa, ...
        1/I0(eta,kappa), 1, 0) - (fwdPrices(i)-strikes(maskPut))/(sigmaATM(i)*sqrt(ttm(i))) ;
    
    % Normalize and store market prices
    pricesMkt = [pricesMkt, optionPrices(2*i-1,maskCall) ./ ...
        (sigmaATM(i)*sqrt(ttm(i))*disc_i)];
    pricesMkt = [pricesMkt, optionPrices(2*i,maskPut) ./ ...
        (sigmaATM(i)*sqrt(ttm(i))*disc_i)];
end

% Create vectorized pricing function that concatenates outputs from all maturities
pricesModelVec = @(p) cell2mat( cellfun(@(f) f(p(1),p(2)), ...
    pricesModel, 'Uni',false) );

% Define objective function (sum of squared pricing errors)
objective = @(p) sum( ( pricesModelVec(p) - pricesMkt ).^2 );

% Initial parameter guess
x0 = [eta0; k0]; 
% Parameter bounds
lb = [-Inf; 1e-6];  % Lower bounds: eta can be negative, kappa must be positive
ub = [Inf; 500];    % Upper bounds: large cap on kappa

% Optimization settings
opts = optimoptions('fmincon','Display','none');

% Run optimization to find parameters
[par,~] = fmincon(objective,x0,[],[],[],[],lb,ub,[],opts);

% Extract and return calibrated parameters
eta = par(1); 
kappa = par(2);
I0=I0(eta,kappa);

% Compute MSE (Mean Squared Error) for each maturity
MSE = zeros(1, last - start + 1); % Preallocate MSE vector
idx = 1; % Index for filling MSE vector

for i = start:last
    % Compute moneyness
    mon = strikes - fwdPrices(i);

    % Select valid out-of-the-money call and put options
    maskCall = ~isnan(optionPrices(2*i-1,:)) & mon > 0 & mon < 30;
    maskPut  = ~isnan(optionPrices(2*i,:)) & mon < 0 & mon > -30;

    % Standardized log-moneyness
    xCall = (strikes(maskCall) - fwdPrices(i)) ./ (sigmaATM(i)*sqrt(ttm(i)));
    xPut  = (strikes(maskPut) - fwdPrices(i)) ./ (sigmaATM(i)*sqrt(ttm(i)));

    % Discount factor for this maturity
    disc_i = disc(i+1);

    % Model prices (normalized)
    modelCall = C_price_Lewis(xCall, 1, eta, kappa, 1/I0, 1, 0);
    modelPut  = C_price_Lewis(xPut, 1, eta, kappa, 1/I0, 1, 0) - ...
                (fwdPrices(i) - strikes(maskPut)) ./ (sigmaATM(i)*sqrt(ttm(i)));

    % Market prices (normalized)
    marketCall = optionPrices(2*i-1, maskCall) ./ (sigmaATM(i)*sqrt(ttm(i))*disc_i);
    marketPut  = optionPrices(2*i, maskPut)     ./ (sigmaATM(i)*sqrt(ttm(i))*disc_i);

    % Total number of options used for this maturity
    nOptions = numel(modelCall) + numel(modelPut);

    % Compute Mean Squared Error (MSE) for this maturity
    mse_i = (sum((modelCall - marketCall).^2) + ...
             sum((modelPut  - marketPut).^2)) / nOptions;

    % Store MSE in output vector
    MSE(idx) = mse_i;
    idx=idx+1;

end
end