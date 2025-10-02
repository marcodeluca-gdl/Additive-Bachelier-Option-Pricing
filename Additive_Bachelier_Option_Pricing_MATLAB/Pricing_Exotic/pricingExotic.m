function [price, CI] = pricingExotic(strike, disc, ttm, fwdPrices, Nsim, increment, eta, kappa, sigmat, alpha, sigmaLog, H, Notional, flag)
% Price exotic option using Monte Carlo with additive models
%
% Inputs:
%   strike     - Strike price
%   disc       - Discount factors
%   ttm        - Time to maturity vector
%   fwdPrices  - Forward prices
%   Nsim       - Number of Monte Carlo simulations
%   increment  - Grid for CDF inversion
%   eta, kappa, sigmat, alpha - Bachelier model parameters
%   sigmaLog, H - Log-additive model parameters
%   Notional   - Contract notional
%   flag       - Model type ('Bac' or 'Log')
%
% Output:
%   price - Option price
%   CI - confidence interval

    rng(42);  % Set random seed for reproducibility
    
    % Initialize arrays
    cumulative = zeros(length(ttm)-1, length(increment));
    U = rand(length(ttm)-1, Nsim);  % Random numbers for simulation
    incr = zeros(length(ttm)-1, Nsim);
    simFwd = zeros(length(fwdPrices)-1, Nsim);
    payoff = zeros(Nsim, 1);
    
    % Bachelier additive model
    if flag == 'Bac'
        % Define psi function for Bachelier model
        psi_fun = @(u, kappa) (1./kappa) .* (1-alpha)./alpha .* ...
            (1 - (1 + (u.*kappa)./(1-alpha)).^alpha);
        
        % Characteristic function for Bachelier additive model
        phiAddBac = @(u, sigma, eta, kappa, t) ...
            exp(psi_fun(1i*u.*eta.*sigma.*sqrt(t) + ...
            0.5*u.^2 .* sigma.^2 .* t, kappa) + 1i*u.*eta.*sigma.*sqrt(t));
        
        % Compute cumulative distributions for each time step
        for i = 1:length(ttm)-1
            if i == 1
                phiAddBac1 = @(u) 1;
            else
                phiAddBac1 = @(u) phiAddBac(u, sigmat(i-1), eta, kappa, ttm(i));
            end
            phiAddBac2 = @(u) phiAddBac(u, sigmat(i), eta, kappa, ttm(i+1));
            phiAddBacQuot = @(u) phiAddBac2(u)./phiAddBac1(u);
            cumulative(i,:) = lewisFormula(phiAddBacQuot, -0.01, 0, increment, ttm(i), ttm(i+1));
        end
    end
    
    % Log-additive model
    if flag == 'Log'
        % Define b function for log-additive model
        b = @(sigma, H, t) (1-exp(-t.*sigma.^(1/H))).^H;
        
        % Characteristic function for log-additive model
        phiAddLog = @(u, sigma, H, t) (1-b(sigmaLog, H, t)) ...
            * beta_complex(1+(1i*u-1)*b(sigmaLog,H,t), 1-1i*u*b(sigmaLog,H,t));
        
        % Compute cumulative distributions for each time step
        for i = 1:length(ttm)-1
            if i == 1
                phiAddLog1 = @(u) 1;
            else
                phiAddLog1 = @(u) phiAddLog(u, sigmaLog, H, ttm(i));
            end
            phiAddLog2 = @(u) phiAddLog(u, sigmaLog, H, ttm(i+1));
            phiAddLogQuot = @(u) phiAddLog2(u)./phiAddLog1(u);
            cumulative(i,:) = lewisFormula(phiAddLogQuot, 0.01, 1, increment, ttm(i), ttm(i+1));
        end
    end

    % truncate the cumulative matrix
    cumulative = max(min(cumulative, 1), 0);
    for r = 1:size(cumulative,1)
        cumulative(r,:) = cummax(cumulative(r,:));  % monotonicità riga per riga
    end
    
    % Generate increments using inverse CDF method
    for i = 1:length(ttm)-1 
        x=cumulative(i,:);
        mask = [true, diff(x) > 0];
        incr(i,:) = interp1(x(mask), increment(mask), U(i,:), 'spline');
    end
    
    % Simulate forward prices
    for i = 1:length(fwdPrices)
        if i == 1
            simFwd(i,:) = fwdPrices(end) + incr(i,:);
        else
            simFwd(i,:) = simFwd(i-1,:) + incr(i,:);
        end
    end
    
    % Calculate average of simulated forwards
    sumSimFwd = mean(simFwd, 1);
    
    % Compute payoff (binary option)
    for i = 1:Nsim
        if strike > sumSimFwd(i)
            payoff(i) = 1;
        end
    end
    
    % Calculate discounted expected payoff
    price = Notional * disc(end) * payoff;
    stdev = std(price);
    CI = [mean(price) - norminv(0.995)*stdev/sqrt(Nsim), mean(price) + norminv(0.995)*stdev/sqrt(Nsim) ];
    price = mean(price);
    
end
