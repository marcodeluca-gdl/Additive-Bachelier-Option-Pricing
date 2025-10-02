function sigma = calibrateAddLogHConstant(strikes, fwdPrices, ttm, optionPrices)
% Calibrate additive logistic model parameters to market option prices
%
% Inputs:
%   strikes      - Strike prices vector
%   fwdPrices    - Forward prices for each maturity
%   ttm          - Time to maturity vector
%   optionPrices - Market option prices matrix
%
% Outputs:
%   sigma - Base volatility parameter

    pricesModel = {};
    pricesMkt = [];
    
    % Loop through each maturity
    for i = 2:length(ttm)
        % Calculate moneyness
        mon = strikes - fwdPrices(i);
        
        % Select valid options (non-NaN prices, reasonable moneyness range)
        validCall = find(~isnan(optionPrices(2*i-1,:)) & mon>0 & mon<30);
        validPut = find(~isnan(optionPrices(2*i,:)) & mon>-30 & mon<0);
        numValidCall = length(validCall);
        numValidPut = length(validPut);
        
        % Build model functions and market prices for valid options
        for j = 1:numValidCall
            pricesMkt = [pricesMkt, optionPrices(2*i-1,validCall(j))];
            pricesModel{end+1} = @(p) callPriceAddLog(p*ttm(i)^0.5, fwdPrices(i), strikes(validCall(j)));
        end

        for j = 1:numValidPut
            pricesMkt = [pricesMkt, optionPrices(2*i,validPut(j))];
            pricesModel{end+1} = @(p) putPriceAddLog(p*ttm(i)^0.5, fwdPrices(i), strikes(validPut(j)));
        end

    end
    
    % Create vectorized model function
    pricesModelVec = @(p) cellfun(@(f) f(p), pricesModel);
    
    % Define objective function (mean absolute percentage error)
    objective = @(p) sum( ( pricesModelVec(p) - pricesMkt ).^2 );
    
    % Set optimization options
    opts = optimoptions('fmincon', 'Display', 'none');
    
    % Initial guess and bounds
    p0 = 5;    
    lb = 0;        % Lower bound
    ub = 20;       % Upper bound
    
    % Optimize parameter
    [pstar, ~] = fmincon(objective, p0, [], [], [], [], lb, ub, [], opts);
    
    % Extract calibrated parameter
    sigma = pstar;
    
end