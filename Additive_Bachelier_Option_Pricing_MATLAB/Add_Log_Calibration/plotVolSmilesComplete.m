function plotVolSmilesComplete(strikes, fwdPrices, marketVols, addBachVols, addLogVols, ttm)
% Plot market vs Additive Bachelier vs Additive Logistic normal-implied volatility smiles
%
% INPUTS:
% strikes - Vector of strike prices
% fwdPrices - Vector of forward prices for each maturity
% marketVols - Matrix of market implied volatilities (maturities x strikes)
% addBachVols - Matrix of Additive Bachelier model implied volatilities (maturities x strikes)
% addLogVols - Matrix of Additive Logistic model implied volatilities (maturities x strikes)
% ttm - Vector of times to maturity
%
% OUTPUT:
% Creates subplot figure comparing market (red) vs Additive Bachelier (blue) vs Additive Logistic (green) vol smiles
% First maturity is skipped in the plots

Mtotal = numel(fwdPrices); % Total number of maturities
Mplots = Mtotal - 1; % Number of plots (skip first maturity)

% Create figure with multiple subplots
figure('Name','Market vs Additive Bachelier vs Additive Logistic Vols', ...
       'NumberTitle','off','Color','w');

% --- fixed 2 rows × 4 columns (for 8 plots) ------------------------------
T = tiledlayout(2, 4, 'Padding','compact', 'TileSpacing','compact');
title(T,'Market (red) vs Additive Bachelier (blue) vs Additive Logistic (green)');

% Loop through maturities (skip first)
for i = 2:Mtotal
    nexttile;
    
    % Calculate moneyness (K - F)
    x = strikes - fwdPrices(i);
    yMkt = marketVols(i,:);
    yAddBach = addBachVols(i,:);
    yAddLog = addLogVols(i,:);
    
    % Plot market volatilities (red line)
    vM = ~isnan(yMkt); % Valid market data points
    if any(vM)
        plot(x(vM), yMkt(vM), 'r-', 'LineWidth',1.2, ...
             'DisplayName','Market');
    end
    hold on
    
    % Plot Additive Bachelier model volatilities (blue line)
    vB = ~isnan(yAddBach); % Valid Additive Bachelier data points
    if any(vB)
        plot(x(vB), yAddBach(vB), 'b-', 'LineWidth',1.2, ...
             'DisplayName','Add Bach');
    end
    
    % Plot Additive Logistic model volatilities (green line)
    vL = ~isnan(yAddLog); % Valid Additive Logistic data points
    if any(vL)
        plot(x(vL), yAddLog(vL), 'g-', 'LineWidth',1.2, ...
             'DisplayName','Add Log');
    end
    
    hold off
    
    % Format subplot
    grid on;
    xlim([-30 30]);
    xlabel('Moneyness (K - F)');
    ylabel('Normal-implied volatility');
    title(sprintf('TTM = %.2f yr', ttm(i)));
    legend('Location','best');
end
end