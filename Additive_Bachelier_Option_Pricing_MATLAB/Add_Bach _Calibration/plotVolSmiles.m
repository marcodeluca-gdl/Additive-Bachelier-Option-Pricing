function plotVolSmiles(strikes, fwdPrices, marketVols, modelVols, ttm)
% Plot market vs Additive Bachelier normal-implied volatility smiles
%
% INPUTS:
%   strikes     - Vector of strike prices
%   fwdPrices   - Vector of forward prices for each maturity
%   marketVols  - Matrix of market implied volatilities (maturities x strikes)
%   modelVols   - Matrix of model implied volatilities (maturities x strikes)
%   ttm         - Vector of times to maturity
%
% OUTPUT:
%   Creates subplot figure comparing market (red) vs model (blue) vol smiles
%   First maturity is skipped in the plots

Mtotal = numel(fwdPrices);  % Total number of maturities
Mplots = Mtotal - 1;        % Number of plots (skip first maturity)

% Create figure with multiple subplots
figure('Name','Market vs Additive Bachelier Vols', ...
       'NumberTitle','off','Color','w');

% --- fixed 2 rows × 4 columns (for 8 plots) ------------------------------
T = tiledlayout(2, 4, 'Padding','compact', 'TileSpacing','compact');
title(T,'Market (red) vs Additive Bachelier (blue)');

% Loop through maturities (skip first)
for i = 2:Mtotal
    nexttile;
    
    % Calculate moneyness (K - F)
    x = strikes - fwdPrices(i);
    yMkt = marketVols(i,:);
    yModel = modelVols(i,:);
    
    % Plot market volatilities (red line)
    vM = ~isnan(yMkt);  % Valid market data points
    if any(vM)
        plot(x(vM), yMkt(vM), 'r-', 'LineWidth',1.2, ...
             'DisplayName','Market');
    end
    hold on
    
    % Plot model volatilities (blue line)
    vP = ~isnan(yModel);  % Valid model data points
    if any(vP)
        plot(x(vP), yModel(vP), 'b-', 'LineWidth',1.2, ...
             'DisplayName','Add Bach');
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