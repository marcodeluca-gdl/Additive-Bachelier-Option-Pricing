function plotPL(PLNoHedging, PLHedging, PLHedgingCost, datesPL)
% plotPL  Plots Profit & Loss over time for three scenarios
%         (no hedge, hedge, hedge net of costs) on one figure.
%
%   plotPL(PLNoHedging, PLHedging, PLHedgingCost, datesPL)
%
%   INPUT
%     PLNoHedging   : Nx1 vector – P/L without hedging
%     PLHedging     : Nx1 vector – P/L with hedging
%     PLHedgingCost : Nx1 vector – P/L with hedging after costs
%     datesPL       : Nx1 vector – dates (datetime or Excel serial)
%
%   Quick example
%       plotPL(noH, hedg, hedgCost, myDates)

    % -------------------------------------------------------------
    % 1) Sanity checks
    % -------------------------------------------------------------
    sameLen = numel(PLNoHedging)==numel(PLHedging) && ...
              numel(PLHedging)==numel(PLHedgingCost) && ...
              numel(PLHedgingCost)==numel(datesPL);
    if ~sameLen
        error('All vectors must have identical length.');
    end

    % Convert dates to datetime if needed
    if ~isdatetime(datesPL)
        datesPL = datetime(datesPL, 'ConvertFrom', 'datenum', ...
                                      'Format', 'dd/MM/yyyy');
    end

    % -------------------------------------------------------------
    % 2) Plot
    % -------------------------------------------------------------
    figure('Color','w', 'Position',[100 100 900 450]);
    hold on; box on; grid on;

    p1 = plot(datesPL, PLNoHedging,   '-',  'LineWidth',1.8);
    p2 = plot(datesPL, PLHedging,     '--', 'LineWidth',1.8);
    p3 = plot(datesPL, PLHedgingCost, ':',  'LineWidth',1.8);

    % Simple soft-tone palette
    p1.Color = [0.25 0.47 0.85];   % blue
    p2.Color = [0.90 0.45 0.13];   % orange
    p3.Color = [0.30 0.70 0.30];   % green

    % -------------------------------------------------------------
    % 3) Formatting
    % -------------------------------------------------------------
    title('P/L evolution comparison', 'FontSize',14, 'FontWeight','bold');
    ylabel('P/L (€)','FontSize',12);
    xlabel('Date','FontSize',12);

    legend([p1 p2 p3], ...
           {'No Hedging', 'Hedging', 'Hedging (cost included)'}, ...
           'Location','best', 'FontSize',10);

    % Readable date ticks (1st of each month)
    ax = gca;
    ax.XAxis.TickLabelFormat = 'dd-MMM-yyyy';
    ax.XAxis.TickValues      = datesPL(1):calmonths(1):datesPL(end);
    ax.XTickLabelRotation    = 30;
    ax.FontSize = 10;

    % Zero reference line
    yline(0,'k-','LineWidth',1);

    hold off;% -------------------------------------------------------------
    % 3) Formatting
    % -------------------------------------------------------------
    title('P/L evolution comparison', 'FontSize',14, 'FontWeight','bold');
    ylabel('P/L (€)','FontSize',12);
    xlabel('Date','FontSize',12);
    
    legend([p1 p2 p3], ...
           {'No Hedging', 'Hedging', 'Hedging (cost included)'}, ...
           'Location','best', 'FontSize',10);
    
    % ---- x-axis ticks -------------------------------------------------------
    ax          = gca;
    numXTicks   = 10;                         % <-- change this to taste
    tickIdx     = round(linspace(1, numel(datesPL), numXTicks));
    ax.XTick    = datesPL(tickIdx);
    ax.XAxis.TickLabelFormat = 'dd-MMM-yyyy';
    ax.XTickLabelRotation    = 30;
    ax.FontSize = 10;
    
    % ---- zero reference line (hidden from legend) --------------------------
    yline(0,'k-','LineWidth',1,'HandleVisibility','off');
    
    hold off;

end

