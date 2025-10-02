%% Additive Bachelier Option Pricing

tic
clear all;
close all;
clc;

addpath('datacalls')
addpath('dataputs')
addpath('Add_Bach _Calibration')
addpath('Bootstrap')
addpath('Pricing_Exotic')
addpath('Add_Log_Calibration')
addpath('Hedging')

%% Get the Data
callDir = 'datacalls';
putDir  = 'dataputs';
target  = 20200602; % our value date
[optionPrices, strikes] = buildOptionPrices(callDir, putDir, target);
formatData ='dd/MM/yyyy'; % pay attention to your computer settings
t0 = datetime('02-Jun-2020');
dates=getDates(t0); % futures expiries

%% Point 1 - Get Discounts
[disc, fwdPrices] = bootstrap(optionPrices, strikes);
zero_rates = zeroRates(dates(2:end), disc(2:end),dates(1)); % get zero rates from discounts
plot_disc(dates,disc)
plot_zerorates(dates,zero_rates)

%% Point 2 - Additive Bachelier Calibration

%% Step 1 - ATM volatilities
warning off

ACT_365 = 3;
ttm = yearfrac(dates(1), dates(2:end), ACT_365);
% extend optionPrices thanks to put call parity
optionPricesNew = put_call_parity(optionPrices(3:end,:), strikes, disc,fwdPrices);

% Get ATM volatilities
sigmaATM = ATMvols(optionPricesNew, strikes, fwdPrices, disc, ttm);
plotATMvolatilities(ttm,sigmaATM)
% Get Market Implied vols
market_vols = marketVols(optionPricesNew, strikes, fwdPrices, disc, ttm);

%% Step 2 - Additive Bachelier parameters
% Fixed model parameters
alpha = 1/3;
a = 1/3;
eta0 = 0.2;    % Initial guess for eta parameter
k0 = 1;        % Initial guess for kappa parameter
% Calibration with constant eta and k
[eta, kappa, I0, Mse_Add] = calibrateAddBach(optionPricesNew, sigmaATM, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, 2, length(ttm));
fprintf('Additive Bachelier\n')
fprintf('eta: %.4f\n', eta)
fprintf('kappa: %.4f\n', kappa)
fprintf('--------------\n')
% Rolling calibration: eta e kappa different for each ttm 
eta_t = zeros(1, length(ttm)-1);
kappa_t = zeros(1, length(ttm)-1);
I0_t = zeros(1, length(ttm)-1);
MSE_fixedttm=zeros(1, length(ttm)-1);
for i = 2:length(ttm)
    [eta_t(i-1), kappa_t(i-1), I0_t(i-1), MSE_fixedttm(i-1)] = calibrateAddBach(optionPricesNew, sigmaATM, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, i, i);
end
plotBachelierParameters(eta, kappa, eta_t, kappa_t, ttm);
%% Plot MSE of Additive Bachelier and Levy Fixed ttm
plotMSE(ttm, Mse_Add, MSE_fixedttm)
%% Plot of market implied vols vs Additive Bachelier implied vols
add_bach_vols = addBachVols(optionPrices, strikes, fwdPrices, disc, ttm, sigmaATM, eta, kappa,I0, a, alpha);
plotVolSmiles(strikes, fwdPrices, market_vols,add_bach_vols, ttm);
%% Point 3 - Additive Logistic Calibration

[sigma, H] = calibrateAddLog(strikes, fwdPrices, ttm, optionPrices); % Calibrate sigma and H
sigma_H_const = calibrateAddLogHConstant(strikes, fwdPrices, ttm, optionPrices); % Calibrate sigma with H=0.5
fprintf('Additive Logistic\n')
fprintf('sigma: %.4f\n', sigma)
fprintf('H: %.4f\n', H)
fprintf('--------------\n')

s = sigma.*ttm.^H;
sigmaATMLog = s.*sqrt((2*pi./ttm))*log(2); % H from calibration
s_const = sigma_H_const.*ttm.^0.5;
sigmaATMLogHConst = s_const.*sqrt((2*pi./ttm))*log(2); % H = 0.5 fixed in advance
figure
plot(ttm, sigmaATM, '*-');
hold on;
plot(ttm, sigmaATMLog, '*-');
plot(ttm, sigmaATMLogHConst, '*-');
grid on;
legend('ATM-Vol-Bachelier', 'ATM-Vol-Logistic', 'ATM-Vol-Logistic-H=0.5')
%% Plot of market implied vols vs Additive Logistic and Additive BAchelier implied vols
add_log_vols = addLogVols(optionPrices, strikes,fwdPrices, disc, ttm,s);
plotVolSmilesComplete(strikes, fwdPrices, market_vols, add_bach_vols, add_log_vols, ttm);
%% Point 4 - Pricing Exotic Option
sigmat = sigmaATM/I0;
ttmlong=[0;ttm];
Nsim=1e7;
increment=-100:0.1:100;
Notional=20e6;
strike=fwdPrices(end);
[priceBac, CIBac] = pricingExotic(strike, disc,ttmlong, fwdPrices, Nsim, increment, eta, kappa,sigmat,alpha, 0 , 0, Notional, 'Bac');
[priceLog, CILog] = pricingExotic(strike, disc,ttmlong, fwdPrices, Nsim, increment, 0, 0,0,0, sigma/100 , H, Notional, 'Log');
fprintf('Exotic price\n')
fprintf('priceBac: %d\n', priceBac)
fprintf('CI Bachelier: [%d; %d]\n',CIBac(1), CIBac(2));
fprintf('priceLog: %d\n', priceLog)
fprintf('CI Logistic: [%d; %d]\n',CILog(1), CILog(2));
fprintf('--------------\n')
%% Point 5 - Hedging

%% Step 1 - Hedge the volatility
[deltaPriceSigma, sigmat1] = deltaPriceVolatility(sigmaATM,optionPricesNew, strikes, strike,fwdPrices,...
                    ttm, ttmlong, disc, alpha, a, eta0, k0,  Nsim, increment,priceBac, Notional);
idxTop = paretoSelection(abs(deltaPriceSigma), 0.8);
deltaSigmaTop = deltaPriceSigma(idxTop); % most significant price changes
callVolatility = callVolatilityBump(idxTop, fwdPrices, strikes, ...
    optionPricesNew, sigmat1, ttm, disc, deltaSigmaTop);
% Now we hedge the delta with futures
futuresInfoVol = hedgeDeltaVolatility(callVolatility, fwdPrices, sigmat, ttm, disc);
%% Step 2 - Hedge the skew
bullInfo = hedgeEtaBullSpread(strike, strikes, optionPricesNew, sigmaATM, ttm,ttmlong, disc, ...
        eta, kappa, sigmat, alpha,fwdPrices, Nsim, increment, Notional, priceBac);
deltaBull = bullInfo(1) * ( deltaOption(bullInfo(2), fwdPrices(end), sigmat(end), ttm(end), disc(end),'Nan') ...
    - deltaOption(bullInfo(3),fwdPrices(end), sigmat(end), ttm(end), disc(end),'Nan') );
nFutureEta  = -deltaBull;
numberFutureEta = ceil(nFutureEta / fwdPrices(end));
%% Step 3 - Hedge the vol-vol
strangleInfo = hedgeKappaStrangle(strike, strikes, optionPricesNew, sigmaATM, ttm,ttmlong, disc, ...
        eta, kappa, sigmat, alpha, fwdPrices, Nsim, increment, Notional,priceBac);
deltaStr = strangleInfo(1) * (  deltaOption(strangleInfo(2),fwdPrices(end),sigmat(end),ttm(end),disc(end),'Nan') ...
                              + deltaOption(strangleInfo(3),fwdPrices(end),sigmat(end),ttm(end),disc(end),'Put')  );
nFutureKappa  = -deltaStr;
numberFutureKappa = ceil(nFutureKappa/fwdPrices(end));
%% Step 4 - Evaluate P&L
CostHedging = 6*1e-4 * (sum(abs(callVolatility(1,:)))+abs(bullInfo(1))+abs(strangleInfo(1)))+...
              1*1e-4 *  (sum(abs(futuresInfoVol(1,:)))+abs(nFutureEta)+abs(nFutureKappa));
ptfValue = priceBac + sum(ceil(callVolatility(1,:))) + sum(ceil(futuresInfoVol(1,:))) + ceil(bullInfo(1)) +...
          ceil(nFutureEta) + ceil(strangleInfo(1)) + ceil(nFutureKappa);
targets = targetPL();
datesPL = getdatesPL();
[PLNoHedging, PLHedging, PLHedgingCost] = evaluatePL(callDir, putDir, targets, datesPL, alpha, ...
                                                   a, eta0, k0, Nsim, increment, Notional, callVolatility,...
                                                   futuresInfoVol, bullInfo, numberFutureEta, strangleInfo,...
                                                   numberFutureKappa, ptfValue, CostHedging, priceBac);
plotPL(PLNoHedging, PLHedging, PLHedgingCost,datesPL)
fprintf('Hedging Quality\n')
fprintf('PL 9-Jun: %d\n', PLHedgingCost(5));
fprintf('PL 16-Jun: %d\n', PLHedgingCost(10));
toc


