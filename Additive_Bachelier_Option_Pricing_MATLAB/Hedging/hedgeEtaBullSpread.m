function bullInfo = hedgeEtaBullSpread( ...
    strike, strikes, optionPricesNew, sigmaATM, ttm,ttmlong, disc, ...
    eta, kappa, sigmat, alpha,...
    fwdPrices, Nsim, increment, Notional, ...
    priceBac)
%----------------------------------------------------------------------
% Hedge the η-(eta) sensitivity of an exotic position using a 20-wide
% Bachelier bull spread (long K-20, short K+20).
%
% OUTPUT (4 × 1 column vector)
% bullInfo(1) = notional (number of bull-spreads to BUY if >0, SELL if <0)
% bullInfo(2) = lower-strike call (K − 20)
% bullInfo(3) = upper-strike call (K + 20)
% bullInfo(4) = number of bull spread contracts

bp=1e-4;

% ---------- 1) Shock η by +1 bp ------------------------------------
etaShocked = eta + 100 * bp;

% ---------- 2) Re-price the exotic payoff ----------------------------
[priceShocked, ~] = pricingExotic(strike, disc, ttmlong, fwdPrices, ...
    Nsim, increment, etaShocked, ...
    kappa, sigmat, alpha, ...
    0, 0, Notional, 'Bac');

deltaPriceEta = priceShocked - priceBac; % € change of exotic for +1 bp

% ---------- 3) Locate the strikes K-20 and K+20 --------------------
distBelow = abs(strikes - strike + 10);
distBelow(isnan(optionPricesNew(end-1,:))) = Inf;
[~, idxBelow] = min(distBelow);
K1 = strikes(idxBelow); % lower-leg strike

distAbove = abs(strikes - strike - 40);
distAbove(isnan(optionPricesNew(end-1,:))) = Inf;
[~, idxAbove] = min(distAbove);
K2 = strikes(idxAbove); % upper-leg strike

% ---------- 4) Helper handles for the Lewis-FFT Bachelier price -------
psi = @(u,kap) (1./kap) .* (1-alpha)./alpha .* ...
    (1 - (1 + (u.*kap)./(1-alpha)).^alpha);

phi = @(u,sig,eta_,kap,t) ...
    exp( psi( 1i*u.*eta_.*sig.*sqrt(t) + 0.5*u.^2 .* sig.^2 .* t, kap) ...
         + 1i*u.*eta_.*sig.*sqrt(t) );

M = 15; % FFT grid size
dz = 0.01; % FFT grid spacing
a = 1/3;

C_Bac = @(x,t,eta_,kap,sig,B,Ra) ...
    B .* ( Ra + (exp(-a*x)/(2*pi)) .* ...
    FFTBachelier(@(u) phi(u,sig,eta_,kap,t), M, dz, x, a) );

I0 = @(eta_,kap) sqrt(2*pi) * C_Bac(0, 1, eta_, kap, 1, 1, 0); % normaliser

% Price the two call legs (base η)
C1 = C_Bac((K1-strike)/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, ...
    1/I0(eta,kappa), 1, 0) ...
    * sigmaATM(end)*sqrt(ttm(end))*disc(end);

C2 = C_Bac((K2-strike)/(sigmaATM(end)*sqrt(ttm(end))), 1, eta, kappa, ...
    1/I0(eta,kappa), 1, 0) ...
    * sigmaATM(end)*sqrt(ttm(end))*disc(end);

bullSpread = C1 - C2; % long K1, short K2

% Price the two call legs (shocked η)
C1s = C_Bac((K1-strike)/(sigmaATM(end)*sqrt(ttm(end))), 1, etaShocked, kappa, ...
    1/I0(etaShocked,kappa), 1, 0) ...
    * sigmaATM(end)*sqrt(ttm(end))*disc(end);

C2s = C_Bac((K2-strike)/(sigmaATM(end)*sqrt(ttm(end))), 1, etaShocked, kappa, ...
    1/I0(etaShocked,kappa), 1, 0) ...
    * sigmaATM(end)*sqrt(ttm(end))*disc(end);

bullSpreadShocked = C1s - C2s;

% ---------- 5) Sensitivity of the bull-spread to η --------------------
deltaBullEta = bullSpreadShocked - bullSpread; % € change for +1 bp

% ---------- 6) Notional required to offset exotic η-risk --------------
nBull = - deltaPriceEta / deltaBullEta; % (+) buy, (-) sell

% ---------- 7) Pack the output ----------------------------------------
bullInfo = [nBull; K1; K2; ceil(nBull/(optionPricesNew(2*length(ttm)-1,idxBelow) + optionPricesNew(2*length(ttm)-1,idxAbove)))];

end