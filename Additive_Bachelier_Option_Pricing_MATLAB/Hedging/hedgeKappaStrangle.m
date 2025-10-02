function strangleInfo = hedgeKappaStrangle( ...
    strike, strikes, optionPricesNew, sigmaATM, ttm,ttmlong, disc, ...
    eta, kappa, sigmat, alpha, ...
    fwdPrices, Nsim, increment, Notional, ...
    priceBac)
%--------------------------------------------------------------------------
% Hedge the κ-(kappa) sensitivity of an exotic position using a 10-wide
% Bachelier strangle (long call K+10, long put K-10).
%
% OUTPUT (4 × 1 column vector)
% strangleInfo(1) = notional (number of strangles to BUY if >0, SELL if <0)
% strangleInfo(2) = strike of the call leg (K + 10)
% strangleInfo(3) = strike of the put leg (K − 10)
% strangleInfo(4) = number of strangle contracts
%
% Assumes identical contract-size for exotic, call, and put.
%--------------------------------------------------------------------------

bp=1e-4;

% ---------- 1) Shock κ by +10 bp ---------------------------------------
kappaShocked = kappa + 100 * bp;

% ---------- 2) Re-price the exotic payoff -------------------------------
[priceShocked, ~] = pricingExotic(strike, disc, ttmlong, fwdPrices, ...
    Nsim, increment, eta, kappaShocked, ...
    sigmat, alpha, 0, 0, Notional, 'Bac');

deltaPriceKappa = priceShocked - priceBac; % € change for +10 bp

% ---------- 3) Locate call (K+10) and put (K-10) strikes -----------------
% Call leg
distCall = abs(strikes - strike - 10);
distCall(isnan(optionPricesNew(end-1,:))) = Inf;
[~, idxCall] = min(distCall);
Kcall = strikes(idxCall);

% Put leg
distPut = abs(strikes - strike + 10);
distPut(isnan(optionPricesNew(end,:))) = Inf;
[~, idxPut] = min(distPut);
Kput = strikes(idxPut);

% ---------- 4) Lewis-FFT Bachelier helpers -------------------------------
psi = @(u,kap) (1./kap) .* (1-alpha)./alpha .* ...
    (1 - (1 + (u.*kap)./(1-alpha)).^alpha);

phi = @(u,sig,eta_,kap,t) ...
    exp( psi( 1i*u.*eta_.*sig.*sqrt(t) + 0.5*u.^2 .* sig.^2 .* t, kap) ...
         + 1i*u.*eta_.*sig.*sqrt(t) );

M = 15; % FFT grid size
dz = 0.01; % FFT grid spacing
a = 1/3; % damping

C_Bac = @(x,t,eta_,kap,sig,B,Ra) ...
    B .* ( Ra + (exp(-a*x)/(2*pi)) .* ...
    FFTBachelier(@(u) phi(u,sig,eta_,kap,t), M, dz, x, a) );

I0 = @(eta_,kap) sqrt(2*pi) * C_Bac(0, 1, eta_, kap, 1, 1, 0);

% ---------- 5) Price the strangle (base κ) -------------------------------
FwdScale = sigmaATM(end)*sqrt(ttm(end)); % for normalised moneyness

call = C_Bac((Kcall-strike)/FwdScale, 1, eta, kappa, ...
    1/I0(eta,kappa), 1, 0) * FwdScale * disc(end);

put = C_Bac((Kput-strike)/FwdScale, 1, eta, kappa, ...
    1/I0(eta,kappa), 1, 0) * FwdScale * disc(end) ...
    - disc(end)*(strike - Kput); % Bachelier put-call parity

strangle = call + put;

% ---------- 6) Price the strangle (shocked κ) ----------------------------
callS = C_Bac((Kcall-strike)/FwdScale, 1, eta, kappaShocked, ...
    1/I0(eta,kappaShocked), 1, 0) * FwdScale * disc(end);

putS = C_Bac((Kput-strike)/FwdScale, 1, eta, kappaShocked, ...
    1/I0(eta,kappaShocked), 1, 0) * FwdScale * disc(end) ...
    - disc(end)*(strike - Kput);

strangleShocked = callS + putS;

% ---------- 7) Sensitivity & required notional --------------------------
deltaStrangleKappa = strangleShocked - strangle; % € change for +10 bp
nStrangle = -deltaPriceKappa / deltaStrangleKappa;

% ---------- 8) Pack and return ------------------------------------------
strangleInfo = [nStrangle; Kcall; Kput; ceil(nStrangle/(optionPricesNew(2*length(ttm)-1,idxCall)+optionPricesNew(2*length(ttm),idxPut)))];

end
