import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from datetime import datetime
from typing import Tuple
import matplotlib.dates as mdates

from Pricing_Exotic import pricing_exotic
from Add_Bach_Calibration import *
from Add_Log_Calibration import calibrate_add_log
from Bootstrap import *

def delta_price_volatility(sigmaATM, optionPricesNew, strikes, strike, fwdPrices, ttm, ttmlong, disc, alpha, a, eta0, k0, Nsim, increment, priceBac, Notional):
    # Define shock size: 1 basis point = 0.0001 = 0.01%
    bp = 1e-4
    
    deltaPriceSigma = np.zeros(len(sigmaATM))
    sigmatnew = np.zeros(len(sigmaATM))

    # Loop through each volatility parameter to calculate individual sensitivities
    for i in range(len(sigmaATM)):
        # Create a copy of the original volatility vector
        sigmaATMShocked = np.copy(sigmaATM)
        
        # Apply shock to the i-th volatility parameter
        sigmaATMShocked[i] += 100 * bp  # 1% shock
        
        # Recalibrate the model with the shocked volatility
        eta1, kappa1, I0_1 = calibrate_add_bach(optionPricesNew, sigmaATMShocked, strikes, fwdPrices, ttm, disc, alpha, a, eta0, k0, 1, len(ttm))
        #print(f"Calibrated eta: {eta1}, kappa: {kappa1}, I0: {I0_1}")
        # Adjust volatility term structure using calibration factor
        sigmat1 = sigmaATMShocked / I0_1
        #print(f"Adjusted sigmat1: {sigmat1}")
        sigmatnew[i] = sigmat1[i]
        
        # Price the exotic option with shocked parameters using the Bachelier model
        temp, _ = pricing_exotic(strike, disc, ttmlong, fwdPrices, Nsim, increment, eta1, kappa1, sigmat1, alpha, 0, 0, Notional, 'Bac')
        #print(temp)
        # Calculate the price sensitivity (delta)
        deltaPriceSigma[i] = temp - priceBac
        #print(f"Delta price sensitivity for sigma[{i}]: {deltaPriceSigma[i]}")

    return deltaPriceSigma, sigmatnew


def pareto_selection(values, threshold=0.8):
    # Pareto selection (80/20 rule) to identify the subset of elements that contribute to the specified cumulative fraction of the total
    sorted_indices = np.argsort(values)[::-1]
    cumulative_sum = np.cumsum(values[sorted_indices])
    total = np.sum(values)
    
    idx = np.where(cumulative_sum / total >= threshold)[0][0]
    return sorted_indices[:idx + 1]



def call_volatility_bump(idxTop, fwdPrices, strikes, optionPricesNew, sigmat1, ttm, disc, deltaSigmaTop):
    """
    Constructs hedging call options for volatility exposure.
    
    Parameters:
        idxTop (array-like): Indices of top volatility risk contributors.
        fwdPrices (np.ndarray): Forward prices vector.
        strikes (np.ndarray): Available strike prices.
        optionPricesNew (np.ndarray): Market option prices matrix.
        sigmat1 (np.ndarray): Volatility term structure.
        ttm (np.ndarray): Time to maturity vector.
        disc (np.ndarray): Discount factors.
        deltaSigmaTop (np.ndarray): Volatility sensitivities to be hedged.
    
    Returns:
        callInfo (np.ndarray): 5 × nCall matrix.
            Row 0 → Call option notional amounts
            Row 1 → Associated strike prices
            Row 2 → Time-to-maturity
            Row 3 → Time-to-maturity indexes
            Row 4 → Number of calls
    """
    nCall = len(idxTop)
    bp = 1e-4

    notionalCall = np.zeros(nCall)
    strikeCall = np.zeros(nCall)
    ttmCall = np.zeros(nCall)
    ttmIdx = np.zeros(nCall)
    numCall = np.zeros(nCall)

    for i in range(nCall):
        idx = idxTop[i]
        fwdPrice = fwdPrices[idx]

        # Find closest strike to fwdPrice, excluding NaN option prices
        dist = np.abs(strikes - fwdPrice)
        dist[np.isnan(optionPricesNew[2 * idx, :])] = np.inf
        strike_idx = np.argmin(dist)
        closestStrike = strikes[strike_idx]

        # Compute Vega
        vegaCall = VegaOption(closestStrike, fwdPrice, sigmat1[idx], ttm[idx], disc[idx + 1]) * 100 * bp

        # Hedge notional
        notionalCall[i] = -deltaSigmaTop[i] / vegaCall

        strikeCall[i] = closestStrike
        ttmCall[i] = ttm[idx]
        ttmIdx[i] = idx
        numCall[i] = notionalCall[i] / optionPricesNew[2 * idx, strike_idx]

    # Stack the results into a 5 × nCall array
    callInfo = np.vstack([
        notionalCall,
        strikeCall,
        ttmCall,
        ttmIdx,
        np.ceil(numCall)
    ])

    return callInfo


def hedge_delta_volatility(callInfo, fwdPrices, sigmat, ttm, disc):
    nCall = callInfo.shape[1]
    futureNotional = np.zeros(nCall)
    futureTTM = callInfo[2, :]
    TtmIdx = callInfo[3, :]
    nFutures = np.zeros(nCall)

    for i in range(nCall):
        N_opt = callInfo[0, i]
        K = callInfo[1, i]
        tau = callInfo[2, i]

        idx = np.where(np.abs(ttm - tau) < 1e-10)[0]
        k = idx[0] if idx.size > 0 else None

        F = fwdPrices[k]
        vol = sigmat[k]
        DF = disc[k]
        
        delta_call = deltaOption(K, F, vol, tau, DF, 'Nan')  # Delta option is a placeholder function
        portfolio_delta = N_opt * delta_call
        
        futureNotional[i] = -portfolio_delta
        nFutures[i] = futureNotional[i] / fwdPrices[int(TtmIdx[i])]

    futureInfo = np.vstack([futureNotional, futureTTM, TtmIdx, np.ceil(nFutures)])
    return futureInfo

def VegaOption(strike, fwdPrice, sigmat, ttm, disc):
    # Calculate the moneyness (y) based on the provided parameters
    y = (strike - fwdPrice) / (sigmat * np.sqrt(ttm))
    
    # Vega is the discount factor multiplied by the normal PDF of the moneyness
    vega = disc * norm.pdf(-y / sigmat)  # Normal PDF from scipy.stats.norm
    return vega

def deltaOption(strike, fwdPrice, sigmat, ttm, disc, flag):
    # Calculate the moneyness (y) based on the provided parameters
    y = (strike - fwdPrice) / (sigmat * np.sqrt(ttm))
    
    # Delta is the discount factor multiplied by the normal CDF of the moneyness
    delta = disc * norm.cdf(-y / sigmat)  # Normal CDF from scipy.stats.norm
    
    # If the option is a put, adjust the delta
    if flag == 'put':
        delta = delta - disc
        
    return delta



def hedge_eta_bull_spread(strike, strikes, option_prices_new, sigma_atm, ttm, ttmlong, disc,
                          eta, kappa, sigmat, alpha, fwd_prices, n_sim, increment, notional, price_bac):
    """
    Hedge the η-(eta) sensitivity of an exotic position using a 20-wide Bachelier bull spread.
    (long K-20, short K+20).

    OUTPUT:
        bullInfo: (4 × 1 column vector)
            bullInfo[0] = notional (number of bull-spreads to BUY if >0, SELL if <0)
            bullInfo[1] = lower-strike call (K − 20)
            bullInfo[2] = upper-strike call (K + 20)
            bullInfo[3] = number of bull spread contracts
    """
    bp = 1e-4  # Shock size (1 basis point)

    # ---------- 1) Shock η by +1 bp ------------------------------------
    eta_shocked = eta + 1 * bp

    # ---------- 2) Re-price the exotic payoff ----------------------------
    price_shocked, _ = pricing_exotic(strike, disc, ttmlong, fwd_prices, n_sim, increment, eta_shocked,
                                      kappa, sigmat, alpha, 0, 0, notional, 'Bac')

    delta_price_eta = price_shocked - price_bac  # € change of exotic for +1 bp

    # ---------- 3) Locate the strikes K-20 and K+20 --------------------
    dist_below = np.abs(strikes - strike + 20)
    dist_below[np.isnan(option_prices_new[-2, :])] = np.inf
    idx_below = np.argmin(dist_below)
    k1 = strikes[idx_below]  # lower-leg strike

    dist_above = np.abs(strikes - strike - 20)
    dist_above[np.isnan(option_prices_new[-2, :])] = np.inf
    idx_above = np.argmin(dist_above)
    k2 = strikes[idx_above]  # upper-leg strike

    # ---------- 4) Helper handles for the Lewis-FFT Bachelier price -------
    def psi(u, kap):
        return (1 / kap) * (1 - alpha) / alpha * (1 - (1 + (u * kap) / (1 - alpha)) ** alpha)

    def phi(u, sig, eta_, kap, t):
        return np.exp(psi(1j * u * eta_ * sig * np.sqrt(t) + 0.5 * u ** 2 * sig ** 2 * t, kap) +
                      1j * u * eta_ * sig * np.sqrt(t))

    M = 15  # FFT grid size
    dz = 0.01  # FFT grid spacing
    a = 1 / 3

    def c_bac(x, t, eta_, kap, sig, B, Ra):
        return B * (Ra + (np.exp(-a * x) / (2 * np.pi)) *
                    fft_bachelier(lambda u: phi(u, sig, eta_, kap, t), M, dz, x, a))

    def i0(eta_, kap):
        return np.sqrt(2 * np.pi) * c_bac(0, 1, eta_, kap, 1, 1, 0)  # Normalizer

    # Price the two call legs (base η)
    c1 = c_bac((k1 - strike) / (sigma_atm[-1] * np.sqrt(ttm[-1])), 1, eta, kappa, 1 / i0(eta, kappa), 1, 0) * \
              sigma_atm[-1] * np.sqrt(ttm[-1]) * disc[-1]

    c2 = c_bac((k2 - strike) / (sigma_atm[-1] * np.sqrt(ttm[-1])), 1, eta, kappa, 1 / i0(eta, kappa), 1, 0) * \
              sigma_atm[-1] * np.sqrt(ttm[-1]) * disc[-1]

    bull_spread = c1 - c2  # Long K1, short K2

    # Price the two call legs (shocked η)
    c1s = c_bac((k1 - strike) / (sigma_atm[-1] * np.sqrt(ttm[-1])), 1, eta_shocked, kappa, 1 / i0(eta_shocked, kappa),
                1, 0) * sigma_atm[-1] * np.sqrt(ttm[-1]) * disc[-1]

    c2s = c_bac((k2 - strike) / (sigma_atm[-1] * np.sqrt(ttm[-1])), 1, eta_shocked, kappa, 1 / i0(eta_shocked, kappa),
                1, 0) * sigma_atm[-1] * np.sqrt(ttm[-1]) * disc[-1]

    bull_spread_shocked = c1s - c2s

    # ---------- 5) Sensitivity of the bull-spread to η --------------------
    delta_bull_eta = bull_spread_shocked - bull_spread  # € change for +1 bp

    # ---------- 6) Notional required to offset exotic η-risk --------------
    n_bull = -delta_price_eta / delta_bull_eta  # (+) buy, (-) sell

    # ---------- 7) Pack the output ----------------------------------------
    bull_info = np.vstack([n_bull, k1, k2, np.ceil(n_bull / (option_prices_new[2 * (len(ttm)-1), idx_below] +
                                                          option_prices_new[2 * (len(ttm)-1), idx_above]))])
    return bull_info

def hedge_kappa_strangle(strike, strikes, option_prices_new, sigma_atm, ttm, ttmlong, disc, eta, kappa, 
                          sigmat, alpha, fwd_prices, n_sim, increment, notional, price_bac):
    """
    Hedge the kappa (κ) sensitivity of an exotic position using a 10-wide
    Bachelier strangle (long call K+10, long put K-10).
    
    Args:
        strike: The strike price of the exotic option.
        strikes: Strike prices of the market options.
        option_prices_new: Market option prices for model calibration.
        sigma_atm: At-the-money volatility.
        ttm: Time to maturity vector.
        ttmlong: Long maturity time vector.
        disc: Discount factors.
        eta: The initial eta (volatility of volatility parameter).
        kappa: The initial kappa (volatility sensitivity parameter).
        sigmat: The model's volatility term structure.
        alpha: The model parameter for the volatility structure.
        fwd_prices: The forward prices.
        n_sim: Number of simulations for pricing.
        increment: The increment step size for simulations.
        notional: The notional amount for pricing.
        price_bac: The base price of the exotic option.

    Returns:
        strangle_info: A 4x1 numpy array with the following values:
            - Notional for the strangle to buy or sell
            - Lower-strike call (K − 10)
            - Upper-strike call (K + 10)
            - Number of strangle contracts
    """
    bp = 1e-4
    
    # 1) Shock kappa by +10 bp
    kappa_shocked = kappa + 10 * bp
    
    # 2) Re-price the exotic payoff
    price_shocked, _ = pricing_exotic(strike, disc, ttmlong, fwd_prices, n_sim, increment, eta, 
                                      kappa_shocked, sigmat, alpha, 0, 0, notional, 'Bac')
    print(price_shocked)
    print(price_bac)
    
    delta_price_kappa = price_shocked - price_bac  # Euro change for +10 bp
    
    # 3) Locate the call (K+10) and put (K-10) strikes
    dist_call = np.abs(strikes - strike - 10)
    dist_call[np.isnan(option_prices_new[-2, :])] = np.inf
    idx_call = np.argmin(dist_call)
    k_call = strikes[idx_call]  # Lower-leg strike
    print(f"Call strike: {k_call}")
    
    dist_put = np.abs(strikes - strike + 10)
    dist_put[np.isnan(option_prices_new[-1, :])] = np.inf
    idx_put = np.argmin(dist_put)
    k_put = strikes[idx_put]  # Upper-leg strike
    print(f"Put strike: {k_put}")
    
    # 4) Lewis-FFT Bachelier helpers
    def psi(u, kap):
        return (1 / kap) * (1 - alpha) / alpha * (1 - (1 + (u * kap) / (1 - alpha))**alpha)
    
    def phi(u, sig, eta_, kap, t):
        return np.exp(psi(1j * u * eta_ * sig * np.sqrt(t) + 0.5 * u**2 * sig**2 * t, kap) + 
                      1j * u * eta_ * sig * np.sqrt(t))
    
    M = 15  # FFT grid size
    dz = 0.01  # FFT grid spacing
    a = 1 / 3  # damping factor
    
    def c_bac(x, t, eta_, kap, sig, B, Ra):
        return B * (Ra + (np.exp(-a * x) / (2 * np.pi)) * fft_bachelier(lambda u: phi(u, sig, eta_, kap, t), M, dz, x, a))
    
    def i0(eta_, kap):
        return np.sqrt(2 * np.pi) * c_bac(0, 1, eta_, kap, 1, 1, 0)
    
    # 5) Price the strangle (base kappa)
    fwd_scale = sigma_atm[-1] * np.sqrt(ttm[-1])  # For normalized moneyness
    
    call = c_bac((k_call - strike) / fwd_scale, 1, eta, kappa, 1 / i0(eta, kappa), 1, 0) * fwd_scale * disc[-1]
    put = c_bac((k_put - strike) / fwd_scale, 1, eta, kappa, 1 / i0(eta, kappa), 1, 0) * fwd_scale * disc[-1] - disc[-1] * (strike - k_put)
    
    strangle = call + put
    print(f"Strangle price (base kappa): {strangle}")

    # 6) Price the strangle (shocked kappa)
    call_s = c_bac((k_call - strike) / fwd_scale, 1, eta, kappa_shocked, 1 / i0(eta, kappa_shocked), 1, 0) * fwd_scale * disc[-1]
    put_s = c_bac((k_put - strike) / fwd_scale, 1, eta, kappa_shocked, 1 / i0(eta, kappa_shocked), 1, 0) * fwd_scale * disc[-1] - disc[-1] * (strike - k_put)
    
    strangle_shocked = call_s + put_s
    print(f"Strangle price (shocked kappa): {strangle_shocked}")
    
    # 7) Sensitivity & required notional
    delta_strangle_kappa = strangle_shocked - strangle  # Euro change for +10 bp
    n_strangle = -delta_price_kappa / delta_strangle_kappa

    print(f"Delta strangle kappa: {delta_strangle_kappa}, Notional required: {n_strangle}")
    print(idx_call, idx_put)
    print(option_prices_new[2 * (len(ttm)-1), idx_call], option_prices_new[2 * (len(ttm)-1) + 1, idx_put])
    # 8) Pack and return
    strangle_info = np.array([n_strangle, k_call, k_put, np.ceil(n_strangle / (option_prices_new[2 * (len(ttm)-1), idx_call] + 
                                                                            option_prices_new[2 * (len(ttm)-1) + 1, idx_put]))])
    return strangle_info

def targetPL():
    return [20200603, 20200604, 20200605, 20200608, 20200609, 20200610,
            20200611, 20200612, 20200615, 20200616, 20200617]

def datesPL():
    return [datetime.strptime(date, "%d-%b-%Y") for date in [
        '03-Jun-2020', '04-Jun-2020', '05-Jun-2020', '08-Jun-2020', '09-Jun-2020',
        '10-Jun-2020', '11-Jun-2020', '12-Jun-2020', '15-Jun-2020', '16-Jun-2020', '17-Jun-2020'
    ]]

def hedgingQuality(callDir, putDir, target, t0, alpha, a, eta0, k0, Nsim, increment,
                   Notional, callVolatility, futuresInfoVol, bullInfo, numberFutureEta,
                   strangleInfo, nFutureKappa, ptfValue) -> Tuple[float, float]:

    # Step 1 - Define characteristic function components
    def psi_fun(u, kappa):
        return (1 / kappa) * ((1 - alpha) / alpha) * (1 - (1 + (u * kappa) / (1 - alpha)) ** alpha)

    
    def phi_fun(u, sigma, eta, kappa, t):
        z = 1j * u * eta * sigma * np.sqrt(t) + 0.5 * (u ** 2) * sigma ** 2 * t
        return np.exp(psi_fun(z, kappa) + 1j * u * eta * sigma * np.sqrt(t))
    
    # Step 2 - FFT parameters
    M = 15
    dz = 0.01

    def C_price_Lewis(x, t, eta, kappa, sigma_t, B, Ra):
        return B * (Ra + np.exp(-a * x) / (2 * np.pi) *
                    fft_bachelier(lambda u: phi_fun(u, sigma_t, eta, kappa, t), M, dz, x, a))

    # Step 3 - Market data and calibration
    optionPrices, strikes = build_option_prices(callDir, putDir, target)
    disc, fwdPrices = bootstrap(optionPrices, strikes)
    dates = getDatesNew(t0)
    ttm = [yearfrac(dates[0], d, 3) for d in dates[1:]]
    optionPricesNew = put_call_parity(optionPrices[2:, :], strikes, disc, fwdPrices)
    sigmaATM = ATMvols(optionPricesNew, strikes, fwdPrices, disc, ttm)
    eta, kappa, I0 = calibrate_add_bach(optionPricesNew, sigmaATM, strikes, fwdPrices, ttm,
                                      disc, alpha, a, eta0, k0, 1, len(ttm))
    sigmat = sigmaATM / I0

    # Step 4 - Price exotic
    ttmlong = np.concatenate(([0], ttm))
    strike = fwdPrices[-1]
    priceBac, _ = pricing_exotic(strike, disc, ttmlong, fwdPrices, Nsim, increment,
                                eta, kappa, sigmat, alpha, 0, 0, Notional, 'Bac')

    # Step 5 - Hedge instruments
    priceCallVol = np.zeros(callVolatility.shape[1])
    for i in range(len(priceCallVol)):
        t = int(callVolatility[3, i])
        x = (callVolatility[1, i] - fwdPrices[t]) / (sigmaATM[t] * np.sqrt(ttm[t]))
        priceCallVol[i] = C_price_Lewis(x, 1, eta, kappa, 1, 1, 0) * sigmaATM[t] * np.sqrt(ttm[t]) * disc[t]
        
        #if np.isnan(sigmaATM[t]) or sigmaATM[t] == 0 or ttm[t] == 0:
        #    print(f"Invalid inputs at t={t}: sigmaATM={sigmaATM[t]}, ttm={ttm[t]}")


    valueFuturesVol = fwdPrices[futuresInfoVol[2, :].astype(int)] * futuresInfoVol[3, :]
    
    # Bull spread value
    x1 = (bullInfo[1] - fwdPrices[-1]) / (sigmaATM[-1] * np.sqrt(ttm[-1]))
    x2 = (bullInfo[2] - fwdPrices[-1]) / (sigmaATM[-1] * np.sqrt(ttm[-1]))
    bullCall1 = C_price_Lewis(x1, 1, eta, kappa, 1, 1, 0) * sigmaATM[-1] * np.sqrt(ttm[-1]) * disc[-1]
    bullCall2 = C_price_Lewis(x2, 1, eta, kappa, 1, 1, 0) * sigmaATM[-1] * np.sqrt(ttm[-1]) * disc[-1]
    bullValue = bullCall1 - bullCall2

    # Strangle value
    x_call = (strangleInfo[1] - fwdPrices[-1]) / (sigmaATM[-1] * np.sqrt(ttm[-1]))
    x_put = (strangleInfo[2] - fwdPrices[-1]) / (sigmaATM[-1] * np.sqrt(ttm[-1]))
    call_strike_value = C_price_Lewis(x_call, 1, eta, kappa, 1, 1, 0) * sigmaATM[-1] * np.sqrt(ttm[-1]) * disc[-1]
    put_approx_value = C_price_Lewis(x_put, 1, eta, kappa, 1, 1, 0) * sigmaATM[-1] * np.sqrt(ttm[-1]) * disc[-1]
    put_call_parity_term = disc[-1] * (fwdPrices[-1] - strangleInfo[2])
    strangleValue = call_strike_value - (put_approx_value - put_call_parity_term)
    
    #print("---------- DEBUG -----------")
    #print("priceBac =", priceBac)
    #print("callVolVol =", np.sum(priceCallVol * callVolatility[4, :]))
    #print("futuresVol =", np.sum(valueFuturesVol))
    #print("bullValue * bullInfo[3] =", bullValue * bullInfo[3])
    #print("numberFutureEta * fwdPrices[-1] =", numberFutureEta * fwdPrices[-1])
    #print("strangleValue =", strangleValue)
    #print("nFutureKappa * fwdPrices[-1] =", nFutureKappa * fwdPrices[-1])
    #print("ptfValue =", ptfValue)

    # Step 6 - Total P&L
    P_L = (
        priceBac +
        np.sum(priceCallVol * callVolatility[4, :]) +
        np.sum(valueFuturesVol) +
        bullValue * bullInfo[3] +
        numberFutureEta * fwdPrices[-1] +
        strangleValue +
        nFutureKappa * fwdPrices[-1] -
        ptfValue
    )

    return P_L, priceBac


def evaluatePL(callDir, putDir, targets, datesPL, alpha, a, eta0, k0, Nsim, increment,
               Notional, callVolatility, futuresInfoVol, bullInfo, numberFutureEta,
               strangleInfo, numberFutureKappa, ptfValue, CostHedging, priceBac):
    
    PLNoHedging = np.zeros(len(targets))
    PLHedging = np.zeros(len(targets))
    priceBacNew = np.zeros(len(targets))

    for i in range(len(targets)):
        PLHedging[i], priceBacNew[i] = hedgingQuality(
            callDir, putDir, targets[i], datesPL[i], alpha, a, eta0, k0, Nsim, increment,
            Notional, callVolatility, futuresInfoVol, bullInfo, numberFutureEta,
            strangleInfo, numberFutureKappa, ptfValue
        )
        PLNoHedging[i] = priceBacNew[i] - priceBac
        print(f"Date: {datesPL[i]} | P&L: {PLHedging[i]} | Exotic Price: {priceBacNew[i]}")

    PLHedgingCost = PLHedging - CostHedging
    return PLNoHedging, PLHedging, PLHedgingCost


def plotPL(PLNoHedging, PLHedging, PLHedgingCost, datesPL):
    """
    Plots P&L over time for three strategies:
    - No hedging
    - Hedging
    - Hedging net of costs
    """
    # Sanity check
    if not (len(PLNoHedging) == len(PLHedging) == len(PLHedgingCost) == len(datesPL)):
        raise ValueError("All input vectors must have the same length.")

    # Plot
    plt.figure(figsize=(12, 5))
    plt.grid(True)
    
    # Plot lines
    p1, = plt.plot(datesPL, PLNoHedging, '-',  linewidth=2, label='No Hedging', color=(0.25, 0.47, 0.85))
    p2, = plt.plot(datesPL, PLHedging, '--', linewidth=2, label='Hedging', color=(0.90, 0.45, 0.13))
    p3, = plt.plot(datesPL, PLHedgingCost, ':', linewidth=2, label='Hedging (cost included)', color=(0.30, 0.70, 0.30))

    # Formatting
    plt.title('P/L evolution comparison', fontsize=14, fontweight='bold')
    plt.xlabel('Date', fontsize=12)
    plt.ylabel('P/L (€)', fontsize=12)
    plt.legend(loc='best', fontsize=10)
    plt.axhline(0, color='black', linewidth=1)

    # X-axis ticks formatting
    ax = plt.gca()
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%b-%Y'))
    plt.xticks(rotation=30)
    ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.show()

