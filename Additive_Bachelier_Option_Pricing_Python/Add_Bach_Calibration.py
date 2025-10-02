import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds, fsolve, brentq
from numpy.fft import fft
from functools import partial


def put_call_parity(option_prices, strikes, discounts, forwards):
    """
    Fill missing call/put prices using put-call parity.
    """
    option_prices = option_prices.copy()
    M, N = option_prices.shape[0] // 2, option_prices.shape[1]

    for i in range(M):
        DF = discounts[i + 1]
        F0 = forwards[i]

        call_row, put_row = zip(*[
            (
                P + DF * (F0 - K) if np.isnan(C) and np.isfinite(P) and (P + DF * (F0 - K)) > 0 else C,
                C - DF * (F0 - K) if np.isnan(P) and np.isfinite(C) and (C - DF * (F0 - K)) > 0 else P
            )
            for C, P, K in zip(option_prices[2 * i], option_prices[2 * i + 1], strikes)
        ])

        option_prices[2 * i] = np.array(call_row)
        option_prices[2 * i + 1] = np.array(put_row)

    return option_prices

def ATMvols(option_prices, strikes, fwd_prices, disc, ttm):
    """
    Compute ATM implied volatilities using Bachelier model.
    """
    n = len(ttm)
    sigma_atm = np.full(n, np.nan)  # inizializza tutto a NaN

    for i in range(n):
        call_prices = option_prices[2 * i]
        valid = np.isfinite(call_prices)

        if np.sum(valid) >= 2:
            K_valid = strikes[valid]
            C_valid = call_prices[valid]

            interp_func = interp1d(K_valid, C_valid, kind='cubic', fill_value="extrapolate", bounds_error=False)
            call_mkt = interp_func(fwd_prices[i])

            def objective(sigma):
                return call_price_bachelier(disc[i + 1], ttm[i], fwd_prices[i], fwd_prices[i], sigma) - call_mkt

            sigma_atm[i] = fsolve(objective, x0=0.2)[0]

    return sigma_atm


def market_vols(option_prices, strikes, fwd_prices, disc, ttm):
    """
    Compute the implied volatility surface using the Bachelier model.
    
    Parameters:
        option_prices (2D np.array): matrix of option prices. Every first row is used for calls.
        strikes (1D np.array): strike prices.
        fwd_prices (1D np.array): forward prices for each maturity.
        disc (1D np.array): discount factors (length should be len(ttm) + 1).
        ttm (1D np.array): times to maturity.
    
    Returns:
        implied_vols (2D np.array): matrix of implied volatilities with same shape as input prices.
    """
    n_maturities = len(ttm)
    n_strikes = len(strikes)
    implied_vols = np.full((n_maturities, n_strikes), np.nan)

    for i in range(n_maturities):
        for j in range(n_strikes):
            price_call = option_prices[2 * i, j]

            # Skip invalid or missing prices
            if not np.isfinite(price_call):
                continue

            # Define the objective function: Bachelier call price minus observed price
            def objective(sigma):
                return call_price_bachelier(disc[i + 1], ttm[i], strikes[j], fwd_prices[i], sigma) - price_call

            # Attempt 1: fsolve with a good initial guess
            try:
                if i < 2:
                    # Heuristic initialization for early maturities
                    x0 = price_call / (disc[i + 1] * np.sqrt(ttm[i]))
                else:
                    x0 = 0.2  # fallback for longer maturities

                vol, info, ier, _ = fsolve(objective, x0, full_output=True)

                if ier == 1 and vol[0] >= 0:
                    implied_vols[i, j] = vol[0]
                    continue  # successful solution, skip to next
            except:
                pass  # silently ignore and move to brentq

            # Attempt 2: brentq method on a predefined interval
            try:
                implied_vols[i, j] = brentq(objective, a=1e-6, b=5.0)
            except:
                pass  # if it fails, keep NaN

    return implied_vols


def call_price_bachelier(B, ttm, K, F0, sigma):
    if sigma == 0 or ttm == 0:
        return max(B * (F0 - K), 0)

    x = K - F0
    y = x / (sigma * np.sqrt(ttm))
    cb = -y * norm.cdf(-y) + norm.pdf(-y)
    return B * sigma * np.sqrt(ttm) * cb


def plot_vol_smiles(strikes, fwd_prices, implied_vols, ttm):
    for i in range(1, len(fwd_prices)):
        x = strikes - fwd_prices[i]
        y = implied_vols[i]
        valid = [j for j, val in enumerate(y) if np.isfinite(val)]

        if len(valid) < 3:
            continue

        xv = x[valid]
        yv = y[valid]
        x_fine = np.linspace(min(xv), max(xv), 400)
        y_fine = interp1d(xv, yv, kind='cubic')(x_fine)

        plt.figure()
        plt.plot(x_fine, y_fine, 'r-', linewidth=1.5)
        plt.grid(True)
        plt.xlim([-30, 30])
        plt.xlabel('Moneyness (K - F)')
        plt.ylabel('Implied Volatility')
        plt.title(f'Smile at TTM = {ttm[i]:.2f} yrs')
        plt.show()


def calibrate_add_bach(option_prices, sigma_atm, strikes, fwd_prices, ttm, disc, alpha, a, eta0, k0, start, last):
    def psi_fun(u, kappa):
        return (1/kappa) * (1-alpha)/alpha * (1 - (1 + u*kappa/(1-alpha))**alpha)

    def phi_fun(u, sigma, eta, kappa, t):
        return np.exp(psi_fun(1j*u*eta*sigma*np.sqrt(t) + 0.5*u**2*sigma**2*t, kappa) + 1j*u*eta*sigma*np.sqrt(t))

    def c_price_lewis(x, t, eta, kappa, sigma_t, B, Ra):
        phi = lambda u: phi_fun(u, sigma_t, eta, kappa, t)
        result = fft_bachelier(phi, 15, 0.01, x, a)
        if np.any(np.isnan(result)):
            print(f"fft_bachelier returned NaN for eta={eta}, kappa={kappa}, x={x}")
            result = np.nan_to_num(result, nan=0.0)
        return B * (Ra + np.exp(-a*x)/(2*np.pi) * result)

    def I0_fn(eta, kappa):
        val = np.sqrt(2*np.pi) * c_price_lewis(0, 1, eta, kappa, 1, 1, 0)
        if not np.isfinite(val):
            print(f"I0_fn returned NaN for eta={eta}, kappa={kappa}")
        return val

    def model_call_fn(x, i0_fn, c_price_fn):
        return lambda eta, kappa: c_price_fn(x, 1, eta, kappa, 1/i0_fn(eta, kappa), 1, 0)

    def model_put_fn(x, fwd, K, sig, t, i0_fn, c_price_fn):
        denom = sig * np.sqrt(t)
        denom = np.where(denom < 1e-6, 1e-6, denom)
        return lambda eta, kappa: c_price_fn(x, 1, eta, kappa, 1/i0_fn(eta, kappa), 1, 0) - (fwd - K) / denom

    prices_model = []
    prices_mkt = []

    for i in range(start, last):
        call_prices = option_prices[2*i, :]
        put_prices = option_prices[2*i+1, :]

        mon = strikes - fwd_prices[i]
        mask_call = np.isfinite(call_prices) & (mon > 0) & (mon < 30)
        mask_put = np.isfinite(put_prices) & (mon < 0) & (mon > -30)

        x_call = (strikes[mask_call] - fwd_prices[i]) / (sigma_atm[i] * np.sqrt(ttm[i]))
        x_put = (strikes[mask_put] - fwd_prices[i]) / (sigma_atm[i] * np.sqrt(ttm[i]))
        disc_i = disc[i + 1]

        prices_model.append(model_call_fn(x_call, I0_fn, c_price_lewis))
        prices_model.append(model_put_fn(x_put, fwd_prices[i], strikes[mask_put], sigma_atm[i], ttm[i], I0_fn, c_price_lewis))

        call_prices_norm = call_prices[mask_call] / (sigma_atm[i]*np.sqrt(ttm[i])*disc_i)
        put_prices_norm = put_prices[mask_put] / (sigma_atm[i]*np.sqrt(ttm[i])*disc_i)

        if np.any(np.isnan(call_prices_norm)) or np.any(np.isnan(put_prices_norm)):
            print(f"NaN in prices_mkt at maturity {i}")

        prices_mkt.extend(call_prices_norm)
        prices_mkt.extend(put_prices_norm)

    def prices_model_vec(p):
        vals = [f(p[0], p[1]) for f in prices_model]
        if any([np.any(np.isnan(v)) for v in vals]):
            print(f"NaN in prices_model_vec at p = {p}")
        return np.concatenate(vals)

    def objective(p):
        return np.sum((prices_model_vec(p) - prices_mkt)**2)

    bounds = Bounds([-np.inf, 1e-6], [np.inf, 500])
    options = {
        'xtol': 1e-8,
        'gtol': 1e-8,
        'barrier_tol': 1e-8,
        'maxiter': 5000,
        'verbose': 0
    }

    res = minimize(
        objective,
        x0=[eta0, k0],
        method='trust-constr',
        bounds=bounds,
        options=options
    )

    eta, kappa = res.x
    return eta, kappa, I0_fn(eta, kappa)



def fft_bachelier(phi, M, dz, query_points, a):
    N = 2**M
    z_1 = -(N-1)/2 * dz
    z = np.linspace(z_1, -z_1, N)
    d_xi = 2 * np.pi / (N * dz)
    xi_1 = -(N-1)/2 * d_xi
    xi = np.linspace(xi_1, -xi_1, N)

    f = phi(xi - 1j*a) / (1j*xi + a)**2
    f_tilde = f * np.exp(-1j * z_1 * d_xi * np.arange(N))

    FFT = fft(f_tilde)
    prefactor = d_xi * np.exp(-1j * xi_1 * z)
    I = prefactor * FFT
    I = np.real(I)

    interp = interp1d(z, I, kind='cubic', fill_value='extrapolate')
    return interp(query_points)

def plot_bachelier_parameters(eta, kappa, eta_t, kappa_t, ttm):
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(ttm[1:], eta_t, 'gs--', linewidth=1.2, markersize=5)
    plt.axhline(eta, color='b', linewidth=1.2)
    plt.xlabel(r'$t$ [years]')
    plt.ylabel(r'$\eta$')
    plt.legend([r'$\eta_t$', r'$\eta$'], loc='upper right')
    plt.grid(True)

    plt.subplot(2,1,2)
    plt.plot(ttm[1:], kappa_t, 'gs--', linewidth=1.2, markersize=5)
    plt.axhline(kappa, color='b', linewidth=1.2)
    plt.xlabel(r'$t$ [years]')
    plt.ylabel(r'$k$')
    plt.legend([r'$k_t$', r'$k$'], loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
