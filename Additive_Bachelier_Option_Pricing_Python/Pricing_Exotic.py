import numpy as np
from scipy.stats import norm
from numpy.fft import fft
from scipy.interpolate import interp1d
import scipy.special
from scipy.special import gamma, betainc
from scipy.stats import norm

def pricing_exotic(strike, disc, ttm, fwd_prices, Nsim, increment, eta, kappa, sigmat, alpha, sigma_log, H, notional, flag):
    """
    Price exotic option using Monte Carlo with additive Bachelier or logistic models.
    """
    np.random.seed(0)  # For reproducibility

    cumulative = np.zeros((len(ttm)-1, len(increment)))
    U = np.random.rand(len(ttm)-1, Nsim)
    incr = np.zeros_like(U)
    sim_fwd = np.zeros((len(fwd_prices)-1, Nsim))
    payoff = np.zeros(Nsim)

    if flag == 'Bac':
        def psi_fun(u, kappa):
            return (1.0 / kappa) * (1 - alpha) / alpha * (1 - (1 + u * kappa / (1 - alpha)) ** alpha)

        def phi_add_bac(u, sigma, eta, kappa, t):
            return np.exp(psi_fun(1j * u * eta * sigma * np.sqrt(t) + 0.5 * u**2 * sigma**2 * t, kappa) + 1j * u * eta * sigma * np.sqrt(t))

        for i in range(len(ttm)-1):
            if i == 0:
                phi1 = lambda u: np.ones_like(u)
            else:
                phi1 = lambda u: phi_add_bac(u, sigmat[i-1], eta, kappa, ttm[i])
            phi2 = lambda u: phi_add_bac(u, sigmat[i], eta, kappa, ttm[i+1])
            phi_quot = lambda u: phi2(u) / phi1(u)
            cumulative[i, :] = lewis_formula(phi_quot, -0.01, 0, increment, ttm[i], ttm[i+1])

    elif flag == 'Log':
        def b(sigma, H, t):
            return (1 - np.exp(-t * sigma**(1 / H))) ** H

        def beta_complex(x, y):
            return gamma_complex(x) * gamma_complex(y) / gamma_complex(x + y)

        def phi_add_log(u, sigma, H, t):
            b_val = b(sigma_log, H, t)
            return (1 - b_val) * beta_complex(1 + (1j*u - 1)*b_val, 1 - 1j*u*b_val)

        for i in range(len(ttm)-1):
            if i == 0:
                phi1 = lambda u: np.ones_like(u)
            else:
                phi1 = lambda u: phi_add_log(u, sigma_log, H, ttm[i])
            phi2 = lambda u: phi_add_log(u, sigma_log, H, ttm[i+1])
            phi_quot = lambda u: phi2(u) / phi1(u)
            cumulative[i, :] = lewis_formula(phi_quot, 0.01, 1, increment, ttm[i], ttm[i+1])

    # Truncate cumulative matrix to [0,1]
    cumulative = np.clip(cumulative, 0, 1)
    for r in range(cumulative.shape[0]):
        cumulative[r, :] = np.maximum.accumulate(cumulative[r, :])  # enforce monotonicity

    # Inverse CDF method for increments
    for i in range(len(ttm)-1):
        x = cumulative[i, :]
        mask = np.concatenate(([True], np.diff(x) > 0))
        interp_func = interp1d(x[mask], increment[mask], kind='cubic', fill_value="extrapolate")
        incr[i, :] = interp_func(U[i, :])

    # Simulate forward prices
    for i in range(len(fwd_prices)-1):
        if i == 0:
            sim_fwd[i, :] = fwd_prices[-1] + incr[i, :]
        else:
            sim_fwd[i, :] = sim_fwd[i-1, :] + incr[i, :]

    sum_sim_fwd = np.mean(sim_fwd, axis=0)

    # Compute binary payoff
    payoff = (strike > sum_sim_fwd).astype(float)

    # Discounted expected payoff
    price = notional * disc[-1] * payoff
    stdev = np.std(price)
    ci_low = np.mean(price) - norm.ppf(0.995) * stdev / np.sqrt(Nsim)
    ci_high = np.mean(price) + norm.ppf(0.995) * stdev / np.sqrt(Nsim)
    CI = [ci_low, ci_high]
    price = np.mean(price)

    return price, CI


def lewis_formula(phi, a, Ra, x, s, t):
    """
    Lewis formula for computing cumulative distribution function using FFT.
    """
    M = 15
    int_fft = integral_fft(phi, M, x, a, s, t)
    cumulative = Ra - np.exp(-a*x) / (2 * np.pi) * int_fft
    return cumulative


def integral_fft(phi, M, query_points, a, s, t):
    """
    Integral calculation for Lewis formula using FFT.
    """
    N = 2 ** M
    xi_1 = -5 * np.sqrt(t - s)
    d_xi = -2 * xi_1 / (N - 1)
    dz = 2 * np.pi / (N * d_xi)
    z_1 = -dz * (N - 1) / 2

    xi = np.linspace(xi_1, -xi_1, N)
    z = np.linspace(z_1, -z_1, N)

    f = phi(xi - 1j * a) / (1j * xi + a)
    f_tilde = f * np.exp(-1j * z_1 * d_xi * np.arange(N))

    FFT = fft(f_tilde)

    prefactor = d_xi * np.exp(-1j * xi_1 * z)
    I = prefactor * FFT
    I = np.real(I)

    interp_func = interp1d(z, I, kind='linear', fill_value="extrapolate")
    I_interp = interp_func(query_points)

    return I_interp


def gamma_complex(z):
    """
    Gamma function for complex inputs using Lanczos approximation.
    """
    #return np.exp(loggamma_complex(z))
    return scipy.special.gamma(z)

def beta_complex(x, y):
    """
    Beta function for complex inputs.
    """
    return gamma_complex(x) * gamma_complex(y) / gamma_complex(x + y)

# Lanczos approximation for log gamma function for complex numbers
def loggamma_complex(z):
    coeff = np.array([
        0.9999999999998099322768470047348,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ])
    g = 7
    z = np.asarray(z, dtype=np.complex128)
    mask = (z.real < 0.5)
    lg = np.zeros(z.shape, dtype=np.complex128)

    # Reflection formula for Re(z)<0.5
    if np.any(mask):
        lg[mask] = np.log(np.pi) - np.log(np.sin(np.pi * z[mask])) - loggamma_complex(1 - z[mask])

    if np.any(~mask):
        z1 = z[~mask] - 1
        x = coeff[0] * np.ones_like(z1)
        for i in range(1, len(coeff)):
            x += coeff[i] / (z1 + i)
        t = z1 + g + 0.5
        lg[~mask] = 0.5*np.log(2*np.pi) + (z1 + 0.5)*np.log(t) - t + np.log(x)

    return lg

