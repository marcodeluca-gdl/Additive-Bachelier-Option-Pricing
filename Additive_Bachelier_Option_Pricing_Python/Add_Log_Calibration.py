import numpy as np
from scipy.optimize import minimize

# Additive Logistic Call Price Function
def call_price_add_log(s, S0, K):
    return s * np.log(1 + np.exp((S0 - K) / s))

# Additive Logistic Put Price Function
def put_price_add_log(s, S0, K):
    return s * np.log(1 + np.exp((K - S0) / s))

def calibrate_add_log(strikes, fwd_prices, ttm, option_prices):
    prices_model = []
    prices_mkt = []

    for i in range(1, len(ttm)):  # Skip the first maturity
        mon = strikes - fwd_prices[i]
        valid_call = np.where(np.isfinite(option_prices[2*i]) & (mon > 0) & (mon < 30))[0]
        valid_put = np.where(np.isfinite(option_prices[2*i+1]) & (mon < 0) & (mon > -30))[0]

        for j in valid_call:
            prices_mkt.append(option_prices[2*i, j])
            prices_model.append(lambda p, i=i, j=j: call_price_add_log(p[0] * ttm[i]**p[1], fwd_prices[i], strikes[j]))

        for j in valid_put:
            prices_mkt.append(option_prices[2*i+1, j])
            prices_model.append(lambda p, i=i, j=j: put_price_add_log(p[0] * ttm[i]**p[1], fwd_prices[i], strikes[j]))

    # Vectorized model function
    def prices_model_vec(p):
        return np.array([f(p) for f in prices_model])

    # Objective function (mean absolute percentage error)
    def objective(p):
        model_vals = prices_model_vec(p)
        return np.sum((model_vals - prices_mkt)**2)

    # Bounds for sigma and H
    bounds = [(0, 20), (0, 1)]  # sigma in [0, 20], H in [0, 1]
    p0 = [5, 0.5]  # Initial guess for sigma and H
    
    options = {
        'disp': False,
        'maxiter': 5000,  # Increase maximum iterations
        'xtol': 1e-4,     # Relax step size tolerance
        'gtol': 1e-4,     # Relax gradient tolerance
        'barrier_tol': 1e-8
    }

    # Optimization using 'trust-constr' method
    result = minimize(objective, p0, bounds=bounds, method='trust-constr', options=options)

    # Extract optimized parameters
    sigma, H = result.x
    return sigma, H

# Function - Calibration of sigma with H = 0.5 (constant H)
def calibrate_add_log_h_constant(strikes, fwd_prices, ttm, option_prices):
    prices_model = []
    prices_mkt = []

    # Loop through each maturity (skip the first maturity)
    for i in range(1, len(ttm)):
        mon = strikes - fwd_prices[i]
        valid_call = np.where(~np.isnan(option_prices[2*i]) & (mon > 0) & (mon < 30))[0]
        valid_put = np.where(~np.isnan(option_prices[2*i+1]) & (mon < 0) & (mon > -30))[0]

        # Add valid call prices and model functions
        for j in valid_call:
            prices_mkt.append(option_prices[2*i, j])
            prices_model.append(lambda p, i=i, j=j: call_price_add_log(p * (ttm[i]**0.5), fwd_prices[i], strikes[j]))

        # Add valid put prices and model functions
        for j in valid_put:
            prices_mkt.append(option_prices[2*i+1, j])
            prices_model.append(lambda p, i=i, j=j: put_price_add_log(p * (ttm[i]**0.5), fwd_prices[i], strikes[j]))

    # Vectorized model function
    def prices_model_vec(p):
        return np.array([f(p) for f in prices_model])

    # Objective function (mean absolute percentage error)
    def objective(p):
        model_vals = prices_model_vec(p)
        return np.sum((model_vals - prices_mkt)**2)

    # Bounds for sigma
    bounds = [(0, 20)]
    p0 = [5]  # Initial guess for sigma

    # Optimization using 'trust-constr' method with relaxed tolerances
    options = {
        'disp': True,
        'maxiter': 10000  # Increase maximum iterations
        #'xtol': 1e-16,      # Relax step size tolerance
        #'gtol': 1e-16,      # Relax gradient tolerance
        #'barrier_tol': 1e-16
    }

    result = minimize(objective, p0, bounds=bounds, method='Nelder-Mead', options=options)

    # Extract optimized parameter
    sigma = result.x[0]
    return sigma
