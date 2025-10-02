import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt


def get_dates():
    """Returns futures expiration dates as datetime objects."""
    return pd.to_datetime([
        '2020-06-02', '2020-08-17', '2020-11-17', '2021-02-17',
        '2021-05-17', '2021-08-17', '2021-11-16', '2022-05-17',
        '2022-11-16'
    ])

def read_csv_data(file1, file2, target):
    """Reads and filters data for a given target date, returns prices as numpy array."""
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    filtered1 = df1[df1.iloc[:, 0] == target].iloc[:, 1:]
    filtered2 = df2[df2.iloc[:, 0] == target].iloc[:, 1:]

    return np.vstack([filtered1.values, filtered2.values])

def build_option_prices(call_folder, put_folder, target):
    """Builds option price matrix and strike vector from calls and puts directories."""
    call_folder = Path(call_folder)
    put_folder = Path(put_folder)

    call_files = sorted(call_folder.glob('*call.csv'))
    if not call_files:
        raise FileNotFoundError(f'No "*call.csv" files found in {call_folder}')

    # Use first file to get strikes
    first_df = pd.read_csv(call_files[0])
    strikes = first_df.columns[1:].astype(float).values
    N = len(strikes)
    M = len(call_files)

    option_prices = np.zeros((2 * M, N))

    for k, call_file in enumerate(call_files):
        put_file = put_folder / call_file.name.replace('call', 'put')
        if not put_file.exists():
            raise FileNotFoundError(f'Missing put file: {put_file}')

        block = read_csv_data(call_file, put_file, target)
        option_prices[2 * k: 2 * k + 2] = block

    return option_prices, strikes

def bootstrap(option_prices, strikes):
    """
    Bootstrap discount factors and forward prices from option prices.

    Parameters:
        option_prices (ndarray): Matrix of option prices (2M x N)
        strikes (ndarray): Vector of strike prices (length N)

    Returns:
        disc (ndarray): Vector of discount factors (length M)
        fwd_prices (ndarray): Vector of forward prices (length M-1)
    """
    M = option_prices.shape[0] // 2
    N = option_prices.shape[1]

    disc = np.ones(M)
    fwd_prices = np.zeros(M - 1)

    for i in range(1, M):  # Start from the second maturity
        call_row = option_prices[2 * i]
        put_row = option_prices[2 * i + 1]

        call_put_diff, K = zip(*[
        (c - p, k) for c, p, k in zip(call_row, put_row, strikes)
        if np.isfinite(c) and np.isfinite(p)
        ])
        call_put_diff = np.array(call_put_diff)
        K = np.array(K)

        G_mean = np.mean(call_put_diff)
        K_mean = np.mean(K)

        # Calculate discount factor as negative of regression slope
        numerator = np.dot(call_put_diff - G_mean, K - K_mean)
        denominator = np.dot(K - K_mean, K - K_mean)
        disc[i] = - numerator / denominator

        # Linear regression: call_put_diff = a*K + b
        A = np.vstack([K, np.ones_like(K)]).T
        a, b = np.linalg.lstsq(A, call_put_diff, rcond=None)[0]

        # Forward price is intercept / discount
        fwd_prices[i - 1] = b / disc[i]

    return disc, fwd_prices

def yearfrac(t0, dates, basis=3):
    """
    Compute year fraction between t0 and each date in dates.

    Parameters:
        t0 (datetime): Starting date
        dates (Series or array): Target dates
        basis (int): Day count convention. 3 = ACT/365

    Returns:
        Array of year fractions
    """
    if basis == 3:
        return (dates - t0).days / 365.0
    else:
        raise NotImplementedError("Only ACT/365 (basis=3) is implemented")

def zero_rates(dates, discounts, t0):
    """
    Compute zero rates from discount factors and dates.

    Parameters:
        dates (Series): Dates corresponding to discount factors
        discounts (ndarray): Discount factors
        t0 (datetime): Reference (value) date

    Returns:
        Array of zero rates
    """
    deltas = yearfrac(t0, dates)
    return -np.log(discounts) / deltas

from datetime import datetime

def getDatesNew(t0):
    """
    Parameters:
    - t0: datetime.datetime, value date

    Returns:
    - list of datetime.datetime: vector of futures expiry dates
    """
    return [
        t0,
        datetime(2020, 8, 17),
        datetime(2020, 11, 17),
        datetime(2021, 2, 17),
        datetime(2021, 5, 17),
        datetime(2021, 8, 17),
        datetime(2021, 11, 16),
        datetime(2022, 5, 17),
        datetime(2022, 11, 16)
    ]

