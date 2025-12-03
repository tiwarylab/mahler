from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import curve_fit
from scipy.stats import expon, kstest

__all__ = [
    "cdf_from_sample",
    "cdf_from_param",
    "fit_poisson",
    "bootstrap",
]


def optimal_bin(t: ArrayLike, nbins: int = 1000) -> NDArray[np.float64]:
    """Return logarithmically spaced bins that cover the data with padding."""
    t_array = np.asarray(t, dtype=float)
    log_t = np.log10(t_array)

    min_logt = float(log_t.min())
    max_logt = float(log_t.max())

    spread_logt = max_logt - min_logt
    shift = np.maximum(spread_logt * 0.1, 1.0)

    bin_lo = min_logt - shift
    bin_hi = max_logt + shift

    return np.logspace(bin_lo, bin_hi, nbins)


def cdf_from_sample(
    t: ArrayLike,
    bins: Sequence[float] | NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute an empirical CDF over the provided bins."""
    if bins is None:
        bins = optimal_bin(t)

    hist, edges = np.histogram(t, bins=bins, density=False)
    cum = np.cumsum(hist) / np.asarray(t).shape[0]
    xx = (edges[1:] + edges[:-1]) / 2

    return np.array([xx, cum], dtype=float)


def cdf_from_param(x: ArrayLike, inv_lambda: float) -> NDArray[np.float64]:
    """Exponential CDF parameterized by the inverse rate."""
    return expon.cdf(x, scale=inv_lambda)


def fit_poisson(t: ArrayLike) -> tuple[float, float, float]:
    """Fit a Poisson-process CDF (exponential waiting time) to empirical data."""
    _cdf = cdf_from_sample(t)

    # Use the sample median as a robust initial guess.
    initial_guess = float(np.median(t))

    args, pcov = curve_fit(cdf_from_param, *_cdf, p0=initial_guess)

    inv_lambda = float(args[0])
    std_args = float(np.sqrt(np.diag(pcov))[0])

    p_val = float(kstest(t, expon.cdf, (0, inv_lambda)).pvalue)

    return inv_lambda, std_args, p_val


def bootstrap(
    data: ArrayLike,
    n: int = 300,
    n_samples: int = 1000,
    random_seed: int = 42,
) -> NDArray[np.float64]:
    """
    Bootstrap resampling of data.

    :param data: array-like
    :param n: number of samples per resampling
    :param n_samples: number of resamplings
    :return: list of resampled arrays
    """
    rng = np.random.default_rng(seed=random_seed)
    results: list[float] = []
    for _ in range(n_samples):
        res, _, _ = fit_poisson(rng.choice(data, size=n, replace=True))
        results.append(res)
    return np.asarray(results, dtype=float)
