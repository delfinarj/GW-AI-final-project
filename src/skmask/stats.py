"""Small statistical helpers shared by the analyses."""
from scipy.stats import beta


def clopper_pearson(k, n, level=0.95):
    """Exact binomial confidence interval for k successes in n trials."""
    lo = beta.ppf((1 - level) / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - (1 - level) / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)
