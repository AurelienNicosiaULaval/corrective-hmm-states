"""EM (Baum-Welch) for a univariate Gaussian HMM, with scaled forward-backward.

Design choices
--------------
- Scaled forward-backward (Rabiner-style) to avoid underflow.
- Multiple random initializations; the best log-likelihood is kept.
- A minimum standard deviation `sigma_min` prevents degenerate components.
- Fitted states are reordered by increasing mean before being returned.
- The EM trace is checked for (numerical) monotonicity of the log-likelihood.

Pure NumPy: no external HMM library is required.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SIGMA_MIN = 0.05  # minimal emission standard deviation


@dataclass
class HMMFit:
    """Result of one EM fit (best over random starts), states sorted by mean."""
    K: int
    loglik: float
    pi: np.ndarray          # (K,) initial distribution
    A: np.ndarray           # (K, K) transition matrix
    means: np.ndarray       # (K,)
    sds: np.ndarray         # (K,)
    gamma: np.ndarray       # (T, K) posterior state probabilities
    viterbi: np.ndarray     # (T,) most likely state path (0-based)
    n_iter: int
    loglik_trace: list = field(default_factory=list)

    def n_params(self) -> int:
        return (self.K - 1) + self.K * (self.K - 1) + 2 * self.K

    def bic(self, T: int) -> float:
        return -2.0 * self.loglik + self.n_params() * np.log(T)

    def icl(self, T: int) -> float:
        """ICL = BIC + 2 * posterior entropy."""
        g = np.clip(self.gamma, 1e-300, 1.0)
        entropy = -float(np.sum(self.gamma * np.log(g)))
        return self.bic(T) + 2.0 * entropy


def _emission_probs(y, means, sds):
    var = sds ** 2
    B = np.exp(-0.5 * (y[:, None] - means[None, :]) ** 2 / var[None, :])
    B /= np.sqrt(2.0 * np.pi * var[None, :])
    return np.maximum(B, 1e-300)


def _fb_core(B, pi, A):
    """Scaled forward-backward recursions on the emission matrix B (T, K)."""
    T, K = B.shape
    alpha = np.zeros((T, K))
    scale = np.zeros(T)
    alpha[0] = pi * B[0]
    scale[0] = alpha[0].sum()
    alpha[0] /= scale[0]
    for t in range(1, T):
        alpha[t] = (alpha[t - 1] @ A) * B[t]
        scale[t] = alpha[t].sum()
        alpha[t] /= scale[t]
    beta = np.zeros((T, K))
    beta[T - 1] = 1.0
    for t in range(T - 2, -1, -1):
        beta[t] = A @ (B[t + 1] * beta[t + 1])
        beta[t] /= scale[t + 1]
    gamma = alpha * beta
    gamma /= gamma.sum(axis=1, keepdims=True)
    obs = B[1:] * beta[1:]
    xi = alpha[:-1, :, None] * A[None, :, :] * obs[:, None, :]
    xi_sum = (xi / xi.sum(axis=(1, 2), keepdims=True)).sum(axis=0)
    loglik = np.sum(np.log(scale))
    return loglik, gamma, xi_sum


def forward_backward(y, pi, A, means, sds):
    """Scaled forward-backward. Returns (loglik, gamma, xi_sum)."""
    B = _emission_probs(y, means, sds)
    loglik, gamma, xi_sum = _fb_core(B, pi, A)
    return float(loglik), gamma, xi_sum


def viterbi_path(y, pi, A, means, sds):
    """Viterbi decoding in log-space."""
    T, K = len(y), len(pi)
    logB = np.log(_emission_probs(y, means, sds))
    logA = np.log(np.maximum(A, 1e-300))
    delta = np.zeros((T, K))
    psi = np.zeros((T, K), dtype=int)
    delta[0] = np.log(np.maximum(pi, 1e-300)) + logB[0]
    for t in range(1, T):
        cand = delta[t - 1][:, None] + logA
        psi[t] = cand.argmax(axis=0)
        delta[t] = cand.max(axis=0) + logB[t]
    path = np.zeros(T, dtype=int)
    path[-1] = int(delta[-1].argmax())
    for t in range(T - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]
    return path


def _em_single(y, K, rng, max_iter, tol, sigma_min):
    """One EM run from a random initialization. Returns dict or None on failure."""
    T = len(y)
    y_sd = float(np.std(y))
    # Random-but-sensible initialization: spread means over quantiles, jittered.
    means = np.quantile(y, np.linspace(0.1, 0.9, K)) + rng.normal(0.0, 0.25 * y_sd, K)
    sds = np.full(K, max(2 * sigma_min, y_sd / 2.0)) * rng.uniform(0.6, 1.4, K)
    pi = rng.dirichlet(np.ones(K))
    A = np.vstack([rng.dirichlet(np.ones(K) + 6.0 * np.eye(K)[i]) for i in range(K)])

    trace = []
    prev = -np.inf
    for it in range(max_iter):
        loglik, gamma, xi_sum = forward_backward(y, pi, A, means, sds)
        trace.append(loglik)
        # Monotonicity check (tiny numerical tolerance allowed)
        if loglik < prev - 1e-6 * (1.0 + abs(prev)):
            raise RuntimeError(
                f"EM log-likelihood decreased at iteration {it}: {prev} -> {loglik}"
            )
        converged = (loglik - prev) < tol * (1.0 + abs(prev))
        prev = loglik
        # M-step
        weights = gamma.sum(axis=0) + 1e-12
        pi = gamma[0] + 1e-10
        pi /= pi.sum()
        A = xi_sum + 1e-10
        A /= A.sum(axis=1, keepdims=True)
        means = (gamma * y[:, None]).sum(axis=0) / weights
        var = (gamma * (y[:, None] - means[None, :]) ** 2).sum(axis=0) / weights
        sds = np.sqrt(np.maximum(var, sigma_min ** 2))
        if converged and it > 0:
            break
    loglik, gamma, _ = forward_backward(y, pi, A, means, sds)
    trace.append(loglik)
    return dict(loglik=loglik, pi=pi, A=A, means=means, sds=sds,
                gamma=gamma, n_iter=len(trace) - 1, trace=trace)


def fit_hmm(y, K, seed, n_starts=20, max_iter=300, tol=1e-6,
            sigma_min=SIGMA_MIN) -> HMMFit:
    """Fit a K-state univariate Gaussian HMM by EM with multiple random starts.

    States in the returned fit are ordered by increasing mean.
    """
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_starts):
        res = _em_single(y, K, rng, max_iter, tol, sigma_min)
        if best is None or res["loglik"] > best["loglik"]:
            best = res
    order = np.argsort(best["means"])
    pi = best["pi"][order]
    A = best["A"][np.ix_(order, order)]
    means = best["means"][order]
    sds = best["sds"][order]
    gamma = best["gamma"][:, order]
    vit = viterbi_path(y, pi, A, means, sds)
    return HMMFit(K=K, loglik=best["loglik"], pi=pi, A=A, means=means,
                  sds=sds, gamma=gamma, viterbi=vit, n_iter=best["n_iter"],
                  loglik_trace=best["trace"])


def validate_fit(fit: HMMFit, T: int) -> None:
    """Assertions on a fitted model; raises AssertionError on failure."""
    assert fit.gamma.shape == (T, fit.K), "gamma has wrong shape"
    assert np.allclose(fit.A.sum(axis=1), 1.0, atol=1e-8), "rows of A must sum to 1"
    assert np.allclose(fit.pi.sum(), 1.0, atol=1e-8), "pi must sum to 1"
    assert np.all(fit.sds > 0), "estimated sds must be positive"
    assert np.allclose(fit.gamma.sum(axis=1), 1.0, atol=1e-8), \
        "posterior probabilities must sum to 1 at each time"
    assert np.all(np.isfinite(fit.gamma)), "gamma contains NaN/Inf"
    assert np.all(np.isfinite([fit.loglik, *fit.means, *fit.sds])), \
        "parameters contain NaN/Inf"
    tr = np.array(fit.loglik_trace)
    assert np.all(np.diff(tr) >= -1e-6 * (1.0 + np.abs(tr[:-1]))), \
        "EM log-likelihood is not non-decreasing (beyond tolerance)"
