"""EM for univariate HMMs with Gaussian-mixture emissions.

The hidden state process is a standard finite-state Markov chain. Conditional
on hidden state k, the observation density is a finite mixture of M Gaussian
components. This is used as an enriched-emission diagnostic: if a corrective
extra Gaussian HMM state disappears after allowing mixtures within aggregate
states, the extra state is more naturally interpreted as an emission correction.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .hmm_em import SIGMA_MIN, _fb_core


@dataclass
class MixtureHMMFit:
    """Best EM fit over random starts, states sorted by emission mean."""
    K: int
    M: int
    loglik: float
    pi: np.ndarray              # (K,)
    A: np.ndarray               # (K, K)
    weights: np.ndarray         # (K, M)
    means: np.ndarray           # (K, M)
    sds: np.ndarray             # (K, M)
    gamma: np.ndarray           # (T, K)
    component_resp: np.ndarray  # (T, K, M)
    viterbi: np.ndarray         # (T,)
    n_iter: int
    loglik_trace: list = field(default_factory=list)

    def state_means(self) -> np.ndarray:
        return np.sum(self.weights * self.means, axis=1)

    def state_vars(self) -> np.ndarray:
        second = np.sum(self.weights * (self.sds ** 2 + self.means ** 2), axis=1)
        return second - self.state_means() ** 2

    def n_params(self) -> int:
        transition = (self.K - 1) + self.K * (self.K - 1)
        emissions = self.K * ((self.M - 1) + 2 * self.M)
        return transition + emissions

    def bic(self, T: int) -> float:
        return -2.0 * self.loglik + self.n_params() * np.log(T)


def _norm_pdf(y: np.ndarray, means: np.ndarray, sds: np.ndarray) -> np.ndarray:
    var = sds ** 2
    out = np.exp(-0.5 * (y[:, None, None] - means[None, :, :]) ** 2 / var[None, :, :])
    out /= np.sqrt(2.0 * np.pi * var[None, :, :])
    return np.maximum(out, 1e-300)


def _emission_probs(y, weights, means, sds):
    comp = _norm_pdf(y, means, sds)
    B = np.sum(weights[None, :, :] * comp, axis=2)
    return np.maximum(B, 1e-300), comp


def forward_backward_mixture(y, pi, A, weights, means, sds):
    """Scaled forward-backward for mixture-emission HMMs."""
    B, comp = _emission_probs(y, weights, means, sds)
    loglik, gamma, xi_sum = _fb_core(B, pi, A)
    numer = gamma[:, :, None] * weights[None, :, :] * comp
    component_resp = numer / np.maximum(B[:, :, None], 1e-300)
    return float(loglik), gamma, xi_sum, component_resp


def viterbi_path_mixture(y, pi, A, weights, means, sds):
    """Viterbi decoding with state-level mixture emission densities."""
    T, K = len(y), len(pi)
    B, _ = _emission_probs(y, weights, means, sds)
    logB = np.log(B)
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


def _initial_components(y, K, M, rng, sigma_min):
    y_sd = max(float(np.std(y)), sigma_min)
    q_state = np.linspace(0.15, 0.85, K)
    state_centers = np.quantile(y, q_state) + rng.normal(0.0, 0.15 * y_sd, K)
    offsets = np.linspace(-0.35, 0.35, M) * y_sd
    means = state_centers[:, None] + offsets[None, :] + rng.normal(0.0, 0.08 * y_sd, (K, M))
    sds = np.full((K, M), max(2 * sigma_min, 0.45 * y_sd)) * rng.uniform(0.7, 1.3, (K, M))
    weights = np.vstack([rng.dirichlet(np.ones(M) + 1.0) for _ in range(K)])
    pi = rng.dirichlet(np.ones(K))
    A = np.vstack([rng.dirichlet(np.ones(K) + 6.0 * np.eye(K)[i]) for i in range(K)])
    return pi, A, weights, means, sds


def _sort_components(weights, means, sds):
    weights = weights.copy()
    means = means.copy()
    sds = sds.copy()
    for k in range(weights.shape[0]):
        order = np.argsort(means[k])
        weights[k] = weights[k, order]
        means[k] = means[k, order]
        sds[k] = sds[k, order]
    return weights, means, sds


def _em_loop(y, pi, A, weights, means, sds, max_iter, tol, sigma_min):
    """Run mixture-emission EM from supplied parameters."""
    trace = []
    prev = -np.inf
    for it in range(max_iter):
        loglik, gamma, xi_sum, comp_resp = forward_backward_mixture(
            y, pi, A, weights, means, sds)
        trace.append(loglik)
        if loglik < prev - 1e-6 * (1.0 + abs(prev)):
            raise RuntimeError(
                f"Mixture EM log-likelihood decreased at iteration {it}: {prev} -> {loglik}"
            )
        converged = (loglik - prev) < tol * (1.0 + abs(prev))
        prev = loglik

        pi = gamma[0] + 1e-10
        pi /= pi.sum()
        A = xi_sum + 1e-10
        A /= A.sum(axis=1, keepdims=True)

        state_weights = gamma.sum(axis=0) + 1e-12
        comp_counts = comp_resp.sum(axis=0) + 1e-12
        weights = comp_counts / state_weights[:, None]
        weights /= weights.sum(axis=1, keepdims=True)
        means = np.sum(comp_resp * y[:, None, None], axis=0) / comp_counts
        var = np.sum(comp_resp * (y[:, None, None] - means[None, :, :]) ** 2, axis=0)
        var /= comp_counts
        sds = np.sqrt(np.maximum(var, sigma_min ** 2))
        weights, means, sds = _sort_components(weights, means, sds)
        if converged and it > 0:
            break

    loglik, gamma, _, comp_resp = forward_backward_mixture(y, pi, A, weights, means, sds)
    trace.append(loglik)
    return dict(loglik=loglik, pi=pi, A=A, weights=weights, means=means, sds=sds,
                gamma=gamma, component_resp=comp_resp, n_iter=len(trace) - 1,
                trace=trace)


def _em_single(y, K, M, rng, max_iter, tol, sigma_min):
    pi, A, weights, means, sds = _initial_components(y, K, M, rng, sigma_min)
    return _em_loop(y, pi, A, weights, means, sds, max_iter, tol, sigma_min)


def fit_mixture_hmm(y, K, M=2, seed=1, n_starts=12, max_iter=250, tol=1e-6,
                    sigma_min=SIGMA_MIN, screen_iter=None, refine_top=None,
                    screen_tol=1e-4) -> MixtureHMMFit:
    """Fit a K-state HMM with M Gaussian components per state."""
    rng = np.random.default_rng(seed)
    best = None
    if screen_iter is not None and refine_top is not None and refine_top < n_starts:
        screened = [
            _em_single(y, K, M, rng, int(screen_iter), screen_tol, sigma_min)
            for _ in range(n_starts)
        ]
        selected = sorted(screened, key=lambda res: res["loglik"], reverse=True)[
            :int(refine_top)
        ]
        for res in selected:
            refined = _em_loop(y, res["pi"], res["A"], res["weights"],
                               res["means"], res["sds"], max_iter, tol,
                               sigma_min)
            if best is None or refined["loglik"] > best["loglik"]:
                best = refined
    else:
        for _ in range(n_starts):
            res = _em_single(y, K, M, rng, max_iter, tol, sigma_min)
            if best is None or res["loglik"] > best["loglik"]:
                best = res

    state_mean = np.sum(best["weights"] * best["means"], axis=1)
    order = np.argsort(state_mean)
    pi = best["pi"][order]
    A = best["A"][np.ix_(order, order)]
    weights = best["weights"][order]
    means = best["means"][order]
    sds = best["sds"][order]
    gamma = best["gamma"][:, order]
    comp_resp = best["component_resp"][:, order, :]
    vit = viterbi_path_mixture(y, pi, A, weights, means, sds)
    return MixtureHMMFit(K=K, M=M, loglik=best["loglik"], pi=pi, A=A,
                         weights=weights, means=means, sds=sds, gamma=gamma,
                         component_resp=comp_resp, viterbi=vit,
                         n_iter=best["n_iter"], loglik_trace=best["trace"])


def validate_mixture_fit(fit: MixtureHMMFit, T: int) -> None:
    """Assertions on a fitted mixture-emission HMM."""
    assert fit.gamma.shape == (T, fit.K), "gamma has wrong shape"
    assert fit.component_resp.shape == (T, fit.K, fit.M), "component_resp has wrong shape"
    assert np.allclose(fit.A.sum(axis=1), 1.0, atol=1e-8), "rows of A must sum to 1"
    assert np.allclose(fit.pi.sum(), 1.0, atol=1e-8), "pi must sum to 1"
    assert np.allclose(fit.weights.sum(axis=1), 1.0, atol=1e-8), \
        "mixture weights must sum to 1 within state"
    assert np.all(fit.sds > 0), "estimated sds must be positive"
    assert np.allclose(fit.gamma.sum(axis=1), 1.0, atol=1e-8), \
        "posterior probabilities must sum to 1 at each time"
    assert np.all(np.isfinite(fit.gamma)), "gamma contains NaN/Inf"
    assert np.all(np.isfinite(fit.component_resp)), "component_resp contains NaN/Inf"
    assert np.all(np.isfinite([fit.loglik, *fit.means.ravel(), *fit.sds.ravel()])), \
        "parameters contain NaN/Inf"
    tr = np.array(fit.loglik_trace)
    assert np.all(np.diff(tr) >= -1e-6 * (1.0 + np.abs(tr[:-1]))), \
        "mixture EM log-likelihood is not non-decreasing (beyond tolerance)"
