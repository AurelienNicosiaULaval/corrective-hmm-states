"""Simulation of the true persistent 2-regime HMM with a mixture emission
in regime 1 (refined sub-states 1a, 1b) and a Gaussian emission in regime 2.

True model
----------
Hidden aggregated chain Z_t in {1, 2} with transition matrix

    Gamma0 = [[0.93, 0.07],
              [0.05, 0.95]]

Emissions:
    Y_t | Z_t = 1 ~ 0.70 N(-1.20, 0.55^2) + 0.30 N(1.10, 0.45^2)
    Y_t | Z_t = 2 ~ N(3.80, 0.65^2)

The refined state R_t in {0, 1, 2} encodes {1a, 1b, 2}.
"""
from __future__ import annotations

import numpy as np

# True parameters (module-level constants so every script uses the same truth)
GAMMA0 = np.array([[0.93, 0.07],
                   [0.05, 0.95]])
MIX_WEIGHTS_REGIME1 = np.array([0.70, 0.30])
MEANS_REFINED = np.array([-1.20, 1.10, 3.80])   # 1a, 1b, 2
SDS_REFINED = np.array([0.55, 0.45, 0.65])
REFINED_LABELS = np.array(["1a", "1b", "2"])


def stationary_distribution(A: np.ndarray) -> np.ndarray:
    """Stationary distribution of a transition matrix (left eigenvector)."""
    vals, vecs = np.linalg.eig(A.T)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.real(vecs[:, idx])
    pi = np.abs(pi)
    return pi / pi.sum()


def simulate_hmm(T: int, seed: int):
    """Simulate the true model.

    Returns
    -------
    y : (T,) observations
    z : (T,) aggregated hidden state in {0, 1}  (regimes 1, 2)
    refined : (T,) refined hidden state in {0, 1, 2}  (1a, 1b, 2)
    """
    rng = np.random.default_rng(seed)
    pi0 = stationary_distribution(GAMMA0)
    z = np.zeros(T, dtype=int)
    refined = np.zeros(T, dtype=int)
    y = np.zeros(T)
    z[0] = rng.choice(2, p=pi0)
    for t in range(1, T):
        z[t] = rng.choice(2, p=GAMMA0[z[t - 1]])
    for t in range(T):
        if z[t] == 0:
            refined[t] = rng.choice(2, p=MIX_WEIGHTS_REGIME1)  # 0 -> 1a, 1 -> 1b
        else:
            refined[t] = 2
    y = rng.normal(MEANS_REFINED[refined], SDS_REFINED[refined])
    return y, z, refined
