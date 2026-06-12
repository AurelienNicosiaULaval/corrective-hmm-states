"""Local validation checks for the Stat/Wiley submission package.

Run from the project root:
    python3 tools/validate_submission.py
"""
from __future__ import annotations

import re
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_FILES = [
    ROOT / "README.md",
    ROOT / "article" / "main.tex",
    ROOT / "article" / "references.bib",
    ROOT / "article" / "stat_submission_notes.md",
    ROOT / "article" / "tables" / "simulation_table.tex",
    ROOT / "supplement" / "supporting_information.tex",
    ROOT / "supplement" / "tables" / "review_summary_table.tex",
    ROOT / "numerics" / "README.md",
    ROOT / "numerics" / "requirements.txt",
    ROOT / "numerics" / "src" / "hmm_em.py",
    ROOT / "numerics" / "src" / "mixture_hmm_em.py",
    ROOT / "numerics" / "src" / "diagnostics.py",
    ROOT / "numerics" / "src" / "plotting.py",
    ROOT / "numerics" / "src" / "simulate.py",
    ROOT / "numerics" / "scripts" / "run_review_simulation.py",
    ROOT / "numerics" / "tests" / "test_numerics.py",
]
REQUIRED_FILES = [
    ROOT / "article" / "main.tex",
    ROOT / "article" / "main.pdf",
    ROOT / "article" / "references.bib",
    ROOT / "article" / "figures" / "simulation_density.pdf",
    ROOT / "article" / "figures" / "simulation_density.png",
    ROOT / "article" / "figures" / "state_decoding_heatmap.pdf",
    ROOT / "article" / "figures" / "state_decoding_heatmap.png",
    ROOT / "article" / "figures" / "posterior_corrective_state.pdf",
    ROOT / "article" / "figures" / "posterior_corrective_state.png",
    ROOT / "article" / "figures" / "bic_boxplot.pdf",
    ROOT / "article" / "figures" / "bic_boxplot.png",
    ROOT / "article" / "figures" / "bic_by_T.pdf",
    ROOT / "article" / "figures" / "bic_by_T.png",
    ROOT / "article" / "figures" / "residual_acf.pdf",
    ROOT / "article" / "figures" / "residual_acf.png",
    ROOT / "article" / "tables" / "simulation_table.tex",
    ROOT / "article" / "tables" / "replication_summary_table.tex",
    ROOT / "supplement" / "supporting_information.tex",
    ROOT / "supplement" / "supporting_information.pdf",
    ROOT / "supplement" / "tables" / "review_summary_table.tex",
    ROOT / "numerics" / "output" / "review_gaussian_fit_summary.csv",
    ROOT / "numerics" / "output" / "review_mixture_fit_summary.csv",
    ROOT / "numerics" / "output" / "review_model_comparison.csv",
    ROOT / "numerics" / "output" / "review_replications.csv",
    ROOT / "numerics" / "output" / "review_summary_by_T.csv",
    ROOT / "numerics" / "output" / "review_metadata.json",
    ROOT / "numerics" / "output" / "residual_acf_fit2.csv",
    ROOT / "numerics" / "output" / "residual_acf_fit3.csv",
]

FORBIDDEN_PATTERNS = [
    r"\[" + "email to be " + "added" + r"\]",
    "To be " + "completed",
    r"\b" + "TO" + "DO" + r"\b",
    r"\b" + "FIX" + "ME" + r"\b",
    r"\b" + "Chat" + "GPT" + r"\b",
    r"\b" + "AI-" + "generated" + r"\b",
    r"\b" + "assist" + "ant" + r"\b",
    r"\b" + "Co" + "dex" + r"\b",
]


def fail(message: str) -> None:
    print(f"FAIL: {message}")
    raise SystemExit(1)


def check_required_files() -> None:
    missing = [path for path in REQUIRED_FILES if not path.exists()]
    if missing:
        fail("missing required files:\n" + "\n".join(str(p) for p in missing))
    print(f"required files: {len(REQUIRED_FILES)} found")


def check_public_text() -> None:
    compiled = [(pat, re.compile(pat, flags=re.IGNORECASE)) for pat in FORBIDDEN_PATTERNS]
    hits: list[str] = []
    for path in PUBLIC_FILES:
        text = path.read_text(errors="ignore")
        for pattern, regex in compiled:
            if regex.search(text):
                hits.append(f"{path.relative_to(ROOT)}: {pattern}")
    if hits:
        fail("public-facing text contains unresolved or non-public markers:\n" + "\n".join(hits))
    print("public-facing text: no unresolved markers found")


def check_numerics_smoke_test() -> None:
    sys.path.insert(0, str(ROOT / "numerics"))
    import numpy as np
    from src.hmm_em import fit_hmm, validate_fit
    from src.mixture_hmm_em import fit_mixture_hmm, validate_mixture_fit
    from src.simulate import GAMMA0, simulate_hmm, stationary_distribution

    pi = stationary_distribution(GAMMA0)
    if not np.allclose(pi @ GAMMA0, pi, atol=1e-10):
        fail("stationary distribution check failed")

    y, _, _ = simulate_hmm(160, seed=20260611)
    fit = fit_hmm(y, 2, seed=20260613, n_starts=2, max_iter=80)
    validate_fit(fit, len(y))
    mix = fit_mixture_hmm(y, 2, M=2, seed=20260614, n_starts=2, max_iter=80)
    validate_mixture_fit(mix, len(y))
    print("numerics smoke test: passed")


def check_review_outputs() -> None:
    meta_path = ROOT / "numerics" / "output" / "review_metadata.json"
    metadata = json.loads(meta_path.read_text())
    if metadata.get("n_replications_per_T") != 200:
        fail("review simulation must contain 200 replications per sample size")
    if metadata.get("n_starts_replications", 0) < 10:
        fail("Gaussian replicated fits must use at least 10 EM starts")
    if metadata.get("n_starts_mixture_replications", 0) < 10:
        fail("mixture replicated fits must use at least 10 EM starts")
    if metadata.get("screen_iter_replications", 0) < 1:
        fail("Gaussian replicated fits must record deterministic start screening")
    if metadata.get("screen_iter_mixture_replications", 0) < 1:
        fail("mixture replicated fits must record deterministic start screening")
    if metadata.get("refine_top_replications", 0) < 1:
        fail("Gaussian replicated fits must refine screened starts")
    if metadata.get("refine_top_mixture_replications", 0) < 1:
        fail("mixture replicated fits must refine screened starts")
    expected_T = [500, 1500, 5000]
    if metadata.get("T_grid") != expected_T:
        fail(f"unexpected review T grid: {metadata.get('T_grid')}")
    rows = metadata.get("selection_summary", [])
    if len(rows) != len(expected_T):
        fail("selection summary must contain one row per T")
    for row in rows:
        if int(row.get("n_replications", -1)) != 200:
            fail("each T must contain 200 replications")
    print("review outputs: metadata has 200 replications per T, >=10 starts and deterministic screening")


def main() -> None:
    check_required_files()
    check_public_text()
    check_review_outputs()
    check_numerics_smoke_test()
    print("submission validation: passed")


if __name__ == "__main__":
    main()
