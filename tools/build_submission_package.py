"""Build a clean local submission package for Stat/Wiley.

Run from the project root:
    python3 tools/build_submission_package.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUBMISSION_DIR = ROOT / "submission"
PACKAGE_DIR = SUBMISSION_DIR / "stat_hmm_submission_package"
ZIP_BASE = SUBMISSION_DIR / "stat_hmm_submission_package"


def copy_file(src_rel: str, dst_rel: str | None = None) -> None:
    src = ROOT / src_rel
    dst = PACKAGE_DIR / (dst_rel or src_rel)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_tree_files(src_rel: str, patterns: tuple[str, ...]) -> None:
    src_dir = ROOT / src_rel
    for pattern in patterns:
        for src in src_dir.glob(pattern):
            if src.is_file():
                rel = src.relative_to(ROOT)
                copy_file(str(rel))


def write_package_readme() -> None:
    text = """# Submission package

Manuscript: Extra Hidden States as Kullback-Leibler Corrections under Emission Misspecification in Hidden Markov Models

Contents:

- `article/main.pdf`: compiled manuscript for peer review.
- `article/main.tex`: LaTeX source.
- `article/references.bib`: BibTeX bibliography.
- `article/figures/`: figures used in the manuscript.
- `article/tables/simulation_table.tex`: generated LaTeX table.
- `article/tables/replication_summary_table.tex`: generated replication table.
- `supplement/supporting_information.pdf`: supporting-information PDF.
- `supplement/supporting_information.tex`: supporting-information source.
- `numerics/`: reproducible simulation code, tests, requirements and generated outputs.

Reproduce the numerical results from the project root with:

```bash
python3 numerics/scripts/run_review_simulation.py
```

The current manuscript uses no external data.
"""
    (PACKAGE_DIR / "README_submission.md").write_text(text)


def main() -> None:
    if PACKAGE_DIR.exists():
        shutil.rmtree(PACKAGE_DIR)
    SUBMISSION_DIR.mkdir(parents=True, exist_ok=True)

    for src in [
        "README.md",
        "article/main.tex",
        "article/main.pdf",
        "article/references.bib",
        "article/tables/simulation_table.tex",
        "article/tables/replication_summary_table.tex",
        "supplement/supporting_information.tex",
        "supplement/supporting_information.pdf",
        "supplement/tables/review_summary_table.tex",
        "numerics/README.md",
        "numerics/requirements.txt",
        "numerics/src/__init__.py",
        "numerics/src/simulate.py",
        "numerics/src/hmm_em.py",
        "numerics/src/mixture_hmm_em.py",
        "numerics/src/diagnostics.py",
        "numerics/src/plotting.py",
        "numerics/scripts/run_review_simulation.py",
        "numerics/tests/test_numerics.py",
        "numerics/output/review_gaussian_fit_summary.csv",
        "numerics/output/review_mixture_fit_summary.csv",
        "numerics/output/review_model_comparison.csv",
        "numerics/output/review_replications.csv",
        "numerics/output/review_summary_by_T.csv",
        "numerics/output/review_metadata.json",
        "numerics/output/crosstab_fit2.csv",
        "numerics/output/crosstab_fit3.csv",
        "numerics/output/residual_acf_fit2.csv",
        "numerics/output/residual_acf_fit3.csv",
    ]:
        copy_file(src)

    copy_tree_files("article/figures", ("*.pdf", "*.png"))
    write_package_readme()

    zip_path = shutil.make_archive(str(ZIP_BASE), "zip", root_dir=PACKAGE_DIR)
    print(f"wrote {PACKAGE_DIR}")
    print(f"wrote {zip_path}")


if __name__ == "__main__":
    main()
