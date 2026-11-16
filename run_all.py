"""Reproduce every figure and table of the report.

Usage:
    uv run python run_all.py
"""

import time

import download_data
from experiments import (
    RESULTS_DIR,
    plot_skew_experiment,
    real_genome,
    skew_experiment,
    write_tsv,
)

SEED = 2026


def main():
    """Download the data, run both experiments, write results/."""
    RESULTS_DIR.mkdir(exist_ok=True)

    start = time.perf_counter()
    paths = download_data.main()
    print(f"data: {time.perf_counter() - start:.0f} s")

    start = time.perf_counter()
    rows = skew_experiment(
        skews=[0.002, 0.005, 0.01, 0.02, 0.05],
        length=500_000,
        repeats=30,
        seed=SEED,
    )
    write_tsv(rows, RESULTS_DIR / "validation.tsv")
    plot_skew_experiment(rows, RESULTS_DIR / "validation.png")
    print(f"validation: {time.perf_counter() - start:.0f} s")

    start = time.perf_counter()
    real = [real_genome(path, path.stem) for path in paths]
    write_tsv(real, RESULTS_DIR / "real_data.tsv")
    print(f"real data: {time.perf_counter() - start:.0f} s")


if __name__ == "__main__":
    main()
