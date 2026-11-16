"""Validation on simulated genomes and the result on the real genome."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

from simulate import simulate_genome
from skew import circular_distance, cumulative_skew, find_origin

plt.switch_backend("Agg")  # write files only, no window

RESULTS_DIR = Path(__file__).parent / "results"
REFERENCE = Path(__file__).parent / "reference" / "oriC_annotation.tsv"


def skew_experiment(skews, length, repeats, seed):
    """Measure the error of the origin estimate on simulated genomes.

    Parameters
    ----------
    skews : list of float
        Skews to simulate.
    length : int
        Genome length in bases.
    repeats : int
        Number of genomes per skew, each with a random origin.
    seed : int
        Seed; every genome gets its own generator derived from it.

    Returns
    -------
    list of dict
        One row per genome: skew, repeat, true_ori, found_ori and
        error_bp (the circular distance between them).
    """
    rows = []
    for i, skew in enumerate(skews):
        for repeat in range(repeats):
            rng = np.random.default_rng([seed, i, repeat])
            ori = int(rng.integers(1, length + 1))
            sequence, _ = simulate_genome(length, ori, skew, rng)
            found, _ = find_origin(sequence)
            rows.append(
                {
                    "skew": skew,
                    "repeat": repeat,
                    "true_ori": ori,
                    "found_ori": found,
                    "error_bp": circular_distance(found, ori, length),
                }
            )
    return rows


def plot_skew_experiment(rows, path):
    """Plot the error of every run and the median per skew.

    Parameters
    ----------
    rows : list of dict
        Rows returned by `skew_experiment`.
    path : pathlib.Path
        Output image file.
    """
    skews = sorted({row["skew"] for row in rows})
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for skew in skews:
        errors = [
            max(row["error_bp"], 1) for row in rows if row["skew"] == skew
        ]
        ax.scatter(
            [skew] * len(errors), errors, color="tab:blue", alpha=0.4, s=15
        )
        ax.scatter(
            skew, np.median(errors), color="tab:orange", marker="_", s=300
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(skews, [str(skew) for skew in skews])
    ax.minorticks_off()
    ax.set_xlabel("planted skew")
    ax.set_ylabel("distance to the true origin (bp)")
    ax.set_title("Simulated genomes, 500 kb (orange: median)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def real_genome(fasta, accession):
    """Estimate the origin of a real genome and compare it with oriC.

    Also plots the cumulative skew with the annotated oriC to
    results/skew_<accession>.png.

    Parameters
    ----------
    fasta : pathlib.Path
        FASTA file with exactly one record.
    accession : str
        Accession of the genome in reference/oriC_annotation.tsv.

    Returns
    -------
    dict
        One row of the results table: accession, length, found_ori,
        found_ter, oriC_start, oriC_end and distance_bp (0 if the
        estimate lies inside the annotated oriC).
    """
    record = SeqIO.read(fasta, "fasta")
    sequence = str(record.seq)
    ori, ter = find_origin(sequence)
    with open(REFERENCE) as handle:
        annotation = next(
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["accession"] == accession
        )
    start, end = int(annotation["start"]), int(annotation["end"])
    if start <= ori <= end:
        distance = 0
    else:
        distance = min(
            circular_distance(ori, start, len(sequence)),
            circular_distance(ori, end, len(sequence)),
        )

    skew = cumulative_skew(sequence)
    step = 1000
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(
        np.arange(0, len(skew), step) / 1e6, skew[::step], color="tab:blue"
    )
    ax.axvline(
        start / 1e6, color="tab:orange", linestyle="--", label="annotated oriC"
    )
    ax.set_xlabel("position (Mb)")
    ax.set_ylabel("cumulative skew (G - C)")
    ax.set_title(accession)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"skew_{accession}.png", dpi=150)
    plt.close(fig)

    return {
        "accession": accession,
        "length": len(sequence),
        "found_ori": ori,
        "found_ter": ter,
        "oriC_start": start,
        "oriC_end": end,
        "distance_bp": distance,
    }


def write_tsv(rows, path):
    """Write rows as a TSV file with a header.

    Parameters
    ----------
    rows : list of dict
        Rows with the same keys; the keys of the first row are the
        columns.
    path : pathlib.Path
        Output file.
    """
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
