"""Simulate a circular chromosome with a known origin and terminus.

From the origin to the terminus (clockwise) the forward strand is the
leading strand, so G is more frequent than C by the factor
(1 + skew) / (1 - skew); from the terminus back to the origin it is the
other way round.

Usage:
    uv run python simulate.py --skew 0.02 --seed 1 -o sim.fasta
"""

import argparse

import numpy as np

BASES = np.array(list("ACGT"))


def base_probabilities(skew, gc):
    """Return the probabilities of the bases on a leading stretch.

    Parameters
    ----------
    skew : float
        Strength of the skew; a negative value gives a lagging stretch.
    gc : float
        GC content, between 0 and 1.

    Returns
    -------
    list of float
        Probabilities of A, C, G and T, in this order.
    """
    return [
        (1 - gc) / 2,
        gc / 2 * (1 - skew),
        gc / 2 * (1 + skew),
        (1 - gc) / 2,
    ]


def simulate_genome(length, ori, skew, rng, gc=0.5):
    """Simulate a circular genome with a planted origin.

    Parameters
    ----------
    length : int
        Genome length in bases.
    ori : int
        1-based position of the origin.
    skew : float
        Strength of the skew, (G - C) / (G + C) on the leading strand.
    rng : numpy.random.Generator
        Source of randomness, so the caller decides the seed.
    gc : float, optional
        GC content, by default 0.5.

    Returns
    -------
    sequence : str
        The genome, forward strand.
    ter : int
        1-based position of the terminus, half a genome away from the
        origin.
    """
    ter = (ori - 1 + length // 2) % length + 1
    leading = base_probabilities(skew, gc)
    lagging = base_probabilities(-skew, gc)
    # Position p (0-based) is on the leading stretch if it lies less
    # than half a genome clockwise from the origin: from ori up to,
    # but not including, ter.
    offset = (np.arange(length) - (ori - 1)) % length
    on_leading = offset < length // 2
    bases = np.where(
        on_leading,
        rng.choice(4, size=length, p=leading),
        rng.choice(4, size=length, p=lagging),
    )
    return "".join(BASES[bases]), ter


def main(argv=None):
    """Read the options, simulate one genome and write it as FASTA.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments; by default those of the program, so
        tests can pass their own.
    """
    parser = argparse.ArgumentParser(
        description="Simulate a genome with a planted origin."
    )
    parser.add_argument(
        "--length", type=int, default=100_000,
        help="genome length in bases (default: 100000)",
    )
    parser.add_argument(
        "--ori", type=int, default=None,
        help="1-based position of the origin (default: random, drawn "
        "with the seed)",
    )
    parser.add_argument(
        "--skew", type=float, default=0.02,
        help="strength of the skew, (G - C) / (G + C) on the G-rich half "
        "(default: 0.02)",
    )
    parser.add_argument(
        "--gc", type=float, default=0.5,
        help="GC content, between 0 and 1 (default: 0.5)",
    )
    parser.add_argument(
        "--seed", type=int, default=1,
        help="random seed; the same seed and options give the same "
        "genome (default: 1)",
    )
    parser.add_argument(
        "-o", "--output", required=True, metavar="FASTA",
        help="output FASTA file",
    )
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    if args.ori is not None:
        ori = args.ori
    else:
        ori = int(rng.integers(1, args.length + 1))
    sequence, ter = simulate_genome(args.length, ori, args.skew, rng, args.gc)
    with open(args.output, "w") as handle:
        handle.write(
            f">simulated ori={ori} ter={ter} "
            f"skew={args.skew} seed={args.seed}\n"
        )
        for start in range(0, len(sequence), 80):
            end = start + 80
            handle.write(sequence[start:end] + "\n")


if __name__ == "__main__":
    main()
