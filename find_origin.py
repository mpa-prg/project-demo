"""Run the core on a FASTA file: the origin and terminus of each record.

Usage:
    uv run python find_origin.py genome.fasta -o result.tsv
"""

import argparse

from Bio import SeqIO

from skew import find_origin


def main(argv=None):
    """Read the options, run the core on every record, write a TSV.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments; by default those of the program, so
        tests can pass their own.
    """
    parser = argparse.ArgumentParser(
        description="Estimate the replication origin from the GC skew."
    )
    parser.add_argument(
        "fasta", metavar="FASTA",
        help="input FASTA, one circular chromosome per record",
    )
    parser.add_argument(
        "-o", "--output", required=True, metavar="TSV",
        help="output TSV file",
    )
    args = parser.parse_args(argv)

    with open(args.output, "w") as out:
        out.write("record\tlength\tori\tter\n")
        for record in SeqIO.parse(args.fasta, "fasta"):
            ori, ter = find_origin(str(record.seq))
            out.write(f"{record.id}\t{len(record)}\t{ori}\t{ter}\n")


if __name__ == "__main__":
    main()
