"""End-to-end test of the two command-line tools."""

import numpy as np

import find_origin
import simulate
from simulate import simulate_genome


def test_simulate_then_find(tmp_path):
    """Simulate a genome into a file, then find its origin in it."""
    fasta = tmp_path / "sim.fasta"
    result = tmp_path / "result.tsv"
    simulate.main(["--length", "50000", "--ori", "1234", "--skew", "0.1",
                   "--seed", "1", "-o", str(fasta)])
    find_origin.main([str(fasta), "-o", str(result)])

    header, row = result.read_text().splitlines()
    assert header == "record\tlength\tori\tter"
    _, length, ori, ter = row.split("\t")
    assert length == "50000"
    assert abs(int(ori) - 1234) < 500
    assert abs(int(ter) - 26_234) < 500  # half a genome from the origin


def test_one_row_per_record(tmp_path):
    """A FASTA file with two genomes gives one row for each of them."""
    fasta = tmp_path / "two.fasta"
    result = tmp_path / "result.tsv"
    planted = {"first": (5_000, 1), "second": (35_000, 2)}  # ori, seed
    with open(fasta, "w") as handle:
        for name, (ori, seed) in planted.items():
            rng = np.random.default_rng(seed)
            sequence, _ = simulate_genome(50_000, ori, 0.1, rng)
            handle.write(f">{name}\n{sequence}\n")
    find_origin.main([str(fasta), "-o", str(result)])

    lines = result.read_text().splitlines()[1:]  # without the header
    rows = [line.split("\t") for line in lines]
    assert [row[0] for row in rows] == ["first", "second"]
    for row, (ori, _) in zip(rows, planted.values()):
        assert abs(int(row[2]) - ori) < 500
