"""Tests of the simulator and of the core on simulated genomes."""

import numpy as np
import pytest

from simulate import simulate_genome
from skew import circular_distance, find_origin


def test_same_seed_same_genome():
    """The same seed gives the same genome."""
    a, _ = simulate_genome(1000, 100, 0.1, np.random.default_rng(7))
    b, _ = simulate_genome(1000, 100, 0.1, np.random.default_rng(7))
    assert a == b


def test_different_seed_different_genome():
    """A different seed gives a different genome."""
    a, _ = simulate_genome(1000, 100, 0.1, np.random.default_rng(7))
    b, _ = simulate_genome(1000, 100, 0.1, np.random.default_rng(8))
    assert a != b


def test_length_and_alphabet():
    """The genome has the requested length and only A, C, G, T."""
    sequence, _ = simulate_genome(999, 1, 0.1, np.random.default_rng(1))
    assert len(sequence) == 999
    assert set(sequence) <= set("ACGT")


def test_terminus_is_opposite():
    """The terminus lies half a genome away, across the end."""
    _, ter = simulate_genome(1000, 900, 0.1, np.random.default_rng(1))
    assert ter == 400


def test_g_rich_after_origin():
    """G > C from the origin to the terminus, C > G after it."""
    rng = np.random.default_rng(1)
    sequence, ter = simulate_genome(100_000, 1, 0.2, rng)
    split = ter - 1
    leading, lagging = sequence[:split], sequence[split:]
    assert leading.count("G") > leading.count("C")
    assert lagging.count("G") < lagging.count("C")


@pytest.mark.parametrize("ori", [1, 50_000, 199_000])
def test_planted_origin_is_found(ori):
    """The core finds a planted origin when the signal is strong."""
    rng = np.random.default_rng(3)
    sequence, ter = simulate_genome(200_000, ori, 0.1, rng)
    found_ori, found_ter = find_origin(sequence)
    assert circular_distance(found_ori, ori, len(sequence)) < 1000
    assert circular_distance(found_ter, ter, len(sequence)) < 1000
