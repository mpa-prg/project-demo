"""Tests of the core in skew.py: toy data, invariants, brute force."""

import numpy as np
import pytest

from skew import circular_distance, cumulative_skew, find_origin


def test_cumulative_skew_toy():
    """Toy data with an answer computed by hand."""
    assert list(cumulative_skew("GGCATCC")) == [0, 1, 2, 1, 1, 1, 0, -1]


def test_lower_case_counts():
    """Lower-case bases count as upper-case ones."""
    assert list(cumulative_skew("ggc")) == list(cumulative_skew("GGC"))


@pytest.mark.parametrize("seed", range(5))
def test_last_value_is_total(seed):
    """Invariant: the last value is the total G minus C."""
    rng = np.random.default_rng(seed)
    sequence = "".join(rng.choice(list("ACGT"), size=1000))
    total = sequence.count("G") - sequence.count("C")
    assert cumulative_skew(sequence)[-1] == total


def brute_force_minimum(sequence):
    """Brute force: a slow, obviously correct version of the origin."""
    n = len(sequence)
    total = sequence.count("G") - sequence.count("C")
    best_value, best_i = 0.0, 0
    for i in range(1, n):
        prefix = sequence[:i]
        value = prefix.count("G") - prefix.count("C") - i * total / n
        if value < best_value:
            best_value, best_i = value, i
    return best_i + 1


@pytest.mark.parametrize("seed", range(20))
def test_origin_matches_brute_force(seed):
    """The fast version agrees with the brute force."""
    rng = np.random.default_rng(seed)
    sequence = "".join(rng.choice(list("ACGT"), size=200))
    assert find_origin(sequence)[0] == brute_force_minimum(sequence)


@pytest.mark.parametrize(
    "a, b, n, expected",
    [
        (1, 1, 10, 0),
        (1, 10, 10, 1),  # neighbours across the end of the circle
        (2, 7, 10, 5),  # exactly opposite
        (3, 5, 10, 2),
    ],
)
def test_circular_distance(a, b, n, expected):
    """Edge cases of the circle, in both directions."""
    assert circular_distance(a, b, n) == expected
    assert circular_distance(b, a, n) == expected


@pytest.mark.parametrize("shift", [1, 37, 150])
def test_rotation_moves_origin_with_it(shift):
    """Invariant: the origin moves with a rotated sequence."""
    rng = np.random.default_rng(0)
    sequence = "".join(rng.choice(list("ACGT"), size=200))
    rotated = sequence[shift:] + sequence[:shift]
    ori, _ = find_origin(sequence)
    ori_rotated, _ = find_origin(rotated)
    assert (ori_rotated + shift - 1) % len(sequence) + 1 == ori
