"""GC skew and the position of the replication origin.

On the leading strand of a bacterial chromosome, G is slightly more
frequent than C. Read along the forward strand, the skew (G - C) is
therefore negative from the terminus to the origin and positive from
the origin to the terminus. The cumulative skew reaches its minimum at
the origin and its maximum at the terminus.
"""

import numpy as np


def cumulative_skew(sequence):
    """Return the cumulative skew of a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence. Lower-case bases count as upper-case ones.

    Returns
    -------
    numpy.ndarray
        Integer array of length len(sequence) + 1. Element i is
        (number of G) - (number of C) in the first i bases, so
        element 0 is always 0.
    """
    bases = np.array(list(sequence.upper()))
    # +1 for G, -1 for C, 0 for anything else.
    steps = (bases == "G").astype(np.int64) - (bases == "C")
    return np.concatenate(([0], np.cumsum(steps)))


def find_origin(sequence):
    """Estimate the origin and the terminus of a circular chromosome.

    Parameters
    ----------
    sequence : str
        The chromosome, forward strand. Lower-case bases count as
        upper-case ones.

    Returns
    -------
    ori : int
        1-based position of the first base after the minimum of the
        cumulative skew.
    ter : int
        1-based position of the first base after its maximum.
        Ties are broken by the leftmost position.

    Notes
    -----
    The chromosome is circular, but the total G - C of a whole genome is
    rarely exactly 0, so the cumulative skew would jump where the
    sequence file starts. Subtracting the linear trend i * total / n
    makes the curve end where it began, so the result does not depend on
    where the file starts.
    """
    skew = cumulative_skew(sequence)
    n = len(sequence)
    circular = skew[:n] - np.arange(n) * skew[n] / n
    ori = int(np.argmin(circular)) + 1
    ter = int(np.argmax(circular)) + 1
    return ori, ter


def circular_distance(a, b, n):
    """Return the distance of two positions on a circle.

    Parameters
    ----------
    a, b : int
        Positions on the circle (both 0-based or both 1-based).
    n : int
        Length of the circle.

    Returns
    -------
    int
        The shorter of the two distances around the circle.
    """
    d = abs(a - b) % n
    return min(d, n - d)
