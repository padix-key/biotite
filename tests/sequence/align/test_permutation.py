# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import numpy as np
import biotite.sequence.align as align


def test_random_permutation_modulo():
    """
    Test if the LCG of :class:`RandomPermutation`, omitting the explicit
    modulo due to automatic bit truncation, gives the same result as
    a slower implementation using explicit modulo and an unlimited
    integer type.
    """
    LCG_A = align.RandomPermutation.LCG_A
    LCG_C = align.RandomPermutation.LCG_C
    LCG_M = int(2**64)
    SEQ_LENGTH = 10_000
    
    np.random.seed(0)
    kmers = np.random.randint(np.iinfo(np.int64).max + 1, size=SEQ_LENGTH)

    ref_order = [
        (LCG_A * kmer.item() + LCG_C) % LCG_M
        for kmer in kmers
    ]

    permutation = align.RandomPermutation()
    test_order = permutation.permute(kmers)
    # Convert back to uint64
    # to revert additional step in RandomPermutation
    test_order = test_order.view(np.uint64)

    assert test_order.tolist() == ref_order


def test_random_permutation_randomness():
    """
    A simple test of randomness from :class:`RandomPermutation`:
    At each point in the permuted sequence , there should be an almost
    equal number of positive and negative integers
    """
    SEQ_LENGTH = 1_000_000

    np.random.seed(0)
    ref_kmers = np.random.randint(np.iinfo(np.int64).max + 1, size=SEQ_LENGTH)
    ref_distribution = _create_distribution(ref_kmers)

    test_kmers = np.arange(int(2**62), int(2**62)+SEQ_LENGTH)
    test_distribution = _create_distribution(test_kmers)

def _create_distribution(kmers):
    BINS = np.arange(-30, 30)

    permutation = align.RandomPermutation()
    order = permutation.permute(kmers)
    cum_sign = np.cumsum(np.sign(order))
    distribution, _ = np.histogram(cum_sign, bins=BINS)
    return distribution / len(kmers)