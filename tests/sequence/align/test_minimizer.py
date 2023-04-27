# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import itertools
import numpy as np
import pytest
import biotite.sequence as seq
import biotite.sequence.align as align


@pytest.mark.parametrize(
    "seed, window, from_sequence",
    itertools.product(
        range(20),
        [None, 2, 5, 10, 25],
        [False, True]
    )
)
def test_minimize(seed, window, from_sequence):
    K = 10
    LENGTH = 1000

    sequence = seq.NucleotideSequence(ambiguous=False)
    np.random.seed(seed)
    sequence.code = np.random.randint(len(sequence.alphabet), size=LENGTH)
    kmer_alph = align.KmerAlphabet(sequence.alphabet, K)
    kmers = kmer_alph.create_kmers(sequence.code)

    if window is None:
        ref_minimizer_pos = np.array([np.argmin(kmers)])
        ref_minimizers = np.array([kmers[np.argmin(kmers)]])
    else:
        # Use an inefficient but simple algorithm for comparison
        ref_minimizer_pos = np.array([
            np.argmin(kmers[i : i + window]) + i
            for i in range(len(kmers) - (window - 1))
        ])
        # Remove duplicates
        ref_minimizer_pos = np.unique(ref_minimizer_pos)
        ref_minimizers = kmers[ref_minimizer_pos]
    
    minimizer = align.Minimizer(kmer_alph)
    if from_sequence:
        test_minimizer_pos, test_minimizers = minimizer.minimize(
            sequence, window
        )
    else:
        test_minimizer_pos, test_minimizers = minimizer.minimize_kmers(
            kmers, window
        )
    
    assert test_minimizer_pos.tolist() == ref_minimizer_pos.tolist()
    assert test_minimizers.tolist() == ref_minimizers.tolist()
