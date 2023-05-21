# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.sequence.align"
__author__ = "Patrick Kunzmann"
__all__ = ["Permutation", "RandomPermutation"]

import abc
import numpy as np


class Permutation(metaclass=abc.ABCMeta):
    """
    Provides an order for *k-mers*, usually used by *k-mer* subset rules
    such as :class:`MinimizerRule`.
    The method how such order is computed depends on the concrete
    subclass of this abstract base class.

    Without a :class:`Permutation` subset rules usually resort to
    the symbol order in the :class:`KmerAlphabet`.
    That order is often the lexicographical order, which is known to
    yield suboptimal *k-mer* selection many cases
    :footcite:`Roberts2004`.
    """

    @abc.abstractmethod
    def permute(self, kmers):
        """
        Give the given *k-mers* a new order.

        Parameters
        ----------
        kmers : ndarray, dtype=np.int64
            The *k-mers* to reorder given as *k-mer* code.
        
        Returns
        -------
        order : ndarray, dtype=np.int64
            The new order, i.e. a *k-mer* A is smaller than *k-mer* B,
            if ``order[A] < order[B]``
        """
        pass


class RandomPermutation(Permutation):
    """
    __init__(kmer_alphabet, seed=None)

    Provide a randomized order for *k-mers* from a given
    :class:`KmerAlphabet`.

    Parameters
    ----------
    kmer_alphabet : KmerAlphabet
        The *k-mer* alphabet that defines the range of possible *k-mers*
        that should be permuted.
    seed : int, optional
        The seed for the random number generator that creates the
        new order.
        By default, the seed is randomly chosen.
    
    Notes
    -----

    This class uses a lookup table to achieve permutation.
    Hence, the memory consumption is :math:`8 n^k` bytes,
    where :math:`n` is the size of the base alphabet and :math:`k` is
    the *k-mer* size.

    Examples
    --------

    >>> kmer_alph = KmerAlphabet(NucleotideSequence.alphabet_unamb, k=2)
    >>> permutation = RandomPermutation(kmer_alph, seed=0)
    >>> # k-mer codes representing the k-mers from 'AA' to 'TT'
    >>> # in lexicographic order
    >>> kmer_codes = np.arange(len(kmer_alph))
    >>> print(kmer_codes)
    [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
    >>> # Shuffle order of these k-mer codes using the permutation
    >>> print(permutation.permute(kmer_codes))
    [ 2 11  3 10  0  4  7  5 14 12  6  9 13  8  1 15]
    """

    def __init__(self, kmer_alphabet, seed=None):
        rng = np.random.default_rng(seed)
        self._permutation_table = rng.permutation(len(kmer_alphabet))
    
    def permute(self, kmers):
        return self._permutation_table[kmers]


class FrequencyPermutation(Permutation):
    """
    __init__(kmer_alphabet, counts)

    Provide an order for *k-mers* from a given
    :class:`KmerAlphabet`, such that less frequent *k-mers* are smaller
    that more frequent *k-mers*.
    The frequency of each *k-mer* can either be given directly via the
    constructor or can be computed from a :class:`KmerTable` via
    :meth:`from_table()`.

    Parameters
    ----------
    kmer_alphabet : KmerAlphabet, length=n
        The *k-mer* alphabet that defines the range of possible *k-mers*
        that should be permuted.
    counts : ndarray, shape=(n,), dtype=np.int64
        The absolute frequency, i.e. the number of occurrences, of each
        *k-mer* in `kmer_alphabet` in the sequence database of interest.
        ``counts[c] = f``, where ``c`` is the *k-mer* code and ``f`` is
        the corresponding frequency.
    
    Notes
    -----

    In actual sequences some sequence patterns appear in high quantity.
    When selecting a subset of *k-mers*, e.g. via
    :class:`MinimizerRule`, it is desireable to select the low-frequency
    *informative* *k-mers* to avoid spurious matches.
    To achieve such selection this class can be used.

    This class uses a lookup table to achieve permutation.
    Hence, the memory consumption is :math:`8 n^k` bytes,
    where :math:`n` is the size of the base alphabet and :math:`k` is
    the *k-mer* size.

    Examples
    --------

    >>> TODO
    """
    def __init__(self, kmer_alphabet, counts):
        if len(kmer_alphabet) != len(counts):
            raise IndexError(
                f"The k-mer alphabet has {len(kmer_alphabet)} k-mers, "
                f"but {len(counts)} counts were given"
            )
        order = np.argsort(counts)
        # 'order' maps a permutation a k-mer
        # Via fancy indexing with itself, the result maps a k-mer
        # to its permutation
        self._permutation_table = order[order]
    
    def permute(self, kmers):
        return self._permutation_table[kmers]