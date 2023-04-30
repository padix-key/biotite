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
    Provides an order for *k-mers*, usually used by *submer* rules
    such as :class:`MinimizerRule`.
    The rule how such order is computed depends on the concrete
    subclass of this abstract base class.

    Without a :class:`Permutation` *submer* rules usually resort to
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
    """

    def __init__(self, kmer_alphabet, seed=None):
        super().__init__(kmer_alphabet)
        rng = np.random.default_rng(seed)
        self._permutation_table = rng.permutation(len(kmer_alphabet))
    
    def permute(self, kmers):
        return self._permutation_table[kmers]