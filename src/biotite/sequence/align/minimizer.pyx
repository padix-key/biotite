# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.sequence.align"
__author__ = "Patrick Kunzmann"
__all__ = ["Minimizer"]

cimport cython
cimport numpy as np

import abc
import numpy as np
from .kmeralphabet import KmerAlphabet
from ..alphabet import AlphabetError


ctypedef np.int64_t int64
ctypedef np.uint32_t uint32


# Obtained from 'np.iinfo(np.int64).max'
DEF MAX_INT_64 = 9223372036854775807


class Minimizer:

    class Permutation(metaclass=abc.ABCMeta):

        @abc.abstractmethod
        def permute(self, kmers):
            pass


    def __init__(self, kmer_alphabet, permutation=None):
        self._kmer_alph = kmer_alphabet
        self._permutation = permutation
    

    def minimize(self, sequence, window=None, bint alphabet_check=True):
        if alphabet_check:
            if not self._kmer_alph.base_alphabet.extends(sequence.alphabet):
                raise ValueError(
                    "The sequence's alphabet does not fit the k-mer alphabet"
                )
        kmers = self._kmer_alph.create_kmers(sequence.code)
        return self.minimize_kmers(kmers, window)
    

    def minimize_kmers(self, kmers, window=None):
        if self._permutation is not None:
            kmers = self._permutation.permute(kmers)

        if window is None:
            min_index = np.argmin(kmers)
            return (
                # Create an array containing a single element
                np.array([min_index], dtype=np.uint32),
                # Also single element as array -> slice of size 1
                kmers[min_index : min_index+1],
            )
        else:
            if window < 2:
                raise ValueError("Window size must be at least 2")
            if len(kmers) < window:
                raise ValueError(
                    "The number of k-mers is smaller than the window size"
                )
            return _minimize(kmers, window)
    

@cython.boundscheck(False)
@cython.wraparound(False)
def _minimize(int64[:] kmers, uint32 window):
    """
    Implementation of the algorithm originally devised by
    Marcel van Herk.

    In this implmentation the frame is chosen differently:
    For a position 'x' the frame ranges from 'x' to 'x + window-1'
    instead of 'x - (window-1)/2' to 'x + (window-1)/2'.
    """
    cdef uint32 seq_i
    
    cdef uint32 n_windows = kmers.shape[0] - (window - 1)
    # Pessimistic array allocation size
    # -> Expect that every window has a new minimizer
    cdef uint32[:] mininizer_pos = np.empty(n_windows, dtype=np.uint32)
    cdef int64[:] minimizers = np.empty(n_windows, dtype=np.int64)
    # Counts the actual number of minimiers for later trimming
    cdef uint32 n_minimizers = 0

    # Variables for the position of the previous cumulative minimum
    # Assign an value that can never occur for the start,
    # as in the beginning there is no previous value
    cdef uint32 prev_argcummin = kmers.shape[0]
    # Variables for the position of the current cumulative minimum
    cdef uint32 combined_argcummin, forward_argcummin, reverse_argcummin
    # Variables for the current cumulative minimum
    cdef int64 combined_cummin, forward_cummin, reverse_cummin
    # Variables for cumulative minima at all positions
    cdef uint32[:] forward_argcummins = _chunk_wise_forward_argcummin(
        kmers, window
    )
    cdef uint32[:] reverse_argcummins = _chunk_wise_reverse_argcummin(
        kmers, window
    )

    for seq_i in range(n_windows):
        forward_argcummin = forward_argcummins[seq_i + window - 1]
        reverse_argcummin = reverse_argcummins[seq_i]
        forward_cummin = kmers[forward_argcummin]
        reverse_cummin = kmers[reverse_argcummin]
        
        # At ties the leftmost position is taken,
        # which stems from the reverse pass
        if forward_cummin < reverse_cummin:
            combined_argcummin = forward_argcummin
            combined_cummin = forward_cummin
        else:
            combined_argcummin = reverse_argcummin
            combined_cummin = reverse_cummin
        
        if combined_argcummin != prev_argcummin:
            # A new minimizer is observed
            # -> append it to return value
            mininizer_pos[n_minimizers] = combined_argcummin
            minimizers[n_minimizers] = combined_cummin
            n_minimizers += 1
            prev_argcummin = combined_argcummin
        # If the same minimizer position was observed before,
        # the duplicate is simply ignored

    return (
        np.asarray(mininizer_pos)[:n_minimizers],
        np.asarray(minimizers)[:n_minimizers]
    )


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _chunk_wise_forward_argcummin(int64[:] values, uint32 chunk_size):
    """
    Argument of the cumulative minimum.
    """
    cdef uint32 seq_i

    cdef uint32 current_min_i = 0
    cdef int64 current_min, current_val
    cdef uint32[:] min_pos = np.empty(values.shape[0], dtype=np.uint32)
    
    # Any actual value will be smaller than this placeholder
    current_min = MAX_INT_64
    for seq_i in range(values.shape[0]):
        if seq_i % chunk_size == 0:
            # New chunk begins
            current_min = MAX_INT_64
        current_val = values[seq_i]
        if current_val < current_min:
            current_min_i = seq_i
            current_min = current_val
        min_pos[seq_i] = current_min_i
    
    return min_pos


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _chunk_wise_reverse_argcummin(int64[:] values, uint32 chunk_size):
    """
    The same as above but starting from the other end and iterating
    backwards.
    Separation into two functions leads to code duplication.
    However, single implemention with reversed `values` as input
    has some disadvantages:

    - Indices must be transformed so that they point to the
      non-reversed `values`
    - There are issues in selecting the leftmost argument
    - An offset is necessary to ensure alignment of chunks with forward
      pass
    
    Hence, a separate 'reverse' variant of the function was implemented.
    """
    cdef uint32 seq_i

    cdef uint32 current_min_i = 0
    cdef int64 current_min, current_val
    cdef uint32[:] min_pos = np.empty(values.shape[0], dtype=np.uint32)
    
    current_min = MAX_INT_64
    for seq_i in reversed(range(values.shape[0])):
        # The chunk beginning is a small difference to forward
        # implementation, as it begins on the left of the chunk border
        if seq_i % chunk_size == chunk_size - 1:
            current_min = MAX_INT_64
        current_val = values[seq_i]
        # The '<=' is a small difference to forward implementation
        # to enure the loftmost argument is selected
        if current_val <= current_min:
            current_min_i = seq_i
            current_min = current_val
        min_pos[seq_i] = current_min_i
    
    return min_pos