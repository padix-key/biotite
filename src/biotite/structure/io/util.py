"""
Common functions used by a number of subpackages.
"""

__name__ = "biotite.structure.io"
__author__ = "Patrick Kunzmann"
__all__ = [
    "number_of_integer_digits",
    "convert_unicode_to_uint32",
    "convert_uint32_to_unicode",
]

from typing import Any
import numpy as np
from biotite.structure.box import vectors_from_unitcell
from biotite.typing import XYZ, N, NDArray1, NDArray2, NDArray3


def number_of_integer_digits(values: NDArray1[Any, np.number]) -> int:
    """
    Get the maximum number of characters needed to represent the
    pre-decimal positions of the given numeric values.

    Parameters
    ----------
    values : ndarray, dtype=number
        The values to be checked.

    Returns
    -------
    n_digits : int
        The maximum number of characters needed to represent the
        pre-decimal positions of the given numeric values.
    """
    if len(values) == 0:
        return 0
    values = values.astype(int, copy=False)
    n_digits = 0
    n_digits = max(n_digits, len(str(np.min(values))))
    n_digits = max(n_digits, len(str(np.max(values))))
    return n_digits


def convert_unicode_to_uint32(
    array: NDArray1[N, np.str_],
) -> NDArray2[N, Any, np.uint32]:
    """
    Convert a unicode string array into a 2D uint32 array.

    The second dimension corresponds to the character position within a
    string.

    Parameters
    ----------
    array : ndarray, dtype=unicode
        The unicode string array to be converted.

    Returns
    -------
    array : ndarray, dtype=uint32
        The converted unicode string array.
        Each element along the second dimension corresponds to an UTF-32 symbol.
    """
    dtype = array.dtype
    if not np.issubdtype(dtype, np.str_):
        raise TypeError("Expected unicode string array")
    length = array.shape[0]
    n_char = dtype.itemsize // 4
    return np.frombuffer(array, dtype=np.uint32).reshape(length, n_char)  # pyright: ignore[reportReturnType]


def convert_uint32_to_unicode(
    array: NDArray2[N, Any, np.uint32],
) -> NDArray1[N, np.str_]:
    """
    Convert a 2D uint32 array into a unicode string array.

    This is the inverse of :func:`convert_unicode_to_uint32`.

    Parameters
    ----------
    array : ndarray, dtype=uint32, shape=(n, m)
        The uint32 array to be converted.
        Each element along the second dimension corresponds to an UTF-32 symbol.

    Returns
    -------
    array : ndarray, dtype=unicode
        The converted unicode string array.
    """
    if array.ndim != 2:
        raise ValueError("Expected 2D array")
    length, n_char = array.shape
    return np.frombuffer(
        array.astype(np.uint32, copy=False).tobytes(), dtype=f"<U{n_char}"
    )


def vectors_from_unitcells(
    cell_lengths: NDArray2[N, XYZ, np.floating],
    cell_angles: NDArray2[N, XYZ, np.floating],
) -> NDArray3[N, XYZ, XYZ, np.float32]:
    """
    Vectorized version of :func:`biotite.structure.vectors_from_unitcell()`
    for multiple unit cells, e.g. for each frame of a trajectory.

    Parameters
    ----------
    cell_lengths : ndarray, shape=(n,3), dtype=float
        The lengths of the three box vectors for each unit cell.
    cell_angles : ndarray, shape=(n,3), dtype=float
        The angles *alpha*, *beta*, *gamma* for each unit cell in radians.

    Returns
    -------
    box : ndarray, shape=(n,3,3), dtype=float32
        The box vectors for each unit cell.
    """
    if len(cell_lengths) == 1:
        # Avoid the vectorization overhead for single frames, e.g. in 'read_iter()'
        return vectors_from_unitcell(*cell_lengths[0], *cell_angles[0])[np.newaxis]
    len_a, len_b, len_c = cell_lengths.T
    alpha, beta, gamma = cell_angles.T
    c_x = len_c * np.cos(beta)
    c_y = len_c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    box = np.zeros((len(cell_lengths), 3, 3), dtype=np.float32)
    box[:, 0, 0] = len_a
    box[:, 1, 0] = len_b * np.cos(gamma)
    box[:, 1, 1] = len_b * np.sin(gamma)
    box[:, 2, 0] = c_x
    box[:, 2, 1] = c_y
    box[:, 2, 2] = np.sqrt(len_c * len_c - c_x * c_x - c_y * c_y)

    # Fix numerical errors, as values, that are actually 0,
    # might not be calculated as such
    tol = 1e-4 * (len_a + len_b + len_c)
    box[np.abs(box) < tol[:, np.newaxis, np.newaxis]] = 0
    return box
