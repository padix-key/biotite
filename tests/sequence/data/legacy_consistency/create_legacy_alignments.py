"""
Create legacy alignment files from the Cython implementations.

This script must be run while the Cython implementations of
``align_optimal()``, ``align_banded()`` and ``align_local_gapped()`` are still
active.  It will refuse to run if any of these functions has already been
replaced by a Rust implementation.
"""

import itertools
import json
from pathlib import Path

import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta

BAND_WIDTH = 100
GAP_PENALTY = (-10, -1)
THRESHOLD = 100

DATA_DIR = Path(__file__).resolve().parent
CAS9_PATH = DATA_DIR.parent / "cas9.fasta"


def _assert_cython(func):
    """
    Raise if `func` is not a Cython function.
    """
    type_name = type(func).__name__
    if type_name != "cython_function_or_method":
        raise RuntimeError(
            "This script must be run with the legacy Cython implementation."
        )


def _middle_aligned_position(alignment):
    """
    Return the ``(seq1_pos, seq2_pos)`` tuple from the middle of the aligned
    (non-gap) positions in the alignment.
    """
    trace = alignment.trace
    aligned_mask = (trace != -1).all(axis=1)
    aligned_positions = trace[aligned_mask]
    mid = aligned_positions[len(aligned_positions) // 2]
    return mid[0].item(), mid[1].item()


def _write_alignment(alignment, seq_names, path):
    """
    Write the alignment to a FASTA file.
    """
    fasta_file = fasta.FastaFile()
    fasta.set_alignment(fasta_file, alignment, seq_names)
    fasta_file.write(path)


def main():
    _assert_cython(align.align_optimal)
    _assert_cython(align.align_banded)
    _assert_cython(align.align_local_gapped)

    fasta_file = fasta.FastaFile.read(CAS9_PATH)
    headers = list(fasta_file.keys())
    sequences = [seq.ProteinSequence(sequence) for sequence in fasta_file.values()]
    matrix = align.SubstitutionMatrix.std_protein_matrix()

    params = {
        "align_optimal": {},
        "align_banded": {},
        "align_local_gapped": {},
    }

    for i, j in itertools.combinations(range(len(sequences)), 2):
        seq1, seq2 = sequences[i], sequences[j]
        seq_names = [headers[i], headers[j]]
        pair_key = f"{i}_{j}"

        # --- align_optimal ---
        optimal_alignment = align.align_optimal(
            seq1, seq2, matrix, gap_penalty=GAP_PENALTY, max_number=1
        )[0]
        _write_alignment(
            optimal_alignment, seq_names, DATA_DIR / f"align_optimal_{i}_{j}.fasta"
        )
        params["align_optimal"][pair_key] = {
            "gap_penalty": list(GAP_PENALTY),
        }

        # --- align_banded ---
        mid_pos1, mid_pos2 = _middle_aligned_position(optimal_alignment)
        diag = mid_pos2 - mid_pos1
        band = (diag - BAND_WIDTH, diag + BAND_WIDTH)
        banded_alignment = align.align_banded(
            seq1, seq2, matrix, band=band, gap_penalty=GAP_PENALTY, max_number=1
        )[0]
        _write_alignment(
            banded_alignment, seq_names,
            DATA_DIR / f"align_banded_{i}_{j}.fasta"
        )
        params["align_banded"][pair_key] = {
            "gap_penalty": list(GAP_PENALTY),
            "band": list(band),
        }

        # --- align_local_gapped ---
        seed = (mid_pos1, mid_pos2)
        local_alignment = align.align_local_gapped(
            seq1, seq2, matrix, seed=seed, threshold=THRESHOLD,
            gap_penalty=GAP_PENALTY, max_number=1
        )[0]
        _write_alignment(
            local_alignment, seq_names,
            DATA_DIR / f"align_local_{i}_{j}.fasta"
        )
        params["align_local_gapped"][pair_key] = {
            "gap_penalty": list(GAP_PENALTY),
            "seed": list(seed),
            "threshold": THRESHOLD,
        }

    with open(DATA_DIR / "params.json", "w") as f:
        json.dump(params, f, indent=2)


if __name__ == "__main__":
    main()
