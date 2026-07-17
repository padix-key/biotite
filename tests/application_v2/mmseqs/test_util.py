# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import numpy as np
import pytest
import biotite.sequence as seq
import biotite.structure as struc
import biotite.structure.alphabet as strucalph
import biotite.structure.io.pdbx as pdbx
from biotite.application_v2.mmseqs import (
    create_database_from_protein_sequences,
    get_alignment_table,
    get_alignments,
    get_sequences_from_database,
)

PROTEIN_ALPHABET = np.array(list("ACDEFGHIKLMNPQRSTVWY"))


def _random_protein_sequence(rng):
    length = int(rng.integers(35, 61))
    return "".join(rng.choice(PROTEIN_ALPHABET, size=length))


def _random_sequence_mapping(rng):
    count = int(rng.integers(1, 5))
    return {
        f"sequence_{index}": _random_protein_sequence(rng) for index in range(count)
    }


@pytest.mark.parametrize("seed", range(5))
def test_sequence_roundtrip_properties(mmseqs_app, seed):
    """
    Arbitrary valid protein mappings survive a database roundtrip unchanged.
    """
    sequences = _random_sequence_mapping(np.random.default_rng(seed))
    database = create_database_from_protein_sequences(mmseqs_app, sequences)

    assert get_sequences_from_database(mmseqs_app, database) == sequences


@pytest.mark.parametrize("seed", range(5))
def test_alignment_result_properties(mmseqs_app, seed):
    """
    Self-search tables are rectangular and reconstruct identity alignments.
    """
    sequence = _random_protein_sequence(np.random.default_rng(seed))
    database = create_database_from_protein_sequences(
        mmseqs_app, {"query": seq.ProteinSequence(sequence)}
    )
    alignment_db = mmseqs_app.search(
        database,
        database,
        a=True,
        k=5,
        e=10,
        threads=1,
    ).result()

    columns = ["query", "target", "qstart", "tstart", "cigar", "evalue"]
    table = get_alignment_table(mmseqs_app, alignment_db, columns)
    assert list(table) == columns
    assert len({len(values) for values in table.values()}) == 1
    assert set(table["query"]) == {"query"}
    assert set(table["target"]) == {"query"}

    alignments = get_alignments(mmseqs_app, alignment_db)
    assert list(alignments) == [("query", "query")]
    alignment = alignments[("query", "query")]
    assert alignment.sequences[0] == alignment.sequences[1]
    assert str(alignment).splitlines()[0] == str(alignment).splitlines()[1]


@pytest.mark.parametrize("as_3di", [False, True])
def test_foldseek_sequence_properties(foldseek_app, structure_path, as_3di):
    """
    Foldseek extraction agrees with Biotite for amino-acid and 3Di sequences.
    """
    atoms = pdbx.get_structure(pdbx.CIFFile.read(structure_path), model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    if as_3di:
        reference_sequences = strucalph.to_3di(atoms)
    else:
        reference_sequences = struc.to_sequence(atoms)
    reference = {
        f"1aki_{atoms.chain_id[start]}": str(sequence).upper()
        for sequence, start in zip(*reference_sequences)
    }

    database = foldseek_app.create_db(
        structure_path, chain_name_mode=1, threads=1
    ).result()
    sequences = get_sequences_from_database(foldseek_app, database, as_3di=as_3di)

    assert sequences == reference


@pytest.mark.parametrize("columns", [[], ["query", "query"]])
def test_invalid_alignment_columns(mmseqs_app, columns):
    """
    Empty or duplicate table schemas are rejected without launching a process.
    """
    with pytest.raises(ValueError):
        get_alignment_table(mmseqs_app, None, columns)
