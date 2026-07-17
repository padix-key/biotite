# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from __future__ import annotations

__name__ = "biotite.application_v2.mmseqs"
__author__ = "Patrick Kunzmann"
__all__ = [
    "link_to_ss_database",
    "create_database_from_protein_sequences",
    "get_sequences_from_database",
    "get_alignment_table",
    "get_alignments",
]

import csv
from collections.abc import Mapping
from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import IO, TypeVar
from biotite.application_v2.mmseqs.app import (
    AlignmentFormatMode,
    DatabaseType,
    FoldseekApp,
    MMseqsApp,
    MMseqsLikeApp,
)
from biotite.application_v2.mmseqs.database import (
    AlignmentDatabase,
    SequenceDatabase,
)
from biotite.sequence import ProteinSequence, Sequence
from biotite.sequence.align import Alignment, read_alignment_from_cigar
from biotite.sequence.io.fasta import FastaFile
from biotite.structure.alphabet import I3DSequence

T = TypeVar("T", bound="MMseqsLikeApp")


def link_to_ss_database(
    app: T, database: SequenceDatabase[T]
) -> SequenceDatabase[T]:
    """
    Link a Foldseek database to its 3Di sub-database.

    Parameters
    ----------
    app : MMseqsApp or FoldseekApp
        The application to use.
    database : SequenceDatabase
        The Foldseek database containing the 3Di sub-database.

    Returns
    -------
    linked_database : SequenceDatabase
        The linked 3Di database.
    """

    app.link_db(
        Path(f"{database.name}_h"),
        Path(f"{database.name}_ss_h"),
    ).result()
    ss_database = SequenceDatabase(app, database.path, suffix="_ss")
    ss_database.set_parent(database)
    return ss_database


def create_database_from_protein_sequences(
    app: T,
    sequences: Mapping[str, str | Sequence],
    prostt5_weights: SequenceDatabase[T] | None = None,
) -> SequenceDatabase[T]:
    """
    Create a database from protein sequences.

    Foldseek predicts corresponding 3Di sequences using ProstT5.

    Parameters
    ----------
    app : MMseqsApp or FoldseekApp
        The application to use.
    sequences : mapping of str to str or Sequence
        Sequence identifiers mapped to protein sequences.
    prostt5_weights : SequenceDatabase, optional
        ProstT5 weights for Foldseek. If omitted, they are downloaded.

    Returns
    -------
    database : SequenceDatabase
        The created database.
    """
    with NamedTemporaryFile("w+", suffix=".fasta") as temp_file:
        FastaFile.write_iter(
            temp_file,
            ((identifier, str(sequence)) for identifier, sequence in sequences.items()),
            chars_per_line=1_000_000,
        )
        temp_file.flush()

        match app:
            case MMseqsApp():
                return app.create_db(
                    [temp_file],
                    db_type=DatabaseType.AMINO_ACID,
                ).result()
            case FoldseekApp():
                if prostt5_weights is None:
                    prostt5_weights = app.databases("ProstT5").result()
                return app.create_db(
                    [temp_file], prostt5_model=prostt5_weights
                ).result()
            case _:
                raise TypeError(f"Unsupported application: '{type(app).__name__}'")


def get_sequences_from_database(
    app: T,
    sequence_db: SequenceDatabase[T],
    as_3di: bool = False,
) -> dict[str, str]:
    """
    Extract the sequences from a database.

    Parameters
    ----------
    app : MMseqsApp or FoldseekApp
        The application to use.
    sequence_db : SequenceDatabase
        The database to read.
    as_3di : bool, optional
        If true, extract Foldseek 3Di instead of amino-acid sequences.

    Returns
    -------
    sequences : dict of str to str
        Sequence identifiers mapped to sequences.
    """

    def read_fasta(file: IO[bytes]) -> dict[str, str]:
        # The application writes through a separate file descriptor, so the
        # retained temporary file's read buffer may still cache the original
        # empty file. Reopening its path observes the application output.
        with open(file.name, "rb") as input_file:
            content = input_file.read().decode()
        with StringIO(content) as text:
            return {
                header.partition(" ")[0]: sequence
                for header, sequence in FastaFile.read_iter(text)
            }

    if as_3di:
        sequence_db = link_to_ss_database(app, sequence_db)
    file = app.convert_to_fasta(sequence_db).result()
    try:
        return read_fasta(file)
    finally:
        file.close()


def get_alignment_table(
    app: T,
    alignment_db: AlignmentDatabase[T],
    columns: list[str],
) -> dict[str, list[str]]:
    """
    Export selected alignment columns.

    Parameters
    ----------
    app : MMseqsApp or FoldseekApp
        The application to use.
    alignment_db : AlignmentDatabase
        The alignment database to read.
    columns : list of str
        The ``convertalis --format-output`` columns to include.

    Returns
    -------
    table : dict of str to list of str
        Columns mapped to their raw string values.
    """
    if len(columns) == 0:
        raise ValueError("At least one alignment column is required")
    if len(set(columns)) != len(columns):
        raise ValueError("Alignment columns must be unique")

    def read_table(file: IO[bytes]) -> dict[str, list[str]]:
        table = {column: [] for column in columns}
        # See `read_fasta()` for why the temporary path is reopened here.
        with open(file.name, "rb") as input_file:
            content = input_file.read().decode()
        with StringIO(content, newline="") as text:
            for row in csv.DictReader(text, delimiter="\t"):
                for column in columns:
                    table[column].append(row[column])
        return table

    file = app.convert_alignments(
        alignment_db,
        format_mode=AlignmentFormatMode.BLAST_TABLE_WITH_HEADERS,
        format_output=",".join(columns),
    ).result()
    try:
        return read_table(file)
    finally:
        file.close()


def get_alignments(
    app: T,
    alignment_db: AlignmentDatabase[T],
    as_3di: bool = False,
) -> dict[tuple[str, str], Alignment]:
    """
    Reconstruct Biotite alignments from an alignment database.

    The search creating `alignment_db` must enable backtracing with ``a=True``.

    Parameters
    ----------
    app : MMseqsApp or FoldseekApp
        The application to use.
    alignment_db : AlignmentDatabase
        The alignment database to read.
    as_3di : bool, optional
        If true, reconstruct alignments from Foldseek 3Di sequences.

    Returns
    -------
    alignments : dict
        ``(query_id, target_id)`` pairs mapped to alignments. Entries are ordered
        by increasing E-value.
    """
    columns = ["query", "target", "qstart", "tstart", "cigar", "evalue"]
    table = get_alignment_table(app, alignment_db, columns)
    query_sequences = get_sequences_from_database(app, alignment_db.query_db, as_3di)
    target_sequences = get_sequences_from_database(app, alignment_db.target_db, as_3di)
    sequence_type = I3DSequence if as_3di else ProteinSequence
    row_indices = sorted(
        range(len(table["query"])),
        key=lambda index: float(table["evalue"][index]),
    )
    alignments = {}
    for index in row_indices:
        query_id = table["query"][index]
        target_id = table["target"][index]
        query = sequence_type(query_sequences[query_id])
        target = sequence_type(target_sequences[target_id])
        query_start = int(table["qstart"][index]) - 1
        target_start = int(table["tstart"][index]) - 1
        alignment = read_alignment_from_cigar(table["cigar"][index], 0, target, query)[
            :, ::-1
        ]
        alignment.trace[alignment.trace[:, 0] != -1, 0] += query_start
        alignment.trace[alignment.trace[:, 1] != -1, 1] += target_start
        alignments[(query_id, target_id)] = alignment
    return alignments
