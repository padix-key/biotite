# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import gc
import subprocess
from pathlib import Path
import numpy as np
import pytest
from biotite.application_v2.localapp import LocalProcessFuture
from biotite.application_v2.mmseqs import (
    AlignmentDatabase,
    AlignmentFormatMode,
    DatabaseType,
    FoldseekApp,
    MMseqsApp,
    MSADatabase,
    MSAFormatMode,
    SequenceDatabase,
    create_database_from_protein_sequences,
)
from biotite.sequence.align.matrix import SubstitutionMatrix
from biotite.sequence.seqtypes import ProteinSequence


@pytest.mark.parametrize("seed", range(5))
def test_create_db_command_properties(seed):
    """
    Command construction retains every input and formats generated options.
    """
    rng = np.random.default_rng(seed)
    input_count = int(rng.integers(1, 6))
    threads = int(rng.integers(1, 129))
    database_types = list(DatabaseType)
    database_type = database_types[int(rng.integers(len(database_types)))]
    app = MMseqsApp("true")
    inputs = [Path(f"input_{index}.fasta") for index in range(input_count)]

    future = app.create_db(inputs, threads=threads, db_type=database_type)

    assert isinstance(future, LocalProcessFuture)
    assert future.command.startswith(f"{app.path} createdb ")
    assert f"--threads {threads}" in future.command
    assert f"--dbtype {database_type.value}" in future.command
    assert all(str(path) in future.command for path in inputs)
    assert isinstance(future.result(), SequenceDatabase)


@pytest.mark.parametrize("seed", range(5))
def test_database_value_properties(seed):
    """
    Every generated database suffix is passed verbatim to the command line.
    """
    rng = np.random.default_rng(seed)
    length = int(rng.integers(0, 9))
    suffix = "".join(rng.choice(list("_abcdefghijklmnopqrstuvwxyz"), size=length))
    app = MMseqsApp("true")
    database = SequenceDatabase(app, suffix=suffix)
    fasta_path = database.path / "output.fasta"

    future = app.convert_to_fasta(database, fasta_path=fasta_path)

    assert str(database.name) in future.command
    assert future.result() == fasta_path


def test_conversion_to_temporary_file():
    """
    Omitting a conversion path returns a readable temporary binary file.
    """
    app = MMseqsApp("true")
    sequence_db = SequenceDatabase(app)
    alignment_db = AlignmentDatabase(app, sequence_db, sequence_db)

    futures = [
        app.convert_to_fasta(sequence_db),
        app.convert_alignments(alignment_db),
    ]
    for future in futures:
        output_file = future.result()
        output_path = Path(output_file.name)

        assert output_file.readable()
        assert "b" in output_file.mode
        assert str(output_path) in future.command
        assert output_path.is_file()

        output_file.close()
        assert not output_path.exists()


def test_create_db_accepts_text_file(tmp_path):
    """
    An open text file can be passed directly as a ``createdb`` input.
    """
    fasta_path = tmp_path / "sequences.fasta"
    fasta_path.write_text(">sequence\nBIQTITE\n")
    app = MMseqsApp("true")

    with fasta_path.open() as fasta_file:
        future = app.create_db(fasta_file)
        assert fasta_file.name in future.command
        assert isinstance(future.result(), SequenceDatabase)


@pytest.mark.parametrize(
    "command_name",
    ["databases", "create_db", "create_sub_db", "search", "result_to_msa"],
)
def test_database_command_rejects_output_path(command_name, tmp_path):
    """
    Database results are persisted with :meth:`Database.move`, not an output option.
    """
    app = MMseqsApp("true")
    sequence_db = SequenceDatabase(app)
    alignment_db = AlignmentDatabase(app, sequence_db, sequence_db)
    args = {
        "databases": ("UniRef50",),
        "create_db": ("input.fasta",),
        "create_sub_db": ("keys.tsv", sequence_db),
        "search": (sequence_db, sequence_db),
        "result_to_msa": (alignment_db,),
    }[command_name]

    with pytest.raises(ValueError, match="not an allowed"):
        getattr(app, command_name)(*args, output_path=tmp_path)


def test_disallowed_option():
    """
    Unknown command options are rejected before a process is launched.
    """
    with pytest.raises(ValueError, match="not an allowed"):
        MMseqsApp("true").create_db("input.fasta", invalid_option=True)


@pytest.mark.parametrize(
    ["command_name", "option"],
    [
        ("databases", {"tsv": True}),
        ("create_db", {"createdb_mode": 1}),
        ("create_sub_db", {"subdb_mode": 1}),
        ("search", {"alignment_output_mode": 5}),
        ("search", {"local_tmp": "/tmp"}),
        ("convert_alignments", {"db_output": True}),
    ],
)
def test_result_changing_option_is_disallowed(command_name, option):
    """
    Options that invalidate a command's result abstraction are rejected.
    """
    app = MMseqsApp("true")
    sequence_db = SequenceDatabase(app)
    alignment_db = AlignmentDatabase(app, sequence_db, sequence_db)
    args = {
        "databases": ("UniRef50",),
        "create_db": ("input.fasta",),
        "create_sub_db": ("keys.tsv", sequence_db),
        "search": (sequence_db, sequence_db),
        "convert_alignments": (alignment_db, "alignments.tsv"),
        "result_to_msa": (alignment_db,),
    }[command_name]

    with pytest.raises(ValueError, match="not an allowed"):
        getattr(app, command_name)(*args, **option)


def test_enum_options_are_typed_parameters():
    """
    Enum-backed options are formatted by their respective command methods.
    """
    app = MMseqsApp("true")
    sequence_db = SequenceDatabase(app)
    alignment_db = AlignmentDatabase(app, sequence_db, sequence_db)

    create_future = app.create_db("input.fasta", db_type=DatabaseType.AMINO_ACID)
    convert_future = app.convert_alignments(
        alignment_db,
        "alignments.tsv",
        format_mode=AlignmentFormatMode.BLAST_TABLE_WITH_HEADERS,
    )
    msa_future = app.result_to_msa(
        alignment_db,
        msa_format_mode=MSAFormatMode.A3M,
    )

    assert f"--dbtype {DatabaseType.AMINO_ACID.value}" in create_future.command
    assert (
        f"--format-mode {AlignmentFormatMode.BLAST_TABLE_WITH_HEADERS.value}"
        in convert_future.command
    )
    assert f"--msa-format-mode {MSAFormatMode.A3M.value}" in msa_future.command

    create_future.result()
    convert_future.result()
    msa_future.result()


def test_flat_msa_format_is_disallowed():
    """
    A flat-file MSA format cannot be represented by :class:`MSADatabase`.
    """
    app = MMseqsApp("true")
    sequence_db = SequenceDatabase(app)
    alignment_db = AlignmentDatabase(app, sequence_db, sequence_db)

    with pytest.raises(ValueError, match="flat file"):
        app.result_to_msa(alignment_db, msa_format_mode=MSAFormatMode.STOCKHOLM)


def test_process_error_is_deferred():
    """
    A launched command reports process errors through its future.
    """
    future = MMseqsApp("false").create_db("input.fasta")
    database_dir = Path(future.command.split()[-1]).parent
    assert isinstance(future, LocalProcessFuture)
    with pytest.raises(subprocess.SubprocessError):
        future.result()
    assert database_dir.is_dir()

    del future
    gc.collect()
    assert not database_dir.exists()


def test_incompatible_database():
    """
    A database cannot cross the MMseqs2/Foldseek application boundary.
    """
    database = SequenceDatabase(FoldseekApp("true"))
    with pytest.raises(ValueError, match="not compatible"):
        MMseqsApp("true").convert_to_fasta(database, "output.fasta")


def test_tmp_dir_is_instance_specific(tmp_path):
    """
    Each application instance has its own configurable temporary root directory.
    """
    tmp_dir = tmp_path / "mmseqs"

    mmseqs_app = MMseqsApp("true")
    foldseek_app = FoldseekApp("true")
    mmseqs_app.tmp_dir = tmp_dir

    assert mmseqs_app.tmp_dir == tmp_dir
    assert foldseek_app.tmp_dir != tmp_dir
    assert tmp_dir.is_dir()


@pytest.mark.parametrize(
    ["app_class", "prefix"],
    [(MMseqsApp, "mmseqs_"), (FoldseekApp, "foldseek_")],
)
def test_search_uses_isolated_temporary_work_directories(app_class, prefix, tmp_path):
    """
    Concurrent searches use distinct work directories below the configured root.
    """
    app = app_class("true")
    app.tmp_dir = tmp_path
    database = SequenceDatabase(app)

    futures = [app.search(database, database) for _ in range(2)]
    work_dirs = [Path(future.command.split()[-1]) for future in futures]

    assert all(work_dir.parent == tmp_path for work_dir in work_dirs)
    assert all(work_dir.name.startswith(prefix) for work_dir in work_dirs)
    assert len(set(work_dirs)) == len(work_dirs)
    assert all(work_dir.is_dir() for work_dir in work_dirs)

    for future in futures:
        future.result()
    assert all(work_dir.is_dir() for work_dir in work_dirs)

    del future
    del futures
    gc.collect()
    assert all(not work_dir.exists() for work_dir in work_dirs)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("affine", [False, True])
def test_search_formats_matrix_and_gap_penalty(seed, affine):
    """
    A custom matrix and affine penalty are converted into MMseqs2 options.
    """
    rng = np.random.default_rng(seed)
    matrix_names = ["BLOSUM45", "BLOSUM62", "PAM30", "PAM250"]
    matrix_name = matrix_names[int(rng.integers(len(matrix_names)))]
    gap_open = -int(rng.integers(5, 20))
    gap_extend = -int(rng.integers(1, 5))
    gap_penalty = (gap_open, gap_extend) if affine else gap_open
    expected_extend = gap_extend if affine else gap_open
    app = MMseqsApp("true")
    database = SequenceDatabase(app)

    future = app.search(
        database,
        database,
        matrix=SubstitutionMatrix(
            ProteinSequence.alphabet,
            ProteinSequence.alphabet,
            matrix_name,
        ),
        gap_penalty=gap_penalty,
    )

    assert f"--gap-open {-gap_open}" in future.command
    assert f"--gap-extend {-expected_extend}" in future.command
    assert "--sub-mat" in future.command
    matrix_path = Path(future.command.split("--sub-mat ", 1)[1].split()[0])
    matrix_text = matrix_path.read_text()
    assert "# Background (precomputed optional):" not in matrix_text
    assert "# Lambda     (precomputed optional):" not in matrix_text
    matrix_lines = [
        line for line in matrix_text.splitlines() if not line.startswith("#")
    ]
    assert matrix_lines[0].split()[-1] == "X"
    assert matrix_lines[-1].split()[0] == "X"
    assert isinstance(future.result(), AlignmentDatabase)


@pytest.mark.parametrize(
    "option",
    [
        {"matrix": SubstitutionMatrix.std_3di_matrix()},
        {"gap_penalty": (-10, -1)},
    ],
)
def test_foldseek_rejects_custom_scoring(option):
    """
    Foldseek does not expose custom matrices or gap penalties.
    """
    app = FoldseekApp("true")
    database = SequenceDatabase(app)

    with pytest.raises(ValueError, match="not an allowed"):
        app.search(database, database, **option)


def test_mmseqs_command_pipeline(mmseqs_app, tmp_path):
    """
    All local MMseqs2 database commands compose through their future results.
    """
    sequences = {
        "first": "CGESCVYIPCTVTALLGCSCKDKVCYKNSLAVN",
        "second": "CGESCVYIPCTVTALLGCSCKDKVCYKNSLAVM",
    }
    sequence_db = create_database_from_protein_sequences(mmseqs_app, sequences)

    fasta_path = tmp_path / "sequences.fasta"
    assert mmseqs_app.convert_to_fasta(sequence_db, fasta_path).result() == fasta_path
    assert fasta_path.is_file()

    keys_path = tmp_path / "keys.txt"
    keys_path.write_text("0\n")
    subset_db = mmseqs_app.create_sub_db(keys_path, sequence_db).result()
    assert isinstance(subset_db, SequenceDatabase)

    alignment_db = mmseqs_app.search(
        sequence_db,
        sequence_db,
        matrix=SubstitutionMatrix.std_protein_matrix(),
        gap_penalty=(-11, -1),
        a=True,
        k=5,
        threads=1,
    ).result()
    assert isinstance(alignment_db, AlignmentDatabase)

    alignment_path = tmp_path / "alignments.tsv"
    assert (
        mmseqs_app.convert_alignments(
            alignment_db,
            alignment_path,
            format_output="query,target",
            threads=1,
        ).result()
        == alignment_path
    )
    assert alignment_path.is_file()

    msa_db = mmseqs_app.result_to_msa(alignment_db, threads=1).result()
    assert isinstance(msa_db, MSADatabase)

    linked_path = tmp_path / "linked_database"
    assert mmseqs_app.link_db(sequence_db.name, linked_path).result() == linked_path
    assert linked_path.is_symlink()


def test_foldseek_command_pipeline(foldseek_app, structure_path, tmp_path):
    """
    Foldseek uses the same future and database abstractions for structure inputs.
    """
    sequence_db = foldseek_app.create_db(structure_path, threads=1).result()
    assert isinstance(sequence_db, SequenceDatabase)

    fasta_path = tmp_path / "structure.fasta"
    foldseek_app.convert_to_fasta(sequence_db, fasta_path).result()
    assert fasta_path.is_file()

    alignment_db = foldseek_app.search(
        sequence_db,
        sequence_db,
        a=True,
        threads=1,
    ).result()
    assert isinstance(alignment_db, AlignmentDatabase)

    alignment_path = tmp_path / "structure.tsv"
    foldseek_app.convert_alignments(
        alignment_db,
        alignment_path,
        format_output="query,target",
        threads=1,
    ).result()
    assert alignment_path.is_file()
