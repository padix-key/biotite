# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from __future__ import annotations

__name__ = "biotite.application_v2.mmseqs"
__author__ = "Patrick Kunzmann"
__all__ = [
    "DatabaseType",
    "AlignmentFormatMode",
    "MSAFormatMode",
    "MMseqsApp",
    "FoldseekApp",
]

from collections.abc import Iterable
from enum import IntEnum
from io import IOBase
from os import PathLike
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory, gettempdir
from typing import IO, Any, ClassVar, overload
import numpy as np
from biotite.application_v2.localapp import (
    CLIArgument,
    CLIOption,
    CommandSetup,
    LocalApp,
    LocalProcessFuture,
    command,
)
from biotite.application_v2.mmseqs.database import (
    AlignmentDatabase,
    Database,
    MSADatabase,
    SequenceDatabase,
)
from biotite.application_v2.msa import resolve_gap_penalty
from biotite.sequence.align.matrix import SubstitutionMatrix
from biotite.sequence.seqtypes import ProteinSequence


class DatabaseType(IntEnum):
    """
    Sequence types accepted by the MMseqs2 ``createdb`` command.
    """

    AUTO = 0
    AMINO_ACID = 1
    NUCLEOTIDE = 2


class AlignmentFormatMode(IntEnum):
    """
    Output modes accepted by the ``convertalis`` command.
    """

    BLAST_TABLE = 0
    SAM = 1
    BLAST_TABLE_WITH_LENGTH = 2
    HTML = 3
    BLAST_TABLE_WITH_HEADERS = 4
    CALPHA_PDB = 5


class MSAFormatMode(IntEnum):
    """
    Output modes accepted by the ``result2msa`` command.
    """

    COMPACT_A3M = 0
    COMPACT_A3M_WITH_CONSENSUS = 1
    FASTA = 2
    FASTA_WITH_INFO = 3
    STOCKHOLM = 4
    A3M = 5
    A3M_WITH_INFO = 6


_COMMON_OPTIONS = ["threads", "compressed", "v"]
_DATABASES_OPTIONS = [
    "force_reuse",
    "remove_tmp_files",
    *_COMMON_OPTIONS,
]
_CREATE_DB_OPTIONS = [
    "shuffle",
    "id_offset",
    "write_lookup",
    "mask",
    "mask_prob",
    "mask_lower_case",
    "mask_n_repeat",
    "mask_bfactor_threshold",
    "input_format",
    "file_include",
    "file_exclude",
    "prostt5_model",
    "chain_name_mode",
    "write_mapping",
    "coord_store_mode",
    "db_extraction_mode",
    "distance_threshold",
    "gpu",
    *_COMMON_OPTIONS,
]
_CREATE_SUB_DB_OPTIONS = ["id_mode", "v"]
_CONVERT_TO_FASTA_OPTIONS = ["use_header_file", "v"]
_SEARCH_OPTIONS = [
    "comp_bias_corr",
    "comp_bias_corr_scale",
    "add_self_matches",
    "seed_sub_mat",
    "s",
    "k",
    "target_search_mode",
    "k_score",
    "alph_size",
    "max_seqs",
    "split",
    "split_mode",
    "split_memory_limit",
    "diag_score",
    "exact_kmer_matching",
    "mask",
    "mask_prob",
    "mask_lower_case",
    "mask_n_repeat",
    "min_ungapped_score",
    "spaced_kmer_mode",
    "spaced_kmer_pattern",
    "disk_space_limit",
    "a",
    "alignment_mode",
    "wrapped_scoring",
    "e",
    "min_seq_id",
    "min_aln_len",
    "seq_id_mode",
    "alt_ali",
    "c",
    "cov_mode",
    "max_rejected",
    "max_accept",
    "score_bias",
    "realign",
    "realign_score_bias",
    "realign_max_seqs",
    "corr_score_weight",
    "zdrop",
    "exhaustive_search_filter",
    "pca",
    "pcb",
    "mask_profile",
    "e_profile",
    "wg",
    "filter_msa",
    "filter_min_enable",
    "max_seq_id",
    "qid",
    "qsc",
    "cov",
    "diff",
    "pseudo_cnt_mode",
    "profile_output_mode",
    "num_iterations",
    "exhaustive_search",
    "lca_search",
    "taxon_list",
    "prefilter_mode",
    "rescore_mode",
    "allow_deletion",
    "min_length",
    "max_length",
    "max_gaps",
    "contig_start_mode",
    "contig_end_mode",
    "orf_start_mode",
    "forward_frames",
    "reverse_frames",
    "translation_table",
    "translate",
    "use_all_table_starts",
    "id_offset",
    "sequence_overlap",
    "sequence_split_mode",
    "headers_split_mode",
    "search_type",
    "start_sens",
    "sens_steps",
    "translation_mode",
    "max_seq_len",
    "db_load_mode",
    "gpu",
    "gpu_server",
    "gpu_server_wait_timeout",
    "mpi_runner",
    "force_reuse",
    "remove_tmp_files",
    "filter_hits",
    "sort_results",
    "create_lookup",
    "chain_alignments",
    "merge_query",
    "strand",
    "sort_by_structure_bits",
    "tmscore_threshold",
    "tmscore_threshold_mode",
    "tmalign_hit_order",
    "tmalign_fast",
    "lddt_threshold",
    "alignment_type",
    "exact_tmscore",
    "cluster_search",
    *_COMMON_OPTIONS,
]
_CONVERT_ALIGNMENTS_OPTIONS = [
    "gap_open",
    "gap_extend",
    "format_output",
    "translation_table",
    "search_type",
    "exact_tmscore",
    "sub_mat",
    "db_load_mode",
    *_COMMON_OPTIONS,
]
_RESULT_TO_MSA_OPTIONS = [
    "comp_bias_corr",
    "comp_bias_corr_scale",
    "gap_open",
    "gap_extend",
    "filter_msa",
    "filter_min_enable",
    "max_seq_id",
    "qid",
    "qsc",
    "cov",
    "diff",
    "allow_deletion",
    "sub_mat",
    "db_load_mode",
    "summary_prefix",
    "skip_query",
    *_COMMON_OPTIONS,
]


class MMseqsLikeApp(LocalApp):
    """
    Shared command implementation for MMseqs2-compatible applications.

    Parameters
    ----------
    path : path-like
        Path to the application executable.
    """

    _tmp_prefix: ClassVar[str]

    def __init__(self, path: PathLike[str] | str) -> None:
        super().__init__(path)
        self._tmp_dir = Path(gettempdir())

    @property
    def tmp_dir(self) -> Path:
        """
        The temporary directory required for some *MMseqs2* and *Foldseek*
        commands.

        Returns
        -------
        tmp_dir : Path
            The parent directory for temporary command workspaces.
        """
        return self._tmp_dir

    @tmp_dir.setter
    def tmp_dir(self, tmp_dir: PathLike[str] | str) -> None:
        """
        Set the temporary directory required for some *MMseqs2* and *Foldseek*
        commands.

        Parameters
        ----------
        tmp_dir : path-like
            The parent directory for temporary command workspaces.
        """
        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(exist_ok=True, parents=True)
        self._tmp_dir = tmp_dir

    def _format_value(self, value: Any) -> str:
        match value:
            case Database():
                if not value.is_compatible_with(self):
                    raise ValueError(
                        f"Database is not compatible with {type(self).__name__}"
                    )
                return str(value.name)
            case tuple() | list():
                return ",".join(self._format_value(element) for element in value)
            case _:
                return super()._format_value(value)

    @command(subcommand="databases", allowed_options=_DATABASES_OPTIONS)
    def databases(
        self,
        name: str,
        **kwargs: Any,
    ) -> CommandSetup[SequenceDatabase[Any]]:
        """
        Download and prepare a named database using ``databases``.

        Parameters
        ----------
        name : str
            The database name understood by the application.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of SequenceDatabase
            A handle to the running download.

        Examples
        --------
        >>> database = MMseqsApp().databases("UniRef50").result()  # doctest: +SKIP
        """
        database = SequenceDatabase(self)
        work_dir = self._create_work_dir()
        return CommandSetup(
            parameters=[
                CLIArgument(name),
                CLIArgument(database),
                CLIArgument(work_dir),
            ],
            evaluate=lambda _out, _err: database,
        )

    @command(subcommand="createdb", allowed_options=_CREATE_DB_OPTIONS)
    def create_db(
        self,
        input_paths: (
            Iterable[PathLike[str] | str | IO[str]] | PathLike[str] | str | IO[str]
        ),
        db_type: DatabaseType | None = None,
        **kwargs: Any,
    ) -> CommandSetup[SequenceDatabase[Any]]:
        """
        Create a sequence or structure database using ``createdb``.

        Parameters
        ----------
        input_paths : path-like, text file or iterable thereof
            FASTA/FASTQ inputs for MMseqs2, or structural inputs for Foldseek.
            A text file must expose its file-system path through its ``name``
            attribute, and pending writes must be flushed before calling this
            method.
        db_type : DatabaseType, optional
            The type of the input sequences. By default, *MMseqs2* detects the
            type automatically. This option is not supported by *Foldseek*.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of SequenceDatabase
            A handle to the running conversion.

        Examples
        --------
        >>> database = MMseqsApp().create_db(path_to_sequences / "prot.fasta", threads=1).result()
        """
        if isinstance(input_paths, (str, PathLike, IOBase)):
            inputs = [input_paths]
        else:
            inputs = list(input_paths)
        if len(inputs) == 0:
            raise ValueError("At least one input is required")
        database = SequenceDatabase(self)
        parameters: list[CLIArgument | CLIOption] = [
            *(
                CLIArgument(
                    Path(input) if isinstance(input, (str, PathLike)) else input
                )
                for input in inputs
            ),
            CLIArgument(database),
        ]
        if db_type is not None:
            parameters.append(CLIOption("dbtype", db_type))
        return CommandSetup(
            parameters=parameters,
            evaluate=lambda _out, _err: database,
        )

    @command(subcommand="createsubdb", allowed_options=_CREATE_SUB_DB_OPTIONS)
    def create_sub_db(
        self,
        subset: Database[Any] | PathLike[str] | str,
        sequence_db: SequenceDatabase[Any],
        **kwargs: Any,
    ) -> CommandSetup[SequenceDatabase[Any]]:
        """
        Create a sequence database containing a selected subset using
        ``createsubdb``.

        Parameters
        ----------
        subset : Database or path-like
            A database or file defining the retained database keys.
        sequence_db : SequenceDatabase
            The source sequence database.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of SequenceDatabase
            A handle to the running subset creation.

        Examples
        --------
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> # Trivial case: Select all entries from the original database
        >>> subset = app.create_sub_db(database, database).result()
        """
        database = SequenceDatabase(self)
        return CommandSetup(
            parameters=[
                CLIArgument(subset),
                CLIArgument(sequence_db),
                CLIArgument(database),
            ],
            evaluate=lambda _out, _err: database,
        )

    @overload
    def convert_to_fasta(
        self,
        sequence_db: SequenceDatabase[Any],
        fasta_path: None = None,
        **kwargs: Any,
    ) -> LocalProcessFuture[IO[bytes]]: ...

    @overload
    def convert_to_fasta(
        self,
        sequence_db: SequenceDatabase[Any],
        fasta_path: PathLike[str] | str,
        **kwargs: Any,
    ) -> LocalProcessFuture[Path]: ...

    @command(subcommand="convert2fasta", allowed_options=_CONVERT_TO_FASTA_OPTIONS)
    def convert_to_fasta(
        self,
        sequence_db: SequenceDatabase[Any],
        fasta_path: PathLike[str] | str | None = None,
        **kwargs: Any,
    ) -> CommandSetup[Path | IO[bytes]]:
        """
        Export a sequence database in FASTA format using ``convert2fasta``.

        Parameters
        ----------
        sequence_db : SequenceDatabase
            The database to export.
        fasta_path : path-like, optional
            The output FASTA path. If omitted, the output is written to a
            temporary binary file.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of Path or binary file
            A handle resolving to the output path, or to the temporary file if
            `fasta_path` is omitted.

        Examples
        --------
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> fasta_file = app.convert_to_fasta(database).result()
        >>> fasta_file.read().startswith(b">")
        True
        >>> fasta_file.close()
        """
        output: Path | IO[bytes]
        if fasta_path is None:
            output = NamedTemporaryFile("w+b")
        else:
            output = Path(fasta_path)
        return CommandSetup(
            parameters=[CLIArgument(sequence_db), CLIArgument(output)],
            evaluate=lambda _out, _err: output,
        )

    @overload
    def convert_alignments(
        self,
        alignment_db: AlignmentDatabase[Any],
        output_path: None = None,
        format_mode: AlignmentFormatMode | None = None,
        **kwargs: Any,
    ) -> LocalProcessFuture[IO[bytes]]: ...

    @overload
    def convert_alignments(
        self,
        alignment_db: AlignmentDatabase[Any],
        output_path: PathLike[str] | str,
        format_mode: AlignmentFormatMode | None = None,
        **kwargs: Any,
    ) -> LocalProcessFuture[Path]: ...

    @command(subcommand="convertalis", allowed_options=_CONVERT_ALIGNMENTS_OPTIONS)
    def convert_alignments(
        self,
        alignment_db: AlignmentDatabase[Any],
        output_path: PathLike[str] | str | None = None,
        format_mode: AlignmentFormatMode | None = None,
        **kwargs: Any,
    ) -> CommandSetup[Path | IO[bytes]]:
        """
        Export pairwise alignments in a human-readable format using ``convertalis``.

        Parameters
        ----------
        alignment_db : AlignmentDatabase
            The alignment database to export.
        output_path : path-like, optional
            The output path. If omitted, the output is written to a temporary
            binary file.
        format_mode : AlignmentFormatMode, optional
            The output format. By default, the application writes a BLAST-like
            table.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of Path or binary file
            A handle resolving to the output path, or to the temporary file if
            `output_path` is omitted.

        Examples
        --------
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> alignments = app.search(database, database, threads=1).result()
        >>> alignment_file = app.convert_alignments(alignments).result()
        >>> len(alignment_file.readline()) > 0
        True
        >>> alignment_file.close()
        """
        output: Path | IO[bytes]
        if output_path is None:
            output = NamedTemporaryFile("w+b")
        else:
            output = Path(output_path)
        parameters: list[CLIArgument | CLIOption] = [
            CLIArgument(alignment_db.query_db),
            CLIArgument(alignment_db.target_db),
            CLIArgument(alignment_db),
            CLIArgument(output),
        ]
        if format_mode is not None:
            parameters.append(CLIOption("format_mode", format_mode))
        return CommandSetup(
            parameters=parameters,
            evaluate=lambda _out, _err: output,
        )

    @command(subcommand="result2msa", allowed_options=_RESULT_TO_MSA_OPTIONS)
    def result_to_msa(
        self,
        alignment_db: AlignmentDatabase[Any],
        msa_format_mode: MSAFormatMode | None = None,
        **kwargs: Any,
    ) -> CommandSetup[MSADatabase[Any]]:
        """
        Convert pairwise alignment results into multiple alignments using
        ``result2msa``.

        Parameters
        ----------
        alignment_db : AlignmentDatabase
            The pairwise alignment results.
        msa_format_mode : MSAFormatMode, optional
            The multiple alignment format.
            The flat-file :attr:`MSAFormatMode.STOCKHOLM` mode is not supported because
            this method returns an :class:`MSADatabase`.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of MSADatabase
            A handle to the running conversion.

        Examples
        --------
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> alignments = app.search(database, database, threads=1).result()
        >>> msa_database = app.result_to_msa(alignments, threads=1).result()
        >>> isinstance(msa_database, MSADatabase)
        True
        """
        if msa_format_mode == MSAFormatMode.STOCKHOLM:
            raise ValueError(
                "The STOCKHOLM format produces a flat file instead of an MSADatabase"
            )
        database = MSADatabase(self)
        parameters: list[CLIArgument | CLIOption] = [
            CLIArgument(alignment_db.query_db),
            CLIArgument(alignment_db.target_db),
            CLIArgument(alignment_db),
            CLIArgument(database),
        ]
        if msa_format_mode is not None:
            parameters.append(CLIOption("msa_format_mode", msa_format_mode))
        return CommandSetup(
            parameters=parameters,
            evaluate=lambda _out, _err: database,
        )

    @command(subcommand="lndb", allowed_options=["v"])
    def link_db(
        self,
        source: PathLike[str] | str,
        destination: PathLike[str] | str,
        **kwargs: Any,
    ) -> CommandSetup[Path]:
        """
        Create a symbolic database link using ``lndb``.

        Parameters
        ----------
        source : path-like
            The source database prefix.
        destination : path-like
            The destination database prefix.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of Path
            A handle resolving to the destination prefix.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> with tempfile.NamedTemporaryFile() as link_file:
        ...     link_path = Path(link_file.name)
        >>> _ = app.link_db(database.name, link_path).result()
        >>> link_path.is_symlink()
        True
        >>> link_path.unlink()
        """
        source = Path(source)
        destination = Path(destination)
        return CommandSetup(
            parameters=[CLIArgument(source), CLIArgument(destination)],
            evaluate=lambda _out, _err: destination,
        )

    def _create_work_dir(self) -> TemporaryDirectory[str]:
        """
        Create an isolated command workspace below the temporary directory.

        Returns
        -------
        TemporaryDirectory
            The created temporary directory.
        """
        tmp_dir = self.tmp_dir
        tmp_dir.mkdir(exist_ok=True, parents=True)
        return TemporaryDirectory(prefix=self._tmp_prefix, dir=tmp_dir)


class MMseqsApp(MMseqsLikeApp):
    """
    A reusable handle to the *MMseqs2* command line program.

    Parameters
    ----------
    path : path-like, optional
        Path to the ``mmseqs`` executable.

    Attributes
    ----------
    tmp_dir : Path
        The parent directory for temporary command workspaces.

    Examples
    --------
    >>> future = MMseqsApp().create_db(path_to_sequences / "prot.fasta", threads=1)
    >>> isinstance(future.result(), SequenceDatabase)
    True
    """

    _tmp_prefix = "mmseqs_"

    def __init__(self, path: PathLike[str] | str = "mmseqs") -> None:
        super().__init__(path)

    @command(subcommand="search", allowed_options=_SEARCH_OPTIONS)
    def search(
        self,
        query_db: SequenceDatabase[Any],
        target_db: SequenceDatabase[Any],
        matrix: SubstitutionMatrix[Any, Any] | None = None,
        gap_penalty: int | tuple[int, int] | None = None,
        **kwargs: Any,
    ) -> CommandSetup[AlignmentDatabase[Any]]:
        """
        Search a target database for matches to every query using ``search``.

        Parameters
        ----------
        query_db, target_db : SequenceDatabase
            The query and target databases.
        matrix : SubstitutionMatrix, optional
            A custom symmetric substitution matrix. By default, the standard
            *MMseqs2* matrix is used.
        gap_penalty : int or tuple of (int, int), optional
            The negative gap penalty. A single value applies to both opening
            and extension; a tuple gives the respective affine penalties.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of AlignmentDatabase
            A handle to the running search.

        Examples
        --------
        >>> app = MMseqsApp()
        >>> database = app.create_db(path_to_sequences / "prot.fasta").result()
        >>> matrix = SubstitutionMatrix.std_protein_matrix()
        >>> future = app.search(
        ...     database, database, matrix=matrix,
        ...     gap_penalty=(-11, -1), threads=1,
        ... )
        >>> isinstance(future.result(), AlignmentDatabase)
        True
        """
        options = []
        if matrix is not None:
            matrix_file = NamedTemporaryFile("w", suffix=".out")
            _write_substitution_matrix(matrix, matrix_file)
            matrix_file.flush()
            options.append(CLIOption("sub_mat", matrix_file))
        gap_open, gap_extend = resolve_gap_penalty(gap_penalty)
        if gap_open is not None and gap_extend is not None:
            options += [
                CLIOption("gap_open", -gap_open),
                CLIOption("gap_extend", -gap_extend),
            ]
        database = AlignmentDatabase(self, query_db, target_db)
        work_dir = self._create_work_dir()
        return CommandSetup(
            parameters=[
                *options,
                CLIArgument(query_db),
                CLIArgument(target_db),
                CLIArgument(database),
                CLIArgument(work_dir),
            ],
            evaluate=lambda _out, _err: database,
        )


class FoldseekApp(MMseqsLikeApp):
    """
    A reusable handle to the *Foldseek* command line program.

    Parameters
    ----------
    path : path-like, optional
        Path to the ``foldseek`` executable.

    Attributes
    ----------
    tmp_dir : Path
        The parent directory for temporary command workspaces.

    Examples
    --------
    >>> structure_path = path_to_structures / "1aki.cif"
    >>> future = FoldseekApp().create_db(structure_path, threads=1)
    >>> isinstance(future.result(), SequenceDatabase)
    True
    """

    _tmp_prefix = "foldseek_"

    def __init__(self, path: PathLike[str] | str = "foldseek") -> None:
        super().__init__(path)

    @command(subcommand="search", allowed_options=_SEARCH_OPTIONS)
    def search(
        self,
        query_db: SequenceDatabase[Any],
        target_db: SequenceDatabase[Any],
        **kwargs: Any,
    ) -> CommandSetup[AlignmentDatabase[Any]]:
        """
        Search a target database for structural matches using ``search``.

        Parameters
        ----------
        query_db, target_db : SequenceDatabase
            The query and target structure databases.
        **kwargs
            Additional command line options.

        Returns
        -------
        future : Future of AlignmentDatabase
            A handle to the running search.

        Examples
        --------
        >>> app = FoldseekApp()
        >>> database = app.create_db(path_to_structures / "1aki.cif").result()
        >>> future = app.search(database, database, threads=1)
        >>> isinstance(future.result(), AlignmentDatabase)
        True
        """
        database = AlignmentDatabase(self, query_db, target_db)
        work_dir = self._create_work_dir()
        return CommandSetup(
            parameters=[
                CLIArgument(query_db),
                CLIArgument(target_db),
                CLIArgument(database),
                CLIArgument(work_dir),
            ],
            evaluate=lambda _out, _err: database,
        )


def _write_substitution_matrix(
    matrix: SubstitutionMatrix[Any, Any], file: IO[str]
) -> None:
    """
    Write a substitution matrix in the MMseqs2 matrix-file format.
    """
    if not matrix.is_symmetric():
        raise ValueError("MMseqs2 requires a symmetric substitution matrix")
    alphabet = matrix.get_alphabet1()
    include_mask = np.ones(len(alphabet), dtype=bool)
    if ProteinSequence.alphabet.extends(alphabet):
        # Ambiguous and stop amino-acid symbols are derived from canonical
        # residues and do not represent independent log-odds states.
        include_mask[ProteinSequence.alphabet.encode("B") :] = False
    indices = np.flatnonzero(include_mask)
    symbols = [
        str(matrix.get_alphabet1().decode(int(index))).upper() for index in indices
    ]
    if any(len(symbol) != 1 for symbol in symbols):
        raise ValueError("MMseqs2 matrix symbols must be single characters")
    if "X" in symbols or len(symbols) > 20:
        raise ValueError("MMseqs2 matrices support up to 20 defined symbols")
    active_scores = matrix.score_matrix()[np.ix_(indices, indices)]

    # MMseqs2 requires an explicit final row and column for its unknown state.
    scores = np.zeros((len(indices) + 1, len(indices) + 1), dtype=np.int32)
    scores[:-1, :-1] = active_scores
    original_symbols = [str(symbol).upper() for symbol in matrix.get_alphabet1()]
    if "X" in original_symbols:
        unknown_index = original_symbols.index("X")
        scores[-1, :-1] = matrix.score_matrix()[unknown_index, indices]
        scores[:-1, -1] = matrix.score_matrix()[indices, unknown_index]
        scores[-1, -1] = matrix.score_matrix()[unknown_index, unknown_index]
    symbols.append("X")

    file.write("# Generated by Biotite\n")
    file.write("    " + " ".join(f"{symbol:>3}" for symbol in symbols) + "\n")
    for symbol, row in zip(symbols, scores, strict=True):
        file.write(f"{symbol:>3} " + " ".join(f"{score:>3d}" for score in row) + "\n")
