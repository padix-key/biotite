# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.application.mmseqs"
__author__ = "Patrick Kunzmann"
__all__ = ["MMseqsSearchApp"]

from tempfile import NamedTemporaryFile, gettempdir
from subprocess import SubprocessError
from ..localapp import LocalApp
from ..application import AppState, requires_state
from ...sequence.io.fasta.file import FastaFile
from ...sequence.io.fasta.convert import set_sequence
from ...sequence.seqtypes import NucleotideSequence, ProteinSequence


class MMseqsSearchApp(LocalApp):
    """
    """
    
    def __init__(self, query, target, bin_path="mmseqs"):
        super().__init__(bin_path)
        if isinstance(query, (tuple, list)):
            self._queries = list(query)
        else:
            self._queries = [query]
        self._target = target

        seq_type = type(self._target)
        if isinstance(self._target, NucleotideSequence):
            self._search_type = 3
        elif isinstance(self._target, ProteinSequence):
            self._search_type = 1
        else:
            raise TypeError(
                f"MMseqs2 can only align nucleotide and protein sequences, "
                f"not '{type(self._target).__name__}'"
            )

        for query in self._queries:
            if not isinstance(query, type(self._target)):
                raise TypeError(
                    "Query sequences must be the same sequence type as the "
                    "target sequence"
                )

        self._query_files = [
            NamedTemporaryFile("w", suffix=".fasta") for _ in self._queries
        ]
        self._target_file = NamedTemporaryFile("w", suffix=".fasta")
        self._out_file = NamedTemporaryFile("r", suffix=".out")

    def run(self):
        for sequence, file in zip(
            self._queries + [self._target],
            self._query_files + [self._target_file]
        ):
            in_file = FastaFile()
            set_sequence(in_file, sequence)
            in_file.write(file)
            file.flush()
        
        self.set_arguments(
            ["easy-search"] +
            [file.name for file in self._query_files] +
            [self._target_file.name, self._out_file.name, gettempdir()] +
            [
                "--search-type", str(self._search_type),
                "--format-output", ",".join([
                    "query",
                    "target",
                    "raw",
                    "pident",
                    "evalue",
                    "qstart",
                    "qend",
                    "tstart",
                    "tend",
                    "qaln",
                    "taln",
                ]),
                "-s", "50.0",
                "-k", "4"
            ]
        )
        
        super().run()
    
    def evaluate(self):
        exit_code = self.get_exit_code()
        if exit_code != 0:
            # MMseqs uses STDOUT for error messages
            err_msg = self.get_stdout()
            raise SubprocessError(
                f"MMseqs2 returned with exit code {exit_code}: "
                f"{err_msg}"
            )

        super().evaluate()

        # Somehow reading the file with
        #
        #     output = self._out_file.read()
        #
        # does not work, the file appears always empty
        with open(self._out_file.name, "r") as file:
            output = file.read()
        print(output)
    
    def clean_up(self):
        super().clean_up()
        for file in self._query_files + [self._target_file, self._out_file]:
            file.close()
    
    @requires_state(AppState.JOINED)
    def get_alignments(self):
        pass
