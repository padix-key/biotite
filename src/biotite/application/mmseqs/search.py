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
from ...sequence.align.alignment import Alignment
from ...sequence.io.fasta.file import FastaFile
from ...sequence.io.fasta.convert import set_sequence
from ...sequence.seqtypes import NucleotideSequence, ProteinSequence


class MMseqsSearchApp(LocalApp):
    """
    """
    
    def __init__(self, query, target, matrix, bin_path="mmseqs"):
        super().__init__(bin_path)
        #if isinstance(query, (tuple, list)):
        #    self._queries = list(query)
        #else:
        #    self._queries = [query]
        self._target = target
        self._query = query
        self._matrix = matrix

        if not matrix.is_symmetric():
            raise ValueError("A symmetric substitution matrix is required")
        self._alph = matrix.get_alphabet1()
        if not self._alph.extends(query.alphabet):
            raise TypeError(
                "The query sequence's alphabet is incompatible with the "
                "substitution matrix alphabet"
            )
        if not self._alph.extends(target.alphabet):
            raise TypeError(
                "The target sequence's alphabet is incompatible with the "
                "substitution matrix alphabet"
            )

        self._query_file = NamedTemporaryFile("w", suffix=".fasta")
        self._target_file = NamedTemporaryFile("w", suffix=".fasta")
        self._out_file = NamedTemporaryFile("r", suffix=".results")
        self._matrix_file = NamedTemporaryFile("w", suffix=".out")

    def run(self):
        for sequence, file in zip(
            (self._query, self._target),
            (self._query_file, self._target_file)
        ):
            in_file = FastaFile()
            set_sequence(in_file, sequence)
            in_file.write(file)
            file.flush()
        #self._matrix_file.write(str(self._matrix))
        #self._matrix_file.flush()
        """
        with open("matrix_test.out", "w") as f:
            matrix_str = "\n".join(
                ["\t".join(line.split())
                 for line in str(self._matrix).splitlines()]
            )
            print(matrix_str)
            f.write(matrix_str)
        """
        
        MATRIX_FILE = "matrix_test.out"
        self.set_arguments(
            [
                "easy-search",
                self._query_file.name,
                self._target_file.name,
                self._out_file.name,
                gettempdir(),
                "-s", "10.0",
                "-k", "5",
                #"--alph-size", str(len(self._alph)),
                #"--seed-sub-mat", self._matrix_file.name,
                #"--sub-mat", self._matrix_file.name,
                "--seed-sub-mat", MATRIX_FILE,
                "--sub-mat", MATRIX_FILE,
                "--dbtype", "1",
                "--search-type", "1",
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
                #"-s", "5.0",
                #"-k", "4",
                #"--gap-open", str(2),
                #"--gap-extend", str(2),
            ]
        )
        
        super().run()
    
    def evaluate(self):
        super().evaluate()

        # Somehow reading the file with
        #
        #     output = self._out_file.read()
        #
        # does not work, the file appears always empty
        with open(self._out_file.name, "r") as file:
            output = file.read()
        print(output)
        print()
        for i, line in enumerate(output.splitlines()):
            splitted = line.split()
            if len(splitted) == 0:
                continue
            elif len(splitted) != 11:
                raise ValueError(
                    f"Line {i+1} has {len(splitted)} values, expected 11"
                )
            query_id        = splitted[0]
            target_id       = splitted[1]
            score           = int(splitted[2])
            identity        = float(splitted[3])
            e_value         = float(splitted[4])
            query_interval  = (int(splitted[5])-1, int(splitted[6]))
            target_interval = (int(splitted[7])-1, int(splitted[8]))
            trace = Alignment.trace_from_strings(
                [splitted[9], splitted[10]]
            )
            trace[:,0] += query_interval[0]
            trace[:,1] += target_interval[0]
            alignment = Alignment(
                [self._query, self._target],
                trace,
                score
            )
            print(alignment)
    
    def clean_up(self):
        super().clean_up()
        for file in [
            self._query_file, self._target_file,
            self._out_file, self._matrix_file
        ]:
            file.close()
    
    @requires_state(AppState.JOINED)
    def get_alignments(self):
        pass
