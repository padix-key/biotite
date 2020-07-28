# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.application.mmseqs"
__author__ = "Patrick Kunzmann"
__all__ = ["MMseqsSearchApp"]

import string
from tempfile import NamedTemporaryFile, gettempdir
from subprocess import SubprocessError
import numpy as np
from ..localapp import LocalApp
from ..application import AppState, requires_state
from ...sequence.alphabet import LetterAlphabet
from ...sequence.align.alignment import Alignment
from ...sequence.align.matrix import SubstitutionMatrix
from ...sequence.io.fasta.file import FastaFile
from ...sequence.io.fasta.convert import set_sequence
from ...sequence.sequence import Sequence
from ...sequence.seqtypes import NucleotideSequence, ProteinSequence


class _MappedSequence(Sequence):

    MAX_ALPH_LENGTH = len(string.ascii_uppercase)

    def __init__(self, alph_length):
        if alph_length > _MappedSequence.MAX_ALPH_LENGTH:
            raise ValueError("Alphabet length is too large")
        self._alphabet = LetterAlphabet(string.ascii_uppercase[:alph_length])
    
    def get_alphabet(self):
        return self._alphabet


class MMseqsSearchApp(LocalApp):
    """
    """
    
    def __init__(self, query, target, matrix, bin_path="mmseqs"):
        super().__init__(bin_path)

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
        
        self._query = MMseqsSearchApp._map_sequence(query, self._alph)
        self._target = MMseqsSearchApp._map_sequence(target, self._alph)
        self._matrix = MMseqsSearchApp._map_matrix(matrix)

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
        
        matrix_str = "\n".join(
            ["\t".join(line.split())
             for line in str(self._matrix).splitlines()]
        )
        self._matrix_file.write(matrix_str)
        self._matrix_file.flush()
        ###
        with open("test_matrix.out", "w") as f:
            f.write(matrix_str)
        ###
        
        self.set_arguments(
            [
                "easy-search",
                self._query_file.name,
                self._target_file.name,
                self._out_file.name,
                gettempdir(),
                "-s", "10.0",
                "-k", "5",
                "--alph-size", str(len(self._alph)),
                "--seed-sub-mat", self._matrix_file.name,
                "--sub-mat", self._matrix_file.name,
                #"--seed-sub-mat", MATRIX_FILE,
                #"--sub-mat", MATRIX_FILE,
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

    
    @staticmethod
    def _map_sequence(sequence, alphabet):
        if len(alphabet) > _MappedSequence.MAX_ALPH_LENGTH:
            # Cannot map into target alphabet if the alphabet
            # has more symbols
            raise TypeError(
                f"The software cannot align sequences of type "
                f"{type(sequence).__name__}: "
                f"Alphabet is too large to be converted into single letters"
            )
        # Mapping is done by simply taking over the sequence
        # code of the original sequence
        mapped_seq = _MappedSequence(len(alphabet))
        mapped_seq.code = sequence.code
        return mapped_seq
    
    @staticmethod
    def _map_matrix(matrix):
        if not matrix.is_symmetric():
            raise ValueError("A symmetric matrix is required")
        # Create a new substitution matrix with the values taken
        # from the original matrix
        # Only exchange the alphabets
        old_alph = matrix.get_alphabet1()
        new_alph = _MappedSequence(len(old_alph)).get_alphabet()
        return SubstitutionMatrix(new_alph, new_alph, matrix.score_matrix())
           