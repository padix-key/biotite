# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.application.openbabel"
__author__ = "Patrick Kunzmann"
__all__ = ["BabelDrawApp"]

from xml.etree import ElementTree
import numpy as np
from tempfile import NamedTemporaryFile
from ..localapp import LocalApp, cleanup_tempfile
from ..application import AppState, requires_state
from ...structure.error import BadStructureError
from ...structure.io.pdb import PDBFile

class BabelDrawApp(LocalApp):
    """
    Get coordinates for a 2D representation of any small molecule.
    """

    def __init__(self, molecule, bin_path="obabel"):
        super().__init__(bin_path)

        if molecule.bonds is None:
            raise BadStructureError("Molecule has no associated BondList")
        self._molecule = molecule

        self._in_file = NamedTemporaryFile(
            "w", suffix=".pdb", delete=False
        )
        self._out_file = NamedTemporaryFile(
            "r", suffix=".cdxml", delete=False
        )

    def run(self):
        pdb_file = PDBFile()
        pdb_file.set_structure(self._molecule)
        pdb_file.lines += _write_bonds(self._molecule.bonds)
        pdb_file.write(self._in_file)
        self._in_file.flush()
        self.set_arguments([
            "-i", "pdb", self._in_file.name,
            "-o", "cdxml", "-O", self._out_file.name
        ])
        super().run()

    def evaluate(self):
        super().evaluate()
        self._coordinates = _parse_coord(self._out_file)

    def clean_up(self):
        super().clean_up()
        cleanup_tempfile(self._in_file)
        cleanup_tempfile(self._out_file)

    @requires_state(AppState.JOINED)
    def get_coordinates(self):
        """
        Get coordinates for a 2D representation of any small molecule.
        """
        return self._coordinates


def _write_bonds(bonds):
    lines = []
    all_bonds, _ = bonds.get_all_bonds()
    for i in range(len(all_bonds)):
        if len(all_bonds[i]) > 0:
            line = f"CONECT {i+1:>4d}"
            for j in all_bonds[i]:
                if j == -1:
                    break
                line += f" {j+1:>4d}"
            lines.append(line)
    return lines


def _parse_coord(file):
    coord = []
    #print(file.read())
    root = ElementTree.parse(file).getroot()
    for atom_node in root.findall("./page/fragment/n"):
        coord.append([float(c) for c in atom_node.attrib["p"].split(" ")])
    return np.array(coord)