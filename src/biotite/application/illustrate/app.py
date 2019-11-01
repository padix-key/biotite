# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"
__all__ = ["IllustrateApp"]

from biotite.structure.io import save_structure
from biotite.application import LocalApp
from biotite import temp_file


class IllustrateApp(LocalApp):
    
    def __init__(self, image_file, structure, colors, radii,
                 bin_path="illustrate"):
        super().__init__(bin_path)
        self.ignore_stdout()

        self._structure = structure
        self._colors = colors
        self._radii = radii

        self._x_rot = 0
        self._y_rot = 0
        self._z_rot = 0

        self._structure_file_name = temp_file("pdb")
        self._command_file_name = temp_file("inp")
        self._image_file_name = image_file
    
    def run(self):
        save_structure(self._structure_file_name, self._structure)
        
        # Write the Illustrate command file
        lines = []

        # Structure file input
        lines.extend(["read", self._structure_file_name])

        # Color and radius selectors
        # In this interface the color and radius is described for each
        # atom -> no wildcards are necessary
        for i in range(len(self._structure)):
            het = "HETATM" if self._structure.hetero[i] else "ATOM"
            descriptor = (
                f"{self._structure.atom_name[i]:4} "
                f"{self._structure.res_name[i]:3} "
                f"{self._structure.chain_id[i]:1}"
            )
            residue = int(self._structure.res_id[i])
            r, g, b = self._colors[i]
            radius = self._radii[i]
            lines.append(
                f"{het:<6}{descriptor} {residue},{residue}, "
                f"{r:3.1f},{g:3.1f},{b:3.1f}, {radius:3.1f}"
            )
        lines.append("END")

        lines.extend(["xrot", f"{self._x_rot:<.1f}"])
        lines.extend(["yrot", f"{self._y_rot:<.1f}"])
        lines.extend(["zrot", f"{self._z_rot:<.1f}"])

        # Preliminary: other options
        lines.extend([
            "center",
            "auto",
            "trans",
            "0.,0.,0.",
            "scale",
            "10.0",
            "wor",
            "1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0",
            "1,0.0023,2.0,1.0,0.2",
            "-30,-30",
            "illustrate",
            "3.0,10.0,4,0.0,5.0",
            "3.0,10.0",
            "3.0,8.0,6000.0"
        ])

        # Output image file
        lines.extend(["calculate", self._image_file_name])
        
        # Pipe command file into process' STDIN
        self.set_stdin("\n".join(lines) + "\n")

        super().run()
    
    def set_rotation(self, x=0, y=0, z=0):
        self._x_rot = x
        self._y_rot = y
        self._z_rot = z