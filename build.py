from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


extensions = cythonize(
    "src/biotite/**/*.pyx",
    include_path=[numpy.get_include()],
    language_level=3
)
# Do not use deprecated NumPy API
# Requires unpublished Cython 3.0
#for ext in extensions:
#    ext.define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]


def build(setup_kwargs):
    setup_kwargs.update({
        'ext_modules': extensions
    })