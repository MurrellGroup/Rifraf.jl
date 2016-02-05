from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'Quiver2',
    ext_modules = cythonize("_quiver2.pyx"),
)

