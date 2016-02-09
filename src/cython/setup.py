from distutils.core import setup, Extension
from Cython.Build import cythonize

gudhi = Extension(
    "gudhi",
    sources = ['gudhi.pyx',],
    language = 'c++',
    extra_compile_args=['-std=c++11'],
    include_dirs = ['..'],
)

setup(
    name = 'gudhi',
    ext_modules = cythonize(gudhi),
)
