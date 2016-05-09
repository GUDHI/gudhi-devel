from distutils.core import setup, Extension
from Cython.Build import cythonize

gudhi = Extension(
    "gudhi",
    sources = ['gudhi.pyx',],
    language = 'c++',
    extra_compile_args=['-std=c++11'],
    include_dirs = ['../include','./src/cpp'],
)

setup(
    name = 'gudhi',
    author='Vincent Rouvreau',
    author_email='gudhi-contact@lists.gforge.inria.fr',
    version='0.1.0',
    url='http://gudhi.gforge.inria.fr/',
    ext_modules = cythonize(gudhi),
)
