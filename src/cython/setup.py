from distutils.core import setup, Extension
from Cython.Build import cythonize

gudhi = Extension(
    "gudhi",
    sources = ['gudhi.pyx',],
    language = 'c++',
    extra_compile_args=['-frounding-math','-std=c++11','-DCGAL_EIGEN3_ENABLED','-DCGAL_USE_GMP','-DCGAL_USE_GMPXX','-DCGAL_USE_MPFR'],
    libraries=['mpfr','gmpxx','gmp','CGAL'],
    library_dirs=['/usr/local/lib/'],
    include_dirs = ['../include','./src/cpp','/usr/local/include/eigen3'],
)

setup(
    name = 'gudhi',
    author='Vincent Rouvreau',
    author_email='gudhi-contact@lists.gforge.inria.fr',
    version='0.1.0',
    url='http://gudhi.gforge.inria.fr/',
    ext_modules = cythonize(gudhi),
)
