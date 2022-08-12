We are pleased to announce the release 3.6.0.rc1 of the GUDHI library.

As a major new feature, the GUDHI library now offers automatic differentiation for the computation of
persistence diagrams, Cubical complex persistence scikit-learn like interface, datasets fetch methods,
and weighted version for alpha complex in any dimension D.

Do not hesitate to [fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel). From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).
For further information, please visit the [GUDHI web site](https://gudhi.inria.fr/).

# GUDHI 3.6.0rc1 Release Notes
Below is a list of changes made since GUDHI 3.5.0:

- TensorFlow 2 models that can handle automatic differentiation for the computation of persistence diagrams:
     - [Cubical complex](https://gudhi.inria.fr/python/latest/cubical_complex_tflow_itf_ref.html)
     - [lower-star persistence on simplex trees](https://gudhi.inria.fr/python/latest/ls_simplex_tree_tflow_itf_ref.html)
     - [Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_tflow_itf_ref.html)

- [Cubical complex](https://gudhi.inria.fr/python/latest/cubical_complex_sklearn_itf_ref.html)
     - Cubical complex persistence scikit-learn like interface

- [Datasets](https://gudhi.inria.fr/python/latest/datasets.html)
     - `datasets.remote.fetch_bunny` and `datasets.remote.fetch_spiral_2d` allows to fetch datasets from [GUDHI-data](https://github.com/GUDHI/gudhi-data)

- [Alpha complex](https://gudhi.inria.fr/python/latest/alpha_complex_user.html)
     - python weighted version for alpha complex is now available in any dimension D.
     - `alpha_complex = gudhi.AlphaComplex(off_file='/data/points/tore3D_300.off')` is deprecated, please use [read_points_from_off_file](https://gudhi.inria.fr/python/latest/point_cloud.html#gudhi.read_points_from_off_file) instead.

- [Edge collapse](https://gudhi.inria.fr/doc/latest/group__edge__collapse.html)
     - rewriting of the module to improve performance

- [Čech complex](https://gudhi.inria.fr/doc/latest/group__edge__collapse.html)
     - rewriting of the module to improve performance

- [Representations](https://gudhi.inria.fr/python/latest/representations.html#gudhi.representations.vector_methods.BettiCurve)
     - A more flexible Betti curve class capable of computing exact curves

- [C++ documentation](https://gudhi.inria.fr/doc/latest/)
     - upgrade and improve performance with new doxygen features

- [Simplex tree](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - `__deepcopy__`, `copy` and copy constructors for python module
     - `expansion_with_blockers` python interface

- Installation
     - Boost &ge; 1.66.0 is now required (was &ge; 1.56.0).
     - Python >= 3.5 and cython >= 0.27 are now required.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.5.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.6.0+is%3Aclosed) is available on GitHub.

## Contributors

- @albert-github
- @gspr
- @Hind-M
- @MathieuCarriere
- @mglisse
- @Soriano-Trigueros
- @VincentRouvreau
