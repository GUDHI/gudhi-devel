Release History
===============

[Release 3.11.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.11.0)
-----------

Release date: February 2025

Below is a list of changes made since GUDHI 3.10.1:

- [Delaunay complex](https://gudhi.inria.fr/python/latest/delaunay_complex_user.html)
     - The Delaunay complex can be equipped with different filtrations:
          * Delaunay complex (no filtration values computed)
          * Delaunay-Čech complex (using minimal enclosing ball)
          * Alpha complex (moved in this new section)
     - The Delaunay-Čech and Alpha complex can output square, or not square, filtration values
     - An incremental version of the Delaunay complex (only in C++)

- [Rips complex persistence scikit-learn like interface](https://gudhi.inria.fr/python/latest/rips_complex_sklearn_itf_ref.html)
     - A binding to [Ripser](https://github.com/Ripser/ripser) when it accelerates the computation

- [Persistence graphical tools](https://gudhi.inria.fr/python/latest/persistence_graphical_tools_user.html)
     - Can now handle scikit-learn like interfaces outputs as inputs

- [Simplex tree](https://gudhi.inria.fr/doc/latest/class_gudhi_1_1_simplex__tree.html)
     - Can now store additionnal data on each simplex (only in C++)
     - Can be const

- Installation
     - CMake &ge; 3.15 is now required (was &ge; 3.8).
     - Python &ge; 3.8 is now required (was &ge; 3.5), because of `importlib.metadata`.
     - Support for Python 3.13 is now available

- Miscellaneous
     - The [list of bugs that were solved](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.11.0+is%3Aclosed)
         is available on GitHub.

[Release 3.10.1](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.10.1)
-----------

Release date: June 2024

Below is a list of changes made since GUDHI 3.10.0:

Only bug fixes have been implemented for this minor version.

The [list of bugs that were solved](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.10.1+is%3Aclosed) is available on GitHub.

[Release 3.10.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.10.0)
-----------

Release date: June 2024

Below is a list of changes made since GUDHI 3.9.0:

- [Persistence matrix](https://gudhi.inria.fr/doc/latest/group__persistence__matrix.html)
     > Matrix API is in a beta version and may change in incompatible ways in the near future.
     - Matrix structure for filtered complexes with multiple functionnalities related to persistence homology, such as
     representative cycles computation or vineyards.

- [Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_sklearn_itf_ref.html)
     - Rips complex persistence scikit-learn like interface

- [Čech complex](https://gudhi.inria.fr/cechcomplex/)
     - A new utility to compute the Delaunay-Čech filtration on a Delaunay triangulation.

- Installation
     - CGAL &ge; 5.1.0 is now required (was &ge; 4.11.0).
     - Eigen3 &ge; 3.3.0 is now required (was &ge; 3.1.0).

- Maintenance
     - Some bug fix for CGAL &ge; 6.0, NumPy &ge; 2.0, Scikit-learn &ge; 1.4, Matplotlib &ge; 3.6 and TensorFlow &ge; 2.16.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.9.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.10.0+is%3Aclosed) is available on GitHub.


[Release 3.9.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.9.0)
-----------

Release date: December 2023

Below is a list of changes made since GUDHI 3.8.0:

- [CubicalPersistence](https://gudhi.inria.fr/python/latest/cubical_complex_sklearn_itf_ref.html)
     - Much faster implementation for the 2d case with input from top-dimensional cells.

- [Simplex_tree](https://gudhi.inria.fr/doc/latest/group__simplex__tree.html)
     - A helper `for_each_simplex` that applies a given function object on each simplex
     - A new method `num_simplices_by_dimension` is now available thanks to this helper.
     - A `clear` method to empty the data stucture.
     - A new argument `ignore_infinite_values` for `initialize_filtration` method to skip infinite values. As a side effect, this change enhances the persistence computation.
     - `Simplex_tree_options_full_featured` has been renamed `Simplex_tree_options_default` and `Simplex_tree_options_python`.
     These are respectively the default options used by the `Simplex_tree` and by the python interface of the `SimplexTree` (as before this version).
     - From GUDHI 3.9.0, `Simplex_tree_options_full_featured` now activates `link_nodes_by_label` and `stable_simplex_handles` (making it slower, except for browsing cofaces).

     | Simplex_tree_options_*  | :warning: full_featured | default | python | minimal |
     | ---- | ---- | ---- | ---- | ---- |
     | store_key              | 1       | 1      | 1      | 0 |
     | store_filtration       | 1       | 1      | 1      | 0 |
     | contiguous_vertices    | 0       | 0      | 0      | 0 |
     | link_nodes_by_label    | ***1*** | 0      | 0      | 0 |
     | stable_simplex_handles | ***1*** | 0      | 0      | 0 |
     | Filtration_value       | double  | double | double |   |

- [Simplex_tree options](https://gudhi.inria.fr/doc/latest/struct_simplex_tree_options.html)
     - A new option `link_nodes_by_label` to speed up cofaces and stars access, when set to true.
     - A new option `stable_simplex_handles` to keep Simplex handles valid even after insertions or removals, when set to true.

- [Čech complex](https://gudhi.inria.fr/doc/latest/group__cech__complex.html)
     - A function `assign_MEB_filtration` that assigns to each simplex a filtration value equal to the squared radius of its minimal enclosing ball (MEB), given a simplicial complex and an embedding of its vertices. Applied on a Delaunay triangulation, it computes the Delaunay-Čech filtration.

- [Edge collapse](https://gudhi.inria.fr/python/latest/edge_collapse.html)
     - A Python function `reduce_graph` to simplify a clique filtration (represented as a sparse weighted graph), while preserving its persistent homology.

- [Mapper/GIC/Nerve complexes](https://gudhi.inria.fr/python/latest/cover_complex_sklearn_isk_ref.html)
     - A new method `save_to_html` to ease the Keppler Mapper visualization

- Installation
     - Boost &ge; 1.71.0 is now required (was &ge; 1.66.0).
     - cython >= 3.0.0 is now supported.
     - Python 3.12 pip package.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.8.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.9.0+is%3Aclosed) is available on GitHub.

[Release 3.8.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.8.0)
-----------

Release date: April 2023

As a major new feature, the GUDHI library now offers Perslay, a Tensorflow model for the representations module, scikit-learn like interfaces for Cover Complexes, a new function to compute persistence of a function on ℝ and the possibility to build a Cubical Complex as a lower-star filtration from vertices.

Below is a list of changes made since GUDHI 3.7.1:

- [Perslay](https://gudhi.inria.fr/python/latest/representations_tflow_itf_ref.html)
     - a TensorFlow layer for persistence diagrams representations.

- [Cover Complex](https://gudhi.inria.fr/python/latest/cover_complex_sklearn_user.html)
     - New classes to compute Mapper, Graph Induced complex and Nerves with a scikit-learn like interface.

- [Persistent cohomology](https://gudhi.inria.fr/doc/latest/group__persistent__cohomology.html)
     - New linear-time `compute_persistence_of_function_on_line`, also available though `CubicalPersistence` in Python.

- [Cubical complex](https://gudhi.inria.fr/doc/latest/group__cubical__complex.html)
     - Add possibility to build a lower-star filtration from vertices instead of top-dimensional cubes.
     - Naming the arguments is now mandatory in CubicalComplex python constructor.
     - Remove `newshape` mechanism from [CubicalPersistence](https://gudhi.inria.fr/python/latest/cubical_complex_sklearn_itf_ref.html)

- [Hera version of Wasserstein distance](https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html#hera)
     - now provides matching in its interface.

- [Subsampling](https://gudhi.inria.fr/doc/latest/group__subsampling.html)
     - New `choose_n_farthest_points_metric` as a faster alternative of `choose_n_farthest_points`.

- [SimplexTree](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - `SimplexTree` can now be used with `pickle`.
     - new `prune_above_dimension` method.

- Installation
     - CMake 3.8 is the new minimal standard to compile the library.
     - Support for oneAPI TBB (instead of deprecated TBB) to take advantage of multicore performance.
     - [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/) is the new sphinx theme of the python documentation.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.7.1](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.8.0+is%3Aclosed) is available on GitHub.

[Release 3.7.1](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.7.1)
-----------

Release date: January 2023

This minor post-release is a bug fix version for python representation module.

The [list of bugs that were solved since GUDHI-3.7.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.7.1+is%3Aclosed) is available on GitHub.


[Release 3.7.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.7.0)
-----------

Release date: December 2022

As a major new feature, the GUDHI library now offers new functions to initialize a Simplex tree.
Universal wheel for OSx pip package and python 3.11 are now available.

Below is a list of changes made since GUDHI 3.6.0:

- [Simplex tree](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - New functions to initialize from a matrix or insert batches of simplices of the same dimension.

- [Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_user.html)
     - Construction now rejects positional arguments, you need to specify `points=X`.

- Installation
     - C++17 is the new minimal standard to compile the library. This implies Visual Studio minimal version is now 2017.
     - OSx ARM pip package is now available thanks to a universal wheel
     - Python 3.11 pip package

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.6.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.7.0+is%3Aclosed) is available on GitHub.


[Release 3.6.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.6.0)
-----------

Release date: August 2022

As a major new feature, the GUDHI library now offers automatic differentiation for the computation of
persistence diagrams, Cubical complex persistence scikit-learn like interface, datasets fetch methods,
and weighted version for alpha complex in any dimension D.

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


[Release 3.5.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.5.0)
-----------

Release date: January 2022

As a major new feature, the GUDHI library now offers Coxeter triangulations and points generators.

Below is a list of changes made since GUDHI 3.4.1:

- [Coxeter triangulation](https://gudhi.inria.fr/doc/latest/group__coxeter__triangulation.html)
     - constructs a piecewise-linear approximation of an m-dimensional smooth manifold embedded in R^d using an ambient triangulation.

- [Datasets generators](https://gudhi.inria.fr/python/latest/datasets_generators.html)
     - the python module `points` enables the generation of points on a sphere or a flat torus.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.4.1](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.5.0+is%3Aclosed) is available on GitHub.


[Release 3.4.1](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.4.1)
-----------

Release date: January 2021

This minor release is a bug fix version to make GUDHI compile with CGAL 5.2.

The [list of bugs that were solved since GUDHI-3.4.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.4.1+is%3Aclosed) is available on GitHub.


[Release 3.4.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.4.0)
-----------

Release date: December 2020

As a major new feature, the GUDHI library now offers dD weighted alpha complex, pip and conda packages for Python 3.9.

Below is a list of changes made since GUDHI 3.3.0:

- [Alpha complex](https://gudhi.inria.fr/doc/latest/group__alpha__complex.html)
     - the C++ weighted version for alpha complex is now available in any dimension D.

- Simplex tree [C++](https://gudhi.inria.fr/doc/latest/class_gudhi_1_1_simplex__tree.html) [Python](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - A new method to reset the filtrations
     - A new method to get the boundaries of a simplex

- [Subsampling](https://gudhi.inria.fr/doc/latest/group__subsampling.html)
     - The C++ function `choose_n_farthest_points()` now takes a distance function instead of a kernel as first argument, users can replace `k` with `k.squared_distance_d_object()` in each call in their code.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.3.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.4.0+is%3Aclosed) is available on GitHub.


[Release 3.3.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.3.0)
-----------

Release date: August 2020

As a major new feature, the GUDHI library now offers a persistence-based clustering algorithm, weighted Rips complex using DTM
and edge collapse.

Below is a list of changes made since GUDHI 3.2.0:

- [DTM density estimator](https://gudhi.inria.fr/python/latest/point_cloud.html#module-gudhi.point_cloud.dtm)
     - Python implementation of a density estimator based on the distance to the empirical measure defined by a point set.

- [DTM Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_user.html#dtm-rips-complex)
     - This Python implementation constructs a weighted Rips complex giving larger weights to outliers, which reduces their impact on the persistence diagram

- [Alpha complex](https://gudhi.inria.fr/python/latest/alpha_complex_user.html) - Python interface improvements
     - 'fast' and 'exact' computations
     - Delaunay complex construction by not setting filtration values
     - Use the specific 3d alpha complex automatically to make the computations faster

- [Clustering](https://gudhi.inria.fr/python/latest/clustering.html)
     - Python implementation of [ToMATo](https://doi.org/10.1145/2535927), a persistence-based clustering algorithm

- [Edge Collapse](https://gudhi.inria.fr/doc/latest/group__edge__collapse.html) of a filtered flag complex
     - This C++ implementation reduces a filtration of Vietoris-Rips complex from its graph to another smaller flag filtration with the same persistence.

- [Bottleneck distance](https://gudhi.inria.fr/python/latest/bottleneck_distance_user.html)
     - Python interface to [hera](https://github.com/grey-narn/hera)'s bottleneck distance

- Persistence representations
     - [Atol](https://gudhi.inria.fr/python/latest/representations.html#gudhi.representations.vector_methods.Atol) is integrated in finite vectorisation methods. This [article](https://www.fujitsu.com/global/about/resources/news/press-releases/2020/0316-01.html) talks about applications using Atol. This module was originally available at [https://github.com/martinroyer/atol](https://github.com/martinroyer/atol)
     - Python interface change: [Wasserstein metrics](https://gudhi.inria.fr/python/latest/representations.html#gudhi.representations.metrics.WassersteinDistance) is now [hera](https://github.com/grey-narn/hera) by default

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.2.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.3.0+is%3Aclosed) is available on GitHub.


[Release 3.2.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.2.0)
-----------

Release date: May 2020

As a major new feature, the GUDHI library now offers Weighted Rips Complex, extended persistence and point cloud utilities new modules.

Below is a list of changes made since Gudhi 3.1.1:

- Point cloud utilities
     - A new module [Time Delay Embedding](https://gudhi.inria.fr/python/latest/point_cloud.html#time-delay-embedding) to embed time-series data in the R^d according to [Takens' Embedding Theorem](https://en.wikipedia.org/wiki/Takens%27s_theorem) and obtain the coordinates of each point.
     - A new module [K Nearest Neighbors](https://gudhi.inria.fr/python/latest/point_cloud.html#k-nearest-neighbors) that wraps several implementations for computing the k nearest neighbors in a point set.
     - A new module [Distance To Measure](https://gudhi.inria.fr/python/latest/point_cloud.html#distance-to-measure) to compute the distance to the empirical measure defined by a point set

- [Persistence representations](https://gudhi.inria.fr/python/latest/representations.html)
     - Interface to Wasserstein distances.

- Rips complex
     - A new module [Weighted Rips Complex](https://gudhi.inria.fr/python/latest/rips_complex_user.html#weighted-rips-complex) to construct a simplicial complex from a distance matrix and weights on vertices.

- [Wassertein distance](https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html)
     - An [another implementation](https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html#hera) comes from Hera (BSD-3-Clause) which is based on [Geometry Helps to Compare Persistence Diagrams](http://doi.acm.org/10.1145/3064175) by Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov.
     - `gudhi.wasserstein.wasserstein_distance` has now an option to return the optimal matching that achieves the distance between the two diagrams.
     - A new module [Barycenters](https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html#barycenters) to estimate the Frechet mean (aka Wasserstein barycenter) between persistence diagrams.

- [Simplex tree](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - Extend filtration method to compute extended persistence
     - Flag and lower star persistence pairs generators
     - A new interface to filtration, simplices and skeleton getters to return an iterator

- [Alpha complex](https://gudhi.inria.fr/doc/latest/group__alpha__complex.html)
     - Improve computations (cache circumcenters computation and point comparison improvement)

- [Persistence graphical tools](https://gudhi.inria.fr/python/latest/persistence_graphical_tools_user.html)
     - New rendering option proposed (use LaTeX style, add grey block, improved positioning of labels, etc.).
     - Can now handle (N x 2) numpy arrays as input

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.1.1](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.2.0+is%3Aclosed) is available on GitHub.


[Release 3.1.1](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.1.1)
-----------

Release date: February 2020

Gudhi-3.1.1 is a bug-fix release. In particular, it fixes the installation of the Python representation module.

The [list of bugs that were solved since gudhi-3.1.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.1.1+is%3Aclosed) is available on GitHub.


[Release 3.1.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.1.0)
-----------

Release date: January 2020

As a major new feature, the GUDHI library now offers 2 new Python modules: Persistence representations and Wasserstein distance.

Below is a list of changes made since Gudhi 3.0.0:

- [Persistence representations](https://gudhi.inria.fr/python/latest/representations.html) (new Python module)
     - Vectorizations, distances and kernels that work on persistence diagrams, compatible with scikit-learn. This module was originally available at https://github.com/MathieuCarriere/sklearn-tda and named sklearn_tda.

- [Wasserstein distance](https://gudhi.inria.fr/python/latest/wasserstein_distance_user.html) (new Python module)
     - The q-Wasserstein distance measures the similarity between two persistence diagrams.

- [Alpha complex](https://gudhi.inria.fr/doc/latest/group__alpha__complex.html) (new C++ interface)
     - Thanks to [CGAL 5.0 Epeck_d](https://doc.cgal.org/latest/Kernel_d/structCGAL_1_1Epeck__d.html) kernel, an exact computation version of Alpha complex dD is available and the default one (even in Python).

- [Persistence graphical tools](https://gudhi.inria.fr/python/latest/persistence_graphical_tools_user.html) (new Python interface)
     - Axes as a parameter allows the user to subplot graphics.
     - Use matplotlib default palette (can be user defined).

- Miscellaneous
     - Python `read_off` function has been renamed `read_points_from_off_file` as it only reads points from OFF files.
     - See the list of [bug fixes](https://github.com/GUDHI/gudhi-devel/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3A3.1.0+).


[Release 3.0.0](https://github.com/GUDHI/gudhi-devel/releases/tag/tags%2Fgudhi-release-3.0.0)
-----------

Release date: August 2019

As a major new feature, the GUDHI library is now released under a MIT license in order to ease the external contributions.

Below is a list of changes made since Gudhi 2.3.0:

- [Persistence graphical tools](https://gudhi.inria.fr/python/latest/persistence_graphical_tools_user.html) (new functionnality)
     - Add a persistence density graphical tool

- [Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_user.html) (new Python interface)
     - Sparse Rips complex is now available in Python.

- [Alpha complex](https://gudhi.inria.fr/doc/latest/group__alpha__complex.html) (new C++ interface)
     - Dedicated Alpha complex for 3d cases. Alpha complex 3d can be standard, weighted, periodic or weighted and periodic.

- Third parties (new dependencies)
     - C++14 is the new standard (instead of C++11 on former versions of GUDHI)
     - boost >= 1.56 is now required (instead of 1.48 on former versions of GUDHI)
     - CGAL >= 4.11 is now required (instead of various requirements on former versions of GUDHI)
     - Eigen >= 3.1.0 is now required (version was not checked)

Release 2.3.0
-----------

Release date: September 2018

As a major new feature, the GUDHI library now offers a Python interface
to the Nerve and Graph Induced Complex. The GUDHI conda package is now available through conda-forge channel.

Below is a list of changes made since Gudhi 2.2.0:
- [Nerve and Graph Induced Complex](https://gudhi.inria.fr/python/latest/nerve_gic_complex_user.html) (new Python interface)
     - Cover complexes, that provably contain topological information about the input data.

- Compilation issue (bug fix)
     - CMake warning with ninja generator.
     - thread_local management on old XCode versions.
     - Boost dependency for Windows Python module.

- [GUDHI conda package](https://gudhi.inria.fr/conda/)
     - The GUDHI conda package is now available through conda-forge channel.

Release 2.2.0
-----------

Release date: June 2018

As a major new feature, the GUDHI library now offers a Čech complex module, a sparse version of the Rips complex and
a utility to build the Rips complex from a correlation matrix (no Python interface yet).

Below is a list of changes made since Gudhi 2.1.0:

- [Čech complex](https://gudhi.inria.fr/doc/latest/group__cech__complex.html) (new package)
     - The Čech complex is a simplicial complex where the set of all simplices is filtered by the radius of their minimal enclosing ball.

- [Rips complex](https://gudhi.inria.fr/doc/latest/group__rips__complex.html) (new functions and interfaces)
     - A sparse version of the Rips complex
     - Rips complex from a correlation matrix utility

- [Dockerfile](https://gudhi.inria.fr/dockerfile/) (new installation process)
     A Dockerfile example is now provided to compile, test and install GUDHI in a container

- CGAL 4.12 compilation issue (bug fix)

- CMake minimal version is now 3.1
     - To take advantage of the latest features and simplify CMakeLists.txt files


Release 2.1.0
-----------

Release date: January 2018

As a major new feature, the GUDHI library now offers persistence representations and cover complex (no Python interface yet).

Below is a list of changes made since Gudhi 2.0.1:

- [Cover complex](https://gudhi.inria.fr/doc/latest/group__cover__complex.html) (new package)
     - Nerves and Graph Induced Complexes are cover complexes, that provably contain topological information about the input data.

- [Representations of persistence diagrams](https://gudhi.inria.fr/doc/latest/group___persistence__representations.html) (new package)
     - It contains implementation of various representations of persistence
     diagrams. It implements basic functionalities which are neccessary to use
     persistence in statistics and machine learning.

- [Simplex tree](https://gudhi.inria.fr/doc/latest/group__simplex__tree.html) (new functions and interfaces)
     - Graph expansion with a blocker oracle.
     - Cech complex implementation example using CGAL mini spheres in fixed dimension.
     - Automatic dimension set mechanism.

- Alpha complex (new function)
     - [Weighted periodic 3D](https://gudhi.inria.fr/doc/latest/_alpha_complex_2weighted_periodic_alpha_complex_3d_persistence_8cpp-example.html) version utility.

- CGAL 4.11 compilation issue (bug fix)

- [Cubical complex](https://gudhi.inria.fr/doc/latest/group__cubical__complex.html) (bug fix)
     - Perseus file read function.
     - Computations of incidence indices between cubes.
     - Missing periodic argument for the Python version.

- Documentation
     - New [file formats](https://gudhi.inria.fr/doc/latest/fileformats.html) section.
     - Bugs fix

- [Utilities](https://gudhi.inria.fr/utils)
     - Separate examples from utilities.

- GUDHI Debian package is available for Debian Testing distribution.

Release 2.0.1
-----------

Release date: September 2017

This minor GUDHI library version is fixing issues and improves persistence graphical tools module. All the new modules comes with their Python interface and a lot of examples, even in Python.

Below is a list of changes made since Gudhi 2.0.0:

- Spatial searching (new function)

     - Spatial searching is now offering a radius search method.

- Persistence graphical tools (interfaces improvement)

     - Add a band boot display mechanism on persistence diagram.
     - Number of points and barcode limitation before the display.
     - Interface modification to read persistence from a file.

- Bottleneck distance (bug fix)

     - Read persistence files with infinity values bug is fixed.

- Persistent cohomology (bug fix)

     - Weighted Alpha complex 3d persistence bug is fixed.

- Simplex tree (dead code)

     - Remove useless global filtration attribute, getter and setter.

- cython

     - Plot persistence functions improvement.
     - Rename cythonize_gudhi.py in setup.py to be conform with Python conventions.
     - Windows python module compilation issue fix.
     - Documentation generation bug fix.
     - Reduce CMake interactions.
     - Rips complex memory leak fix.

- Doxygen

     - Using MathJax.js to generate LaTeX instead of png.

- CMake

     - Boost dependencies improvement.
     - Conda compilation issues fix.
     - Modules activation/desactivation mechanism for compilation and test.

- Data points generator

     - Move data/points/generator into src/common/utilities.

Release 2.0.0
-----------

Release date: April 2017

As a major new feature, the GUDHI library now offers an interface with Python. All the new modules comes with their Python interface and a lot of examples, even in Python.

Below is a list of changes made since Gudhi 1.3.1:

- Bottleneck distance (new package)

     - Bottleneck distance measures the similarity between two persistence diagrams.

- cython (new package)

     - A Cython package allows to compile a Python interface with the GUDHI library.

- Spatial searching (new package)

     - Spatial searching is a wrapper around [CGAL dD spatial searching](http://doc.cgal.org/latest/Spatial_searching/index.html)
     algorithms that provides a simplified API to perform (approximate) neighbor searches.

- Subsampling (new package)

     - Subsampling offers methods to subsample a set of points.

- Tangential complex (new package)

     - A Tangential Delaunay complex is a non filtered simplicial complex designed to
     reconstruct a k-dimensional manifold embedded in d-dimensional Euclidean space.

- Witness complex (new relaxed version with a new interface)

     - Witness complex Wit(W,L) is a simplicial complex defined on two sets of
     points in Rd. The new relaxed version is filtrated. The new interface eases the
     simplicial complex construction.

- Alpha complex (new interface)

     - Alpha complex is a simplicial complex data structure constructed from the
     finite cells of a Delaunay Triangulation. The new interface eases the simplicial
     complex construction.

- Rips complex (new interface)

     - The Rips complex is a simplicial complex constructed from an expanded one-skeleton
     graph. The new interface eases the simplicial complex construction and allows the Rips
     complex to be build from a distance matrix.


Release 1.3.1
-----------

Release date: September 2016

Below is a list of changes made since Gudhi 1.3.0:

- Bug fix
As Simplex_handle default type was an 'int' in the simplex tree data structure, persistence cohomology computation was segmentation faulting
when the number of simplices was going further to the 'int' maximum value (2.7 billions of simplices on a classical modern machine and OS).
Simplex_handle default type is now an 'std::uint32_t' and can go up to about 4 billions of simplices.

- CMake
     - The messages from CMake has been rewritten to be more consistent.

- Coding conventions
     - CMake projects names and C++ namespaces have been homogenized.

- Documentation
     - Mandatory and optional third party libraries have been separated in the documentation.

- Data sets
     - in data/points/generator : thanks to [Aurélien Alvarez](http://www.aurelienalvarez.org/), aurelien_alvarez_surfaces_in_R8.py is
     a script to generate points on a surface in R8.


Release 1.3.0
-----------

Release date: April 2016

Below is a list of changes made since Gudhi 1.2.0:

- Alpha complex (new package)
     - Alpha_complex is a simplicial complex data structure constructed from the
     finite cells of a Delaunay Triangulation.

- Cubical complex (new package)
     - The cubical complex is an example of a structured complex useful in
     computational mathematics (specially rigorous numerics) and image
     analysis.

- Witness complex (new package)
     - Witness complex Wit(W,L) is a simplicial complex defined on two sets of
     points in Rd. The data structure is described in the paper "Jean-Daniel
     Boissonnat and Clément Maria. The Simplex Tree: An Efficient Data
     Structure for General Simplicial Complexes. Algorithmica, pages 1–22,
     2014."

- Persistent cohomology (new examples)
     - alpha_complex_persistence: How to compute persistent homology from an
       alpha complex.
     - periodic_alpha_complex_3d_persistence.cpp: How to compute persistent
       homology on a periodic alpha complex in 3D.
     - Bitmap_cubical_complex.cpp: How to compute persistent homology from an
       cubical complex.
     - Bitmap_cubical_complex_periodic_boundary_conditions.cpp: How to compute
       persistent homology from a periodic cubical complex.

- Documentation
     - a section on examples.
     - a summary of modules on the main page.

- Data sets
     - tore3D_300.off: 300 random points on a 3D torus.
     - grid_10_10_10_in_0_1.off: Points every 0.1 in each 3 dimensions in
       [0, 1].
     - in data/bitmaps : examples of Perseus style file for Cubical complex

Release 1.2.0
-----------

Release date: November 2015

Below is a list of changes made since GUDHI 1.1.0:

- GudhUI (new package)
This package provides a User Interface to have a quick overview of GUDHI features.
GudhUI allows to load OFF points or meshes files that can be visualized on a 3D viewer.
A tool suite is available to perform edge contraction, point cloud persistence, and so on.

- Persistent cohomology (new examples)
     - persistence_from_simple_simplex_tree: How to compute persistent homology after inserting simplices and their filtration values
     - plain_homology: How to compute (non-persistent) homology on simplices
     - alpha_shapes_persistence.cpp: How to compute persistent homology after constructing a CGAL 3D alpha shape simplices on a point cloud file

- Simplex tree (new example)
     - mini_simplex_tree: How to fill a minimalist (in terms of memory usage) simplex tree

- Skeleton blocker (new examples)
     - Skeleton_blocker_from_simplices: How to fill a skeleton blocker with simplices
     - Skeleton_blocker_link: How to find a link in a skeleton blocker

- Compilation
     - Fix CMake issues.
     - Fix Clang compilation.
     - Fix Visual Studio 2013 compilation.

- Data sets
     - spiral_3d_10k.off: 10000 points on a 3D spiral.
     - spiral_4d_10k.off: 10000 points on a 4D spiral.
     - tore3D_1307.off: 1307 points on a 3D torus.
     - images: pictures taken from a camera rotating around an object (asian mug, little car, lucky cat and money pig) and their off files associated.

Release 1.1.0
-----------

Release date: December 2014

This version contains :

- Skeleton blocker data structure
- Edge contraction
- Persistence from CGAL alpha shapes in 3D

Release 1.0.2
-----------

Release date: June 2014

First publicly available release of the GUDHI library.
This version contains :

- Simplex tree data structure
- Persistence cohomology
- Persistence from Rips
