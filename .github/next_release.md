We are pleased to announce the release 3.8.0 of the GUDHI library.

As a major new feature, the GUDHI library now offers Perslay, a Tensorflow model for the representations module, scikit-learn like interfaces for Cover Complexes, a new function to compute persistence of a function on ℝ and the possibility to build a Cubical Complex as a lower-star filtration from vertices.

We are now using GitHub to develop the GUDHI library, do not hesitate to [fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel). From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).

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

- Installation
     - CMake 3.8 is the new minimal standard to compile the library.
     - Support for oneAPI TBB (instead of deprecated TBB) to take advantage of multicore performance.

- [Python documentation](https://gudhi.inria.fr/python/latest/installation.html)
     - [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/) is the new sphinx theme of the python documentation.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.7.1](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.8.0+is%3Aclosed) is available on GitHub.

All modules are distributed under the terms of the MIT license.
However, there are still GPL dependencies for many modules. We invite you to check our [license dedicated web page](https://gudhi.inria.fr/licensing/) for further details.

We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use of the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.

We provide [bibtex entries](https://gudhi.inria.fr/doc/latest/_citation.html) for the modules of the User and Reference Manual, as well as for publications directly related to the GUDHI library. 

Feel free to [contact us](https://gudhi.inria.fr/contact/) in case you have any questions or remarks.

For further information about downloading and installing the library ([C++](https://gudhi.inria.fr/doc/latest/installation.html) or [Python](https://gudhi.inria.fr/python/latest/installation.html)), please visit the [GUDHI web site](https://gudhi.inria.fr/).

## Contributors

- @Hind-M
- @MathieuCarriere
- @mglisse
- @wreise
- @VincentRouvreau
