We are pleased to announce the release 3.X.X of the GUDHI library.

As a major new feature, the GUDHI library now offers ...

We are now using GitHub to develop the GUDHI library, do not hesitate to [fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel). From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).

Below is a list of changes made since GUDHI 3.5.0:

- [Alpha complex](https://gudhi.inria.fr/python/latest/alpha_complex_user.html)
     - the python weighted version for alpha complex is now available in any dimension D.
     - `alpha_complex = gudhi.AlphaComplex(off_file='/data/points/tore3D_300.off')` is deprecated, please use [read_points_from_off_file](https://gudhi.inria.fr/python/latest/point_cloud.html#gudhi.read_points_from_off_file) instead.

- [Representations](https://gudhi.inria.fr/python/latest/representations.html#gudhi.representations.vector_methods.BettiCurve)
     - A more flexible Betti curve class capable of computing exact curves

- [Simplex tree](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - `__deepcopy__`, `copy` and copy constructors

- Installation
     - Boost &ge; 1.66.0 is now required (was &ge; 1.56.0).
     - Python >= 3.5 and cython >= 0.27 are now required.

- [Module](link)
     - ...

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.5.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.6.0+is%3Aclosed) is available on GitHub.

All modules are distributed under the terms of the MIT license.
However, there are still GPL dependencies for many modules. We invite you to check our [license dedicated web page](https://gudhi.inria.fr/licensing/) for further details.

We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use of the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.

We provide [bibtex entries](https://gudhi.inria.fr/doc/latest/_citation.html) for the modules of the User and Reference Manual, as well as for publications directly related to the GUDHI library. 

Feel free to [contact us](https://gudhi.inria.fr/contact/) in case you have any questions or remarks.

For further information about downloading and installing the library ([C++](https://gudhi.inria.fr/doc/latest/installation.html) or [Python](https://gudhi.inria.fr/python/latest/installation.html)), please visit the [GUDHI web site](https://gudhi.inria.fr/).

