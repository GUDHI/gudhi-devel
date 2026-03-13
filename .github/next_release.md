We are pleased to announce the release 3.12.0 of the GUDHI library.

As a major new feature, the GUDHI library now offers  **...**

The GUDHI library is mainly developped using GitHub, do not hesitate to
[fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel).
From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).

Below is a list of changes:

- [Module](link)
     - **...**

- [Module](link)
     - **...**

- [Alpha complex dD](https://gudhi.inria.fr/doc/latest/class_gudhi_1_1alpha__complex_1_1_alpha__complex.html)
     - **API break:** The simplicial complex for the Alpha complex concept has been changed to
       `dimension_simplex_range` (that must returns a range of simplices of a given dimension) instead of
       `skeleton_simplex_range` (that was returning a range of simplices lower or equal to a given dimension)

- [Simplex_tree](https://gudhi.inria.fr/doc/latest/class_gudhi_1_1_simplex__tree.html)
     - A new iterator over the simplices of the simplicial complex that match a given dimension

- [Representations](https://gudhi.inria.fr/python/latest/representations.html)
     - in metrics, `BottleneckDistance` argument `epsilon` is deprecrated and renamed `e` to be consistent with `bottleneck_distance` and `pairwise_persistence_diagram_distances`

- Installation
     - Pip package is now available for OSx &ge; 13.0 (was &ge; 12.0).
     - Support Eigen 3.4.1 and 5.X.X
     
- Miscellaneous
     - The [list of bugs that were solved](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.12.0+is%3Aclosed) is available on GitHub.

All modules are distributed under the terms of the MIT license.
However, there are still GPL dependencies for many modules. We invite you to check our [license dedicated web page](https://gudhi.inria.fr/licensing/) for further details.

We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use of
the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.

We provide [bibtex entries](https://gudhi.inria.fr/doc/latest/_citation.html) for the modules of the User and Reference
Manual, as well as for publications directly related to the GUDHI library.

Feel free to [contact us](https://gudhi.inria.fr/contact/) in case you have any questions or remarks.

For further information about downloading and installing the library
([C++](https://gudhi.inria.fr/doc/latest/installation.html) or
[Python](https://gudhi.inria.fr/python/latest/installation.html)), please visit the
[GUDHI web site](https://gudhi.inria.fr/).

## Contributors

- **...**
- **...**
