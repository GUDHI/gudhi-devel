We are pleased to announce the release 3.9.0 of the GUDHI library.

We are now using GitHub to develop the GUDHI library, do not hesitate to [fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel). From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).

Below is a list of changes made since GUDHI 3.8.0:

- [CubicalPersistence](https://gudhi.inria.fr/python/latest/cubical_complex_sklearn_itf_ref.html)
     - Much faster implementation for the 2d case with input from top-dimensional cells.

- [Simplex_tree](https://gudhi.inria.fr/doc/latest/group__simplex__tree.html)
     - A helper `for_each_simplex` that applies a given function object on each simplex
     - `make_filtration_non_decreasing` has been rewritten, and a new method `num_simplices_by_dimension` is now available thanks to this helper.
     - A `clear` method to empty the data stucture.
     - A new argument `ignore_infinite_values` for `initialize_filtration` method to skip infinite values. As a side effect, this change enhances the persistence computation.
     - `Simplex_tree_options_full_featured` has been renamed `Simplex_tree_options_default` and `Simplex_tree_options_python`.
     These are respectively the options used by all the utilities and by the python interface of the `SimplexTree` (as before this version).
     - From GUDHI 3.9.0, `Simplex_tree_options_full_featured` now activates `link_nodes_by_label` and `stable_simplex_handles`.

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
     - A Python function `collapse_edges` to simplify a clique filtration (represented as a sparse weighted graph), while preserving its persistent homology.

- [Mapper/GIC/Nerve complexes](https://gudhi.inria.fr/python/latest/cover_complex_sklearn_isk_ref.html)
     - A new method `save_to_html` to ease the Keppler Mapper visualization

- Installation
     - Boost &ge; 1.71.0 is now required (was &ge; 1.66.0).
     - cython >= 3.0.0 is now supported.
     - Python 3.12 pip package.

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.8.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.9.0+is%3Aclosed) is available on GitHub.

All modules are distributed under the terms of the MIT license.
However, there are still GPL dependencies for many modules. We invite you to check our [license dedicated web page](https://gudhi.inria.fr/licensing/) for further details.

We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use of the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.

We provide [bibtex entries](https://gudhi.inria.fr/doc/latest/_citation.html) for the modules of the User and Reference Manual, as well as for publications directly related to the GUDHI library. 

Feel free to [contact us](https://gudhi.inria.fr/contact/) in case you have any questions or remarks.

For further information about downloading and installing the library ([C++](https://gudhi.inria.fr/doc/latest/installation.html) or [Python](https://gudhi.inria.fr/python/latest/installation.html)), please visit the [GUDHI web site](https://gudhi.inria.fr/).

## Contributors

- @DavidLapous
- @hschreiber
- @MathieuCarriere
- @martinroyer
- @mglisse
- @VincentRouvreau
