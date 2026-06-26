We are pleased to announce the release 3.13.0 of the GUDHI library.

The GUDHI library is mainly developped using GitHub, do not hesitate to
[fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel).
From a user point of view, we recommend to download GUDHI user version (gudhi.3.X.X.tar.gz).

Below is a list of changes:

- Simplex tree [`[C++]`](https://gudhi.inria.fr/doc/latest/class_gudhi_1_1_simplex__tree.html) and [`[Python`](https://gudhi.inria.fr/python/latest/simplex_tree_ref.html)
     - A new `euler_characteristic` method, and a Python binding to `num_simplices_by_dimension`.
     - Preliminary support for multiple filtration values stored in the Simplex tree (available only in `C++`).
     - Serialization has been upgraded for multiple filtration values support, but also for a `SERIALIZATION_VERSION` number in order to track modification and compatibility between serialization files. It means the serialization files, including pickled SimplexTree, generated from a pre-3.13.0 gudhi version won't be supported (an exception is thrown).

- [Reproducibility](https://gudhi.inria.fr/doc/latest/group__reproducibility.html) `[C++]`
     - New random functions, with a default random generator that is used internally, and where the user can set the global seed.

- [Tangential complex](https://gudhi.inria.fr/python/latest/tangential_complex_user.html) `[Python]`
     - `intrisic_dim` constructor argument is deprecated (typo), please use `intrinsic_dim` instead.

- [Module](link)
     - **...**
     
- [Module](link)
     - **...**
     
- Miscellaneous
     - The [list of bugs that were solved](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.13.0+is%3Aclosed) is available on GitHub.

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

- @hschreiber
- @mglisse
- @ryan-charette
- @VincentRouvreau
