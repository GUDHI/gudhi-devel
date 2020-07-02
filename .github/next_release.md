We are pleased to announce the release 3.3.0 of the GUDHI library.

As a major new feature, the GUDHI library now offers a persistence-based clustering algorithm and weighted Rips complex using DTM.

We are now using GitHub to develop the GUDHI library, do not hesitate to [fork the GUDHI project on GitHub](https://github.com/GUDHI/gudhi-devel). From a user point of view, we recommend to download GUDHI user version (gudhi.3.3.0.tar.gz).

Below is a list of changes made since GUDHI 3.2.0:

- [DTM density estimator](https://gudhi.inria.fr/python/latest/point_cloud.html#module-gudhi.point_cloud.dtm)
     - Python implementation of a density estimator based on the distance to the empirical measure defined by a point set.

- [DTM Rips complex](https://gudhi.inria.fr/python/latest/rips_complex_user.html#dtm-rips-complex)
     - This Python implementation constructs a weighted Rips complex giving larger weights to outliers,
     which reduces their impact on the persistence diagram

- [Alpha complex](https://gudhi.inria.fr/python/latest/alpha_complex_user.html) - Python interface improvements
     - 'fast' and 'exact' computations
     - Delaunay complex construction by not setting filtration values
     - Use the specific 3d alpha complex automatically to make the computations faster

- [Clustering](https://gudhi.inria.fr/python/latest/clustering.html)
     - Python implementation of [ToMATo](https://doi.org/10.1145/2535927), a persistence-based clustering algorithm

- [Bottleneck distance](https://gudhi.inria.fr/python/latest/bottleneck_distance_user.html)
     - Python interface to [hera](https://github.com/grey-narn/hera)'s bottleneck distance

- Representations - Python interface change
     - [Wasserstein metrics](https://gudhi.inria.fr/python/latest/representations.html#gudhi.representations.metrics.WassersteinDistance)
     is now [hera](https://github.com/grey-narn/hera) by default

- Miscellaneous
     - The [list of bugs that were solved since GUDHI-3.2.0](https://github.com/GUDHI/gudhi-devel/issues?q=label%3A3.3.0+is%3Aclosed)
     is available on GitHub.

All modules are distributed under the terms of the MIT license.
However, there are still GPL dependencies for many modules. We invite you to check our [license dedicated web page](https://gudhi.inria.fr/licensing/) for further details.

We kindly ask users to cite the GUDHI library as appropriately as possible in their papers, and to mention the use of the GUDHI library on the web pages of their projects using GUDHI and provide us with links to these web pages.

We provide [bibtex entries](https://gudhi.inria.fr/doc/latest/_citation.html) for the modules of the User and Reference Manual, as well as for publications directly related to the GUDHI library. 

Feel free to [contact us](https://gudhi.inria.fr/contact/) in case you have any questions or remarks.

For further information about downloading and installing the library ([C++](https://gudhi.inria.fr/doc/latest/installation.html) or [Python](https://gudhi.inria.fr/python/latest/installation.html)), please visit the [GUDHI web site](https://gudhi.inria.fr/).

