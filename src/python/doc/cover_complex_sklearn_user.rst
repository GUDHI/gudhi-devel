:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Mapper/GIC/Nerve complexes user manual
======================================

Definition
----------

We provide two types of classes for computing cover complexes: one is comprised of the :class:`~gudhi.cover_complex.NerveComplex`, :class:`~gudhi.cover_complex.GraphInducedComplex` and :class:`~gudhi.cover_complex.MapperComplex` classes, which correspond to the Nerve, Graph Induced and Mapper complexes respectively, and are written in a scikit-learn format,
while the other type is :class:`~gudhi.CoverComplex`, which only computes Nerve and Graph Induced complexes, but is a bit more flexible for input types (it can read inputs from paths to files for instance,
while :class:`~gudhi.cover_complex.NerveComplex`, :class:`~gudhi.cover_complex.GraphInducedComplex` and :class:`~gudhi.cover_complex.MapperComplex` need the inputs to be stored in memory).

Key differences between Mapper, Nerve and Graph Induced complexes (GIC) are: Mapper nodes are defined with given input clustering method while GIC nodes are defined with given input graph and Nerve nodes are defined with cover elements, GIC accepts partitions instead of covers while Mapper and Nerve require cover elements to overlap. Also, note that when the cover is functional (i.e., preimages of filter functions), GIC only accepts one scalar-valued filter with gain < 0.5. On the other hand, Mapper complexes accept resolutions and gains with any length.

These classes can print output files, which can then be visualized with either
neato (from `graphviz <http://www.graphviz.org/>`_),
`geomview <http://www.geomview.org/>`_,
`KeplerMapper <https://github.com/scikit-tda/kepler-mapper>`_.

Covers
------

Mapper, Nerve and Graph Induced complexes require a cover C of the input point cloud P,
that is, a set of subsets of P whose union is P itself.
Very often, this cover is obtained from the preimage of a family of hypercubes covering
the image of some function f sending P to Euclidean space. This family is parameterized
by its resolution, which is the number of hypercubes,
and its gain, which is the overlap percentage between consecutive hypercubes (ordered by their lower corners).

Class for cover complexes
=========================

.. list-table::
   :widths: 40 30 30
   :header-rows: 0

   * - :Since: GUDHI 3.5.0
     - :License: MIT
     - :Requires: `Scikit-learn <installation.html#scikit-learn>`__

We provide classes for computing cover complexes, i.e., Nerve, Graph Induced and Mapper complexes. Detailed examples on how to use these classes in practice are available
in the following `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-cover-complex.ipynb>`_.

Example of Mapper cover complex computed from a point cloud
-----------------------------------------------------------

.. testcode::

    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    from gudhi.cover_complex import MapperComplex

    X = np.array([[1,1],[1,1.5],[1,2],[1,2.5],[1,3],[1.5,2],[1.5,3],[2,2],[2,2.5],[2,3],[2,3.5],[2,4]])
    F = np.array([[1,1,1,1,1,1.5,1.5,2,2,2,2,2],[1,1.5,2,2.5,3,2,3,2,2.5,3,3.5,4]]).T


    mapper = MapperComplex(
        input_type="point cloud",
        filter_bnds=np.array([[0.5, 2.5], [0.5, 4.5]]),
        resolutions=np.array([2, 4]),
        gains=np.array([0.3, 0.3]),
        clustering=AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=0.6),
    )

    mapper.fit(X, filters=F)

    print([s for s,_ in mapper.simplex_tree_.get_simplices()])

.. testoutput::

    [[0, 1], [0], [1, 2], [1, 3], [1], [2, 4], [2], [3, 4], [3], [4, 5], [4], [5]]
