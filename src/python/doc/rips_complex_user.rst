:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Rips complex user manual
=========================
Definition
----------

====================================================================  ================================  ======================
:Authors: Clément Maria, Pawel Dlotko, Vincent Rouvreau, Marc Glisse  :Since: GUDHI 2.0.0               :License: GPL v3
====================================================================  ================================  ======================

+-------------------------------------------+----------------------------------------------------------------------+
| :doc:`rips_complex_user`                  | :doc:`rips_complex_ref`                                              |
+-------------------------------------------+----------------------------------------------------------------------+

The `Rips complex <https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex>`_ is a simplicial complex that
generalizes proximity (:math:`\varepsilon`-ball) graphs to higher dimensions. The vertices correspond to the input
points, and a simplex is present if and only if its diameter is smaller than some parameter α.  Considering all
parameters α defines a filtered simplicial complex, where the filtration value of a simplex is its diameter.
The filtration can be restricted to values α smaller than some threshold, to reduce its size.  Beware that some
people define the Rips complex using a bound of 2α instead of α, particularly when comparing it to an ambient
Čech complex.  They end up with the same combinatorial object, but filtration values which are half of ours.

The input discrete metric space can be provided as a point cloud plus a distance function, or as a distance matrix.

When creating a simplicial complex from the graph, :doc:`RipsComplex <rips_complex_ref>` first builds the graph and
inserts it into the data structure. It then expands the simplicial complex (adds the simplices corresponding to cliques)
when required. The expansion can be stopped at dimension `max_dimension`, by default 1.

A vertex name corresponds to the index of the point in the given range (aka. the point cloud).

.. figure::
    ../../doc/Rips_complex/rips_complex_representation.png
    :align: center

    Rips-complex one skeleton graph representation

On this example, as edges (4,5), (4,6) and (5,6) are in the complex, simplex (4,5,6) is added with the filtration value
set with :math:`max(filtration(4,5), filtration(4,6), filtration(5,6))`. And so on for simplex (0,1,2,3).

If the :doc:`RipsComplex <rips_complex_ref>` interfaces are not detailed enough for your need, please refer to
rips_persistence_step_by_step.cpp C++ example, where the graph construction over the Simplex_tree is more detailed.

A Rips complex can easily become huge, even if we limit the length of the edges
and the dimension of the simplices. One easy trick, before building a Rips
complex on a point cloud, is to call :func:`~gudhi.sparsify_point_set` which removes points
that are too close to each other. This does not change its persistence diagram
by more than the length used to define "too close".

A more general technique is to use a sparse approximation of the Rips
introduced by Don Sheehy :cite:`sheehy13linear`. We are using the version
described in :cite:`buchet16efficient` (except that we multiply all filtration
values by 2, to match the usual Rips complex). :cite:`cavanna15geometric` proves
a :math:`\frac{1}{1-\varepsilon}`-interleaving, although in practice the
error is usually smaller.  A more intuitive presentation of the idea is
available in :cite:`cavanna15geometric`, and in a video
:cite:`cavanna15visualizing`. Passing an extra argument `sparse=0.3` at the
construction of a :class:`~gudhi.RipsComplex` object asks it to build a sparse Rips with
parameter :math:`\varepsilon=0.3`, while the default `sparse=None` builds the
regular Rips complex.


Point cloud
-----------

Example from a point cloud
^^^^^^^^^^^^^^^^^^^^^^^^^^

This example builds the neighborhood graph from the given points, up to max_edge_length.
Then it creates a :doc:`SimplexTree <simplex_tree_ref>` with it.

Finally, it is asked to display information about the simplicial complex.

.. testcode::

    import gudhi
    rips_complex = gudhi.RipsComplex(points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]],
                                     max_edge_length=12.0)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    fmt = '%s -> %.2f'
    for filtered_value in simplex_tree.get_filtration():
        print(fmt % tuple(filtered_value))

When launching (Rips maximal distance between 2 points is 12.0, is expanded
until dimension 1 - one skeleton graph in other words), the output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    [0] -> 0.00
    [1] -> 0.00
    [2] -> 0.00
    [3] -> 0.00
    [4] -> 0.00
    [5] -> 0.00
    [6] -> 0.00
    [2, 3] -> 5.00
    [4, 5] -> 5.39
    [0, 2] -> 5.83
    [0, 1] -> 6.08
    [1, 3] -> 6.32
    [1, 2] -> 6.71
    [5, 6] -> 7.28
    [2, 4] -> 8.94
    [0, 3] -> 9.43
    [4, 6] -> 9.49
    [3, 6] -> 11.00

Notice that if we use

.. code-block:: python

    rips_complex = gudhi.RipsComplex(points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]],
                                     max_edge_length=12.0, sparse=2)

asking for a very sparse version (theory only gives some guarantee on the meaning of the output if `sparse<1`),
2 to 5 edges disappear, depending on the random vertex used to start the sparsification.

Example from OFF file
^^^^^^^^^^^^^^^^^^^^^

This example builds the :doc:`RipsComplex <rips_complex_ref>` from the given
points in an OFF file, and max_edge_length value.
Then it creates a :doc:`SimplexTree <simplex_tree_ref>` with it.

Finally, it is asked to display information about the Rips complex.


.. testcode::

    import gudhi
    off_file = gudhi.__root_source_dir__ + '/data/points/alphacomplexdoc.off'
    point_cloud = gudhi.read_points_from_off_file(off_file = off_file)
    rips_complex = gudhi.RipsComplex(points=point_cloud, max_edge_length=12.0)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    fmt = '%s -> %.2f'
    for filtered_value in simplex_tree.get_filtration():
        print(fmt % tuple(filtered_value))

the program output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    [0] -> 0.00
    [1] -> 0.00
    [2] -> 0.00
    [3] -> 0.00
    [4] -> 0.00
    [5] -> 0.00
    [6] -> 0.00
    [2, 3] -> 5.00
    [4, 5] -> 5.39
    [0, 2] -> 5.83
    [0, 1] -> 6.08
    [1, 3] -> 6.32
    [1, 2] -> 6.71
    [5, 6] -> 7.28
    [2, 4] -> 8.94
    [0, 3] -> 9.43
    [4, 6] -> 9.49
    [3, 6] -> 11.00

Distance matrix
---------------

Example from a distance matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example builds the one skeleton graph from the given distance matrix, and max_edge_length value.
Then it creates a :doc:`SimplexTree <simplex_tree_ref>` with it.

Finally, it is asked to display information about the simplicial complex.

.. testcode::

    import gudhi
    rips_complex = gudhi.RipsComplex(distance_matrix=[[],
                                                      [6.0827625303],
                                                      [5.8309518948, 6.7082039325],
                                                      [9.4339811321, 6.3245553203, 5],
                                                      [13.0384048104, 15.6524758425, 8.94427191, 12.0415945788],
                                                      [18.0277563773, 19.6468827044, 13.152946438, 14.7648230602, 5.3851648071],
                                                      [17.88854382, 17.1172427686, 12.0830459736, 11, 9.4868329805, 7.2801098893]],
                                     max_edge_length=12.0)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    fmt = '%s -> %.2f'
    for filtered_value in simplex_tree.get_filtration():
        print(fmt % tuple(filtered_value))

When launching (Rips maximal distance between 2 points is 12.0, is expanded
until dimension 1 - one skeleton graph in other words), the output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    [0] -> 0.00
    [1] -> 0.00
    [2] -> 0.00
    [3] -> 0.00
    [4] -> 0.00
    [5] -> 0.00
    [6] -> 0.00
    [2, 3] -> 5.00
    [4, 5] -> 5.39
    [0, 2] -> 5.83
    [0, 1] -> 6.08
    [1, 3] -> 6.32
    [1, 2] -> 6.71
    [5, 6] -> 7.28
    [2, 4] -> 8.94
    [0, 3] -> 9.43
    [4, 6] -> 9.49
    [3, 6] -> 11.00

Example from csv file
^^^^^^^^^^^^^^^^^^^^^

This example builds the :doc:`RipsComplex <rips_complex_ref>` from the given
distance matrix in a csv file, and max_edge_length value.
Then it creates a :doc:`SimplexTree <simplex_tree_ref>` with it.

Finally, it is asked to display information about the Rips complex.


.. testcode::

    import gudhi
    distance_matrix = gudhi.read_lower_triangular_matrix_from_csv_file(csv_file=gudhi.__root_source_dir__ + \
        '/data/distance_matrix/full_square_distance_matrix.csv')
    rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=12.0)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    fmt = '%s -> %.2f'
    for filtered_value in simplex_tree.get_filtration():
        print(fmt % tuple(filtered_value))

the program output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    [0] -> 0.00
    [1] -> 0.00
    [2] -> 0.00
    [3] -> 0.00
    [4] -> 0.00
    [5] -> 0.00
    [6] -> 0.00
    [2, 3] -> 5.00
    [4, 5] -> 5.39
    [0, 2] -> 5.83
    [0, 1] -> 6.08
    [1, 3] -> 6.32
    [1, 2] -> 6.71
    [5, 6] -> 7.28
    [2, 4] -> 8.94
    [0, 3] -> 9.43
    [4, 6] -> 9.49
    [3, 6] -> 11.00

Correlation matrix
------------------

Example from  a correlation matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Analogously to the case of distance matrix, Rips complexes can be also constructed based on correlation matrix.
Given a correlation matrix M, comportment-wise 1-M is a distance matrix.
This example builds the one skeleton graph from the given corelation matrix and threshold value.
Then it creates a :doc:`SimplexTree <simplex_tree_ref>` with it.

Finally, it is asked to display information about the simplicial complex.

.. testcode::

    import gudhi
    import numpy as np

    # User defined correlation matrix is:
    # |1     0.06    0.23    0.01    0.89|
    # |0.06  1       0.74    0.01    0.61|
    # |0.23  0.74    1       0.72    0.03|
    # |0.01  0.01    0.72    1       0.7 |
    # |0.89  0.61    0.03    0.7     1   |
    correlation_matrix=np.array([[1., 0.06, 0.23, 0.01, 0.89],
                                [0.06, 1., 0.74, 0.01, 0.61],
                                [0.23, 0.74, 1., 0.72, 0.03],
                                [0.01, 0.01, 0.72, 1., 0.7],
                                [0.89, 0.61, 0.03, 0.7, 1.]], float)

    distance_matrix = 1 - correlation_matrix
    rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    fmt = '%s -> %.2f'
    for filtered_value in simplex_tree.get_filtration():
        print(fmt % tuple(filtered_value))

When launching (Rips maximal distance between 2 points is 12.0, is expanded
until dimension 1 - one skeleton graph in other words), the output is:

.. testoutput::

    Rips complex is of dimension 1 - 15 simplices - 5 vertices.
    [0] -> 0.00
    [1] -> 0.00
    [2] -> 0.00
    [3] -> 0.00
    [4] -> 0.00
    [0, 4] -> 0.11
    [1, 2] -> 0.26
    [2, 3] -> 0.28
    [3, 4] -> 0.30
    [1, 4] -> 0.39
    [0, 2] -> 0.77
    [0, 1] -> 0.94
    [2, 4] -> 0.97
    [0, 3] -> 0.99
    [1, 3] -> 0.99

.. note::
    If you compute the persistence diagram and convert distances back to correlation values,
    points in the persistence diagram will be under the diagonal, and
    bottleneck distance and persistence graphical tool will not work properly,
    this is a known issue.
