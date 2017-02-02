=========================
Rips complex user manual
=========================
Definition
----------

=======================================================  =====================================  =====================================
:Authors: Clément Maria, Pawel Dlotko, Vincent Rouvreau  :Introduced in: GUDHI 2.0.0            :Copyright: GPL v3
=======================================================  =====================================  =====================================

+-------------------------------------------+----------------------------------------------------------------------+
| :doc:`rips_complex_user`                  | :doc:`rips_complex_ref`                                              |
+-------------------------------------------+----------------------------------------------------------------------+

`Rips complex <https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex>`_ is a one skeleton graph that allows to
construct a simplicial complex from it. The input can be a point cloud with a given distance function, or a distance
matrix.

The filtration value of each edge is computed from a user-given distance function, or directly from the distance
matrix.

All edges that have a filtration value strictly greater than a given threshold value are not inserted into the complex.

When creating a simplicial complex from this one skeleton graph, Rips inserts the one skeleton graph into the data
structure, and then expands the simplicial complex when required.

Vertex name correspond to the index of the point in the given range (aka. the point cloud).

.. figure::
    img/rips_complex_representation.png
    :align: center

    Rips-complex one skeleton graph representation

On this example, as edges (4,5), (4,6) and (5,6) are in the complex, simplex (4,5,6) is added with the filtration value
set with :math:`max(filtration(4,5), filtration(4,6), filtration(5,6))`. And so on for simplex (0,1,2,3).

If the Rips_complex interfaces are not detailed enough for your need, please refer to rips_persistence_step_by_step.cpp
example, where the graph construction over the Simplex_tree is more detailed.

Point cloud and distance function
---------------------------------

Example from a point cloud and a distance function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example builds the one skeleton graph from the given points, threshold value, and distance function. Then it
creates a :doc:`Simplex_tree <simplex_tree_ref>` with it.

Then, it is asked to display information about the simplicial complex.

.. testcode::

    import gudhi
    rips_complex = gudhi.RipsComplex(points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]],
        max_edge_length=12.0)

    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    for filtered_value in simplex_tree.get_filtered_tree():
        print(filtered_value)

The output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    ([0], 0.0)
    ([1], 0.0)
    ([2], 0.0)
    ([3], 0.0)
    ([4], 0.0)
    ([5], 0.0)
    ([6], 0.0)
    ([2, 3], 5.0)
    ([4, 5], 5.385164807134504)
    ([0, 2], 5.830951894845301)
    ([0, 1], 6.082762530298219)
    ([1, 3], 6.324555320336759)
    ([1, 2], 6.708203932499369)
    ([5, 6], 7.280109889280518)
    ([2, 4], 8.94427190999916)
    ([0, 3], 9.433981132056603)
    ([4, 6], 9.486832980505138)
    ([3, 6], 11.0)

Example from OFF file
^^^^^^^^^^^^^^^^^^^^^

This example builds the :doc:`Rips_complex <rips_complex_ref>` from the given points in an OFF file, threshold value,
and distance function. Then it creates a :doc:`Simplex_tree <simplex_tree_ref>` with it.

Then, it is asked to display information about the Rips complex.


.. testcode::

    import gudhi
    rips_complex = gudhi.RipsComplex(off_file='alphacomplexdoc.off', max_edge_length=12.0)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        repr(simplex_tree.num_vertices()) + ' vertices.'
    print(result_str)
    for filtered_value in simplex_tree.get_filtered_tree():
        print(filtered_value)

the program output is:

.. testoutput::

    Rips complex is of dimension 1 - 18 simplices - 7 vertices.
    ([0], 0.0)
    ([1], 0.0)
    ([2], 0.0)
    ([3], 0.0)
    ([4], 0.0)
    ([5], 0.0)
    ([6], 0.0)
    ([2, 3], 5.0)
    ([4, 5], 5.385164807134504)
    ([0, 2], 5.830951894845301)
    ([0, 1], 6.082762530298219)
    ([1, 3], 6.324555320336759)
    ([1, 2], 6.708203932499369)
    ([5, 6], 7.280109889280518)
    ([2, 4], 8.94427190999916)
    ([0, 3], 9.433981132056603)
    ([4, 6], 9.486832980505138)
    ([3, 6], 11.0)
