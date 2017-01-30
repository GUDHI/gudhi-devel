==============================
Tangential complex user manual
==============================
.. include:: tangential_complex_sum.rst

Definition
----------

A Tangential Delaunay complex is a simplicial complex designed to reconstruct a
:math:`k`-dimensional smooth manifold embedded in :math:`d`-dimensional
Euclidean space. The input is a point sample coming from an unknown manifold,
which means that the points lie close to a structure of "small" intrinsic
dimension. The running time depends only linearly on the extrinsic dimension
:math:`d` and exponentially on the intrinsic dimension :math:`k`.

An extensive description of the Tangential complex can be found in
:cite:`tangentialcomplex2014`).

What is a Tangential Complex?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us start with the description of the Tangential complex of a simple
example, with :math:`k = 1` and :math:`d = 2`. The input data is 4 points
:math:`P` located on a curve embedded in 2D.

.. figure:: img/tc_example_01.png
    :alt: The input
    :figclass: align-center
    The input

For each point :math:`p`, estimate its tangent subspace :math:`T_p` (e.g.
using PCA). 

.. figure:: img/tc_example_02.png
    :alt: The estimated normals
    :figclass: align-center
    The estimated normals

Let us add the Voronoi diagram of the points in orange. For each point
:math:`p`, construct its star in the Delaunay triangulation of :math:`P`
restricted to :math:`T_p`.

.. figure:: img/tc_example_03.png
    :alt: The Voronoi diagram
    :figclass: align-center
    The Voronoi diagram

The Tangential Delaunay complex is the union of those stars.

In practice, neither the ambient Voronoi diagram nor the ambient Delaunay
triangulation is computed. Instead, local :math:`k`-dimensional regular
triangulations are computed with a limited number of points as we only need the
star of each point. More details can be found in :cite:`tangentialcomplex2014`.

Inconsistencies
^^^^^^^^^^^^^^^
Inconsistencies between the stars can occur. An inconsistency occurs when a
simplex is not in the star of all its vertices.

Let us take the same example.

.. figure:: img/tc_example_07_before.png
    :alt: Before
    :figclass: align-center
    Before

Let us slightly move the tangent subspace :math:`T_q`

.. figure:: img/tc_example_07_after.png
    :alt: After
    :figclass: align-center
    After

Now, the star of :math:`Q` contains :math:`QP`, but the star of :math:`P` does
not contain :math:`QP`. We have an inconsistency.

.. figure:: img/tc_example_08.png
    :alt: After
    :figclass: align-center
    After

One way to solve inconsistencies is to randomly perturb the positions of the
points involved in an inconsistency. In the current implementation, this
perturbation is done in the tangent subspace of each point. The maximum
perturbation radius is given as a parameter to the constructor.

In most cases, we recommend to provide a point set where the minimum distance
between any two points is not too small. This can be achieved using the
functions provided by the Subsampling module. Then, a good value to start with
for the maximum perturbation radius would be around half the minimum distance
between any two points. The Example with perturbation below shows an example of
such a process.

In most cases, this process is able to dramatically reduce the number of
inconsistencies, but is not guaranteed to succeed.

Output
^^^^^^
The result of the computation is exported as a Simplex_tree. It is the union of
the stars of all the input points. A vertex in the Simplex Tree is the index of
the point in the range provided by the user. The point corresponding to a
vertex can also be obtained through the Tangential_complex::get_point function.
Note that even if the positions of the points are perturbed, their original
positions are kept (e.g. Tangential_complex::get_point returns the original
position of the point).

The result can be obtained after the computation of the Tangential complex
itself and/or after the perturbation process.


Simple example
--------------

This example builds the Tangential complex of point set read in an OFF file.

.. testcode::

    import gudhi
    tc = gudhi.TangentialComplex(off_file='alphacomplexdoc.off')
    result_str = 'Tangential contains ' + repr(tc.num_simplices()) + \
        ' simplices - ' + repr(tc.num_vertices()) + ' vertices.'
    print(result_str)

    st = tc.create_simplex_tree()
    result_str = 'Simplex tree is of dimension ' + repr(st.dimension()) + \
        ' - ' + repr(st.num_simplices()) + ' simplices - ' + \
        repr(st.num_vertices()) + ' vertices.'
    print(result_str)
    for filtered_value in st.get_filtered_tree():
        print(filtered_value)

The output is:

.. testoutput::

    Tangential contains 12 simplices - 7 vertices.
    Simplex tree is of dimension 1 - 15 simplices - 7 vertices.
    ([0], 0.0)
    ([1], 0.0)
    ([0, 1], 0.0)
    ([2], 0.0)
    ([0, 2], 0.0)
    ([1, 2], 0.0)
    ([3], 0.0)
    ([1, 3], 0.0)
    ([4], 0.0)
    ([2, 4], 0.0)
    ([5], 0.0)
    ([4, 5], 0.0)
    ([6], 0.0)
    ([3, 6], 0.0)
    ([5, 6], 0.0)


Example with perturbation
-------------------------

This example builds the Tangential complex of a point set, then tries to solve
inconsistencies by perturbing the positions of points involved in inconsistent
simplices.

.. testcode::

   import gudhi
   tc = gudhi.TangentialComplex(points=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
   result_str = 'Tangential contains ' + repr(tc.num_vertices()) + ' vertices.'
   print(result_str)

   if tc.num_inconsistent_simplices() > 0:
       print('Tangential contains inconsistencies.')

   tc.fix_inconsistencies_using_perturbation(10, 60)
   if tc.num_inconsistent_simplices() == 0:
       print('Inconsistencies has been fixed.')

The output is:

.. testoutput::

    Tangential contains 4 vertices.
    Inconsistencies has been fixed.

.. include:: biblio.rst
