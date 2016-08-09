=========================
Alpha complex user manual
=========================
Definition
----------

.. include:: alpha_complex_sum.rst

Alpha_complex is constructing a :doc:`Simplex_tree <simplex_tree_sum>` using
`Delaunay Triangulation  <http://doc.cgal.org/latest/Triangulation/index.html#Chapter_Triangulations>`_ 
:cite:`cgal:hdj-t-15b` from `CGAL <http://www.cgal.org/>`_ (the Computational Geometry Algorithms Library
:cite:`cgal:eb-15b`).

Remarks
^^^^^^^
When Alpha_complex is constructed with an infinite value of :math:`\alpha`, the complex is a Delaunay complex.

Example from points
-------------------

This example builds the Delaunay triangulation from the given points, and initializes the alpha complex with it:

.. testcode::

   import gudhi
   alpha_complex = gudhi.AlphaComplex(points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]],
                                      max_alpha_square=60.0)
   result_str = 'Alpha complex is of dimension ' + repr(alpha_complex.dimension()) + ' - ' + \
       repr(alpha_complex.num_simplices()) + ' simplices - ' + \
       repr(alpha_complex.num_vertices()) + ' vertices.'
   print(result_str)
   for fitered_value in alpha_complex.get_filtered_tree():
       print(fitered_value)

The output is:

.. testoutput::

   Alpha complex is of dimension 2 - 25 simplices - 7 vertices.
   ([0], 0.0)
   ([1], 0.0)
   ([2], 0.0)
   ([3], 0.0)
   ([4], 0.0)
   ([5], 0.0)
   ([6], 0.0)
   ([2, 3], 6.25)
   ([4, 5], 7.25)
   ([0, 2], 8.5)
   ([0, 1], 9.25)
   ([1, 3], 10.0)
   ([1, 2], 11.25)
   ([1, 2, 3], 12.5)
   ([0, 1, 2], 12.995867768595042)
   ([5, 6], 13.25)
   ([2, 4], 20.0)
   ([4, 6], 22.736686390532547)
   ([4, 5, 6], 22.736686390532547)
   ([3, 6], 30.25)
   ([2, 6], 36.5)
   ([2, 3, 6], 36.5)
   ([2, 4, 6], 37.24489795918368)
   ([0, 4], 59.710743801652896)
   ([0, 2, 4], 59.710743801652896)


Algorithm
---------

Data structure
^^^^^^^^^^^^^^

In order to build the alpha complex, first, a Simplex tree is built from the cells of a Delaunay Triangulation.
(The filtration value is set to NaN, which stands for unknown value):

.. image::
    img/alpha_complex_doc.png
    :align: center
    :alt: Simplex tree structure construction example

Filtration value computation algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
    \begin{algorithm}
    \caption{Filtration value computation algorithm}\label{alpha}
    \begin{algorithmic}
    \For{i : dimension $\rightarrow$ 0}
    \ForAll{$\sigma$ of dimension i}
    \If {filtration($\sigma$) is NaN}
    \State filtration($\sigma$) = $\alpha^2(\sigma)$
    \EndIf
    \ForAll{$\tau$ face of $\sigma$} \Comment{propagate alpha filtration value}
    \If {filtration($\tau$) is not NaN} 
    \State filtration($\tau$) = min (filtration($\tau$), filtration($\sigma$))
    \Else
    \If {$\tau$ is not Gabriel for $\sigma$} 
    \State filtration($\tau$) = filtration($\sigma$)
    \EndIf
    \EndIf
    \EndFor
    \EndFor
    \EndFor
    \State make\_filtration\_non\_decreasing()
    \State prune\_above\_filtration()
    \end{algorithmic}
    \end{algorithm}

Dimension 2
^^^^^^^^^^^

From the example above, it means the algorithm looks into each triangle ([0,1,2], [0,2,4], [1,2,3], ...),
computes the filtration value of the triangle, and then propagates the filtration value as described
here:

.. image::
    img/alpha_complex_doc_420.png
    :align: center
    :alt: Filtration value propagation example

Dimension 1
^^^^^^^^^^^

Then, the algorithm looks into each edge ([0,1], [0,2], [1,2], ...),
computes the filtration value of the edge (in this case, propagation will have no effect).

Dimension 0
^^^^^^^^^^^

Finally, the algorithm looks into each vertex ([0], [1], [2], [3], [4], [5] and [6]) and
sets the filtration value (0 in case of a vertex - propagation will have no effect).

Non decreasing filtration values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As the squared radii computed by CGAL are an approximation, it might happen that these alpha squared values do not
quite define a proper filtration (i.e. non-decreasing with respect to inclusion).
We fix that up by calling `Simplex_tree::make_filtration_non_decreasing()` (cf.
`C++ version <http://gudhi.gforge.inria.fr/doc/latest/index.html>`_).

Prune above given filtration value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplex tree is pruned from the given maximum alpha squared value (cf. `Simplex_tree::prune_above_filtration()`
int he `C++ version <http://gudhi.gforge.inria.fr/doc/latest/index.html>`_).
In the following example, the value is given by the user as argument of the program.


Example from OFF file
^^^^^^^^^^^^^^^^^^^^^

This example builds the Delaunay triangulation from the points given by an OFF file, and initializes the alpha complex
with it.


Then, it is asked to display information about the alpha complex:

.. testcode::

   import gudhi
   alpha_complex = gudhi.AlphaComplex(points=[[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]],
                                      max_alpha_square=59.0)
   result_str = 'Alpha complex is of dimension ' + repr(alpha_complex.dimension()) + ' - ' + \
       repr(alpha_complex.num_simplices()) + ' simplices - ' + \
       repr(alpha_complex.num_vertices()) + ' vertices.'
   print(result_str)
   for fitered_value in alpha_complex.get_filtered_tree():
       print(fitered_value)

the program output is:

.. testoutput::

   Alpha complex is of dimension 2 - 23 simplices - 7 vertices.
   ([0], 0.0)
   ([1], 0.0)
   ([2], 0.0)
   ([3], 0.0)
   ([4], 0.0)
   ([5], 0.0)
   ([6], 0.0)
   ([2, 3], 6.25)
   ([4, 5], 7.25)
   ([0, 2], 8.5)
   ([0, 1], 9.25)
   ([1, 3], 10.0)
   ([1, 2], 11.25)
   ([1, 2, 3], 12.5)
   ([0, 1, 2], 12.995867768595042)
   ([5, 6], 13.25)
   ([2, 4], 20.0)
   ([4, 6], 22.736686390532547)
   ([4, 5, 6], 22.736686390532547)
   ([3, 6], 30.25)
   ([2, 6], 36.5)
   ([2, 3, 6], 36.5)
   ([2, 4, 6], 37.24489795918368)
