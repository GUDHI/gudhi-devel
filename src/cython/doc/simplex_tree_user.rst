========================
Simplex tree user manual
========================
Definition
----------

.. include:: simplex_tree_sum.rst

A simplicial complex :math:`\mathbf{K}` on a set of vertices :math:`V = \{1, \cdots ,|V|\}` is a collection of
simplices :math:`\{\sigma\}`, :math:`\sigma \subseteq V` such that
:math:`\tau \subseteq \sigma \in \mathbf{K} \rightarrow \tau \in \mathbf{K}`. The dimension :math:`n=|\sigma|-1` of
:math:`\sigma` is its number of elements minus `1`.

A filtration of a simplicial complex is a function :math:`f:\mathbf{K} \rightarrow \mathbb{R}` satisfying
:math:`f(\tau)\leq f(\sigma)` whenever :math:`\tau \subseteq \sigma`. Ordering the simplices by increasing filtration
values (breaking ties so as a simplex appears after its subsimplices of same filtration value) provides an indexing
scheme.


Implementation
--------------

There are two implementation of complexes. The first on is the Simplex_tree data structure.
The simplex tree is an efficient and flexible data structure for representing general (filtered) simplicial complexes.
The data structure is described in :cite`boissonnatmariasimplextreealgorithmica`.

The second one is the Hasse_complex. The Hasse complex is a data structure representing explicitly all co-dimension 1
incidence relations in a complex. It is consequently faster when accessing the boundary of a simplex, but is less
compact and harder to construct from scratch.

Example
-------

.. testcode::

   import gudhi
   st = gudhi.SimplexTree()
   if st.insert_simplex([0, 1]):
       print("[0, 1] inserted")
   if st.insert_simplex_and_subfaces([0, 1, 2], filtration=4.0):
       print("[0, 1, 2] inserted")
   if st.find([0, 1]):
       print("[0, 1] found")
   result_str = 'num_vertices=' + repr(st.num_vertices())
   print(result_str)
   result_str = 'num_simplices=' + repr(st.num_simplices())
   print(result_str)
   print("skeleton_tree(2) =")
   for sk_value in st.get_skeleton_tree(2):
       print(sk_value)


The output is:

.. testoutput::

   [0, 1] inserted
   [0, 1, 2] inserted
   [0, 1] found
   num_vertices=3
   num_simplices=7
   skeleton_tree(2) =
   ([0, 1, 2], 4.0)
   ([0, 1], 0.0)
   ([0, 2], 4.0)
   ([0], 0.0)
   ([1, 2], 4.0)
   ([1], 4.0)
   ([2], 4.0)
