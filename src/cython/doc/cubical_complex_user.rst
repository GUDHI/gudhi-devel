===========================
Cubical complex user manual
===========================
Definition
----------

=====================================  =====================================  =====================================
:Author: Pawel Dlotko                  :Introduced in: GUDHI PYTHON 1.4.0     :Copyright: GPL v3
=====================================  =====================================  =====================================

+---------------------------------------------+----------------------------------------------------------------------+
| :doc:`cubical_complex_user`                 | * :doc:`cubical_complex_ref`                                         |
|                                             | * :doc:`periodic_cubical_complex_ref`                                |
+---------------------------------------------+----------------------------------------------------------------------+

The cubical complex is an example of a structured complex useful in computational mathematics (specially rigorous
numerics) and image analysis.

An *elementary interval* is an interval of a form :math:`[n,n+1]`, or :math:`[n,n]`, for :math:`n \in \mathcal{Z}`.
The first one is called *non-degenerate*, while the second one is a *degenerate* interval. A
*boundary of a elementary interval* is a chain  :math:`\partial [n,n+1] = [n+1,n+1]-[n,n]` in case of
non-degenerated elementary interval and :math:`\partial [n,n] = 0` in case of degenerate elementary interval. An
*elementary cube* :math:`C` is a product of elementary intervals, :math:`C=I_1 \times \ldots \times I_n`.
*Embedding dimension* of a cube is n, the number of elementary intervals (degenerate or not) in the product.
A *dimension of a cube* :math:`C=I_1 \times ... \times I_n` is the number of non degenerate elementary
intervals in the product. A *boundary of a cube* :math:`C=I_1 \times \ldots \times I_n` is a chain obtained
in the following way:

.. math::

    \partial C = (\partial I_1 \times \ldots \times I_n) + (I_1 \times \partial I_2 \times \ldots \times I_n) +
    \ldots + (I_1 \times I_2 \times \ldots \times \partial I_n).

A *cubical complex* :math:`\mathcal{K}` is a collection of cubes closed under operation of taking boundary
(i.e. boundary of every cube from the collection is in the collection). A cube :math:`C` in cubical complex
:math:`\mathcal{K}` is *maximal* if it is not in a boundary of any other cube in :math:`\mathcal{K}`. A
*support* of a cube :math:`C` is the set in :math:`\mathbb{R}^n` occupied by :math:`C` (:math:`n` is the embedding
dimension of :math:`C`).

Cubes may be equipped with a filtration values in which case we have filtered cubical complex. All the cubical
complexes considered in this implementation are filtered cubical complexes (although, the range of a filtration may
be a set of two elements).

For further details and theory of cubical complexes, please consult :cite:`kaczynski2004computational` as well as the
following paper :cite:`peikert2012topological`.

Data structure.
---------------

The implementation of Cubical complex provides a representation of complexes that occupy a rectangular region in
:math:`\mathbb{R}^n`. This extra assumption allows for a memory efficient way of storing cubical complexes in a form
of so called bitmaps. Let
:math:`R = [b_1,e_1] \times \ldots \times [b_n,e_n]`, for :math:`b_1,...b_n,e_1,...,e_n \in \mathbb{Z}`,
:math:`b_i \leq d_i` be the considered rectangular region and let :math:`\mathcal{K}` be a filtered
cubical complex having the rectangle :math:`R` as its support. Note that the structure of the coordinate system gives
a way a lexicographical ordering of cells of :math:`\mathcal{K}`. This ordering is a base of the presented
bitmap-based implementation. In this implementation, the whole cubical complex is stored as a vector of the values
of filtration. This, together with dimension of :math:`\mathcal{K}` and the sizes of :math:`\mathcal{K}` in all
directions, allows to determine, dimension, neighborhood, boundary and coboundary of every cube
:math:`C \in \mathcal{K}`.

.. image::
    img/Cubical_complex_representation.png
    :align: center
    :alt: Cubical complex.

Note that the cubical complex in the figure above is, in a natural way, a product of one dimensional cubical
complexes in :math:`\mathbb{R}`. The number of all cubes in each direction is equal :math:`2n+1`, where :math:`n` is
the number of maximal cubes in the considered direction. Let us consider a cube at the position :math:`k` in the
bitmap.
Knowing the sizes of the bitmap, by a series of modulo operation, we can determine which elementary intervals are
present in the product that gives the cube :math:`C`. In a similar way, we can compute boundary and the coboundary of
each cube. Further details can be found in the literature.

Input Format.
-------------

In the current implantation, filtration is given at the maximal cubes, and it is then extended by the lower star
filtration to all cubes. There are a number of constructors that can be used to construct cubical complex by users
who want to use the code directly. They can be found in the :doc:`cubical_complex_ref`.
Currently one input from a text file is used. It uses a format used already in
`Perseus software <http://www.sas.upenn.edu/~vnanda/perseus/>`_ by Vidit Nanda.
Below we are providing a description of the format. The first line contains a number d begin the dimension of the
bitmap (2 in the example below). Next d lines are the numbers of top dimensional cubes in each dimensions (3 and 3
in the example below). Next, in lexicographical order, the filtration of top dimensional cubes is given (1 4 6 8
20 4 7 6 5 in the example below).

.. image::
    img/exampleBitmap.png
    :align: center
    :alt: Example of a input data.

The input file for the following complex is:

.. literalinclude:: cubicalcomplexdoc.txt

.. centered:: cubicalcomplexdoc.txt

.. testcode::

    import gudhi
    cubical_complex = gudhi.CubicalComplex(perseus_file='cubicalcomplexdoc.txt')
    result_str = 'Cubical complex is of dimension ' + repr(cubical_complex.dimension()) + ' - ' + \
        repr(cubical_complex.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::
    
    Cubical complex is of dimension 2 - 49 simplices.

Periodic boundary conditions.
-----------------------------

Often one would like to impose periodic boundary conditions to the cubical complex (cf.
:doc:`periodic_cubical_complex_ref`).
Let :math:`I_1\times ... \times I_n` be a box that is decomposed with a cubical complex :math:`\mathcal{K}`.
Imposing periodic boundary conditions in the direction i, means that the left and the right side of a complex
:math:`\mathcal{K}` are considered the same. In particular, if for a bitmap :math:`\mathcal{K}` periodic boundary
conditions are imposed in all directions, then complex :math:`\mathcal{K}` became n-dimensional torus. One can use
various constructors from the file Bitmap_cubical_complex_periodic_boundary_conditions_base.h to construct cubical
complex with periodic boundary conditions. One can also use Perseus style input files. To indicate periodic boundary
conditions in a given direction, then number of top dimensional cells in this direction have to be multiplied by -1.
For instance:

.. literalinclude:: periodiccubicalcomplexdoc.txt

.. centered:: periodiccubicalcomplexdoc.txt

Indicate that we have imposed periodic boundary conditions in the direction x, but not in the direction y.

.. testcode::

    import gudhi
    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file='periodiccubicalcomplexdoc.txt')
    result_str = 'Periodic cubical complex is of dimension ' + repr(periodic_cc.dimension()) + ' - ' + \
        repr(periodic_cc.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::
    
    Periodic cubical complex is of dimension 2 - 42 simplices.

Examples.
---------

End user programs are available in cython/example/ folder.
