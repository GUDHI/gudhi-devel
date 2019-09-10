:orphan:

.. To get rid of WARNING: document isn't included in any toctree

Cubical complex user manual
===========================
Definition
----------

=====================================  =====================================  =====================================
:Author: Pawel Dlotko                  :Introduced in: GUDHI PYTHON 2.0.0     :Copyright: GPL v3
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

.. figure::
    ../../doc/Bitmap_cubical_complex/Cubical_complex_representation.png
    :alt: Cubical complex.
    :figclass: align-center

    Cubical complex.

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
Currently one input from a text file is used. It uses a format inspired from the Perseus software
`Perseus software <http://www.sas.upenn.edu/~vnanda/perseus/>`_ by Vidit Nanda.

.. note::
    While Perseus assume the filtration of all maximal cubes to be non-negative, over here we do not enforce this and
    we allow any filtration values. As a consequence one cannot use ``-1``'s to indicate missing cubes. If you have
    missing cubes in your complex, please set their filtration to :math:`+\infty` (aka. ``inf`` in the file).

The file format is described in details in :ref:`Perseus file format` file format section.

.. testcode::

    import gudhi
    cubical_complex = gudhi.CubicalComplex(perseus_file=gudhi.__root_source_dir__ + \
        '/data/bitmap/cubicalcomplexdoc.txt')
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
complex with periodic boundary conditions.

One can also use Perseus style input files (see :doc:`Perseus <fileformats>`) for the specific periodic case:

.. testcode::

    import gudhi
    periodic_cc = gudhi.PeriodicCubicalComplex(perseus_file=gudhi.__root_source_dir__ + \
        '/data/bitmap/periodiccubicalcomplexdoc.txt')
    result_str = 'Periodic cubical complex is of dimension ' + repr(periodic_cc.dimension()) + ' - ' + \
        repr(periodic_cc.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::
    
    Periodic cubical complex is of dimension 2 - 42 simplices.

Or it can be defined as follows:

.. testcode::

    from gudhi import PeriodicCubicalComplex as pcc
    periodic_cc = pcc(dimensions=[3,3],
         top_dimensional_cells= [0, 0, 0, 0, 1, 0, 0, 0, 0],
         periodic_dimensions=[True, False])
    result_str = 'Periodic cubical complex is of dimension ' + repr(periodic_cc.dimension()) + ' - ' + \
        repr(periodic_cc.num_simplices()) + ' simplices.'
    print(result_str)

the program output is:

.. testoutput::

    Periodic cubical complex is of dimension 2 - 42 simplices.

Examples.
---------

End user programs are available in python/example/ folder.

Bibliography
============

.. bibliography:: ../../biblio/bibliography.bib
   :filter: docnames
   :style: unsrt
