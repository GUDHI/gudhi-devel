.. GUDHI documentation master file, created by
   sphinx-quickstart on Thu Jun 30 09:55:51 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

GUDHI's documentation
#####################

.. image:: img/Gudhi_banner.png
   :align: center

Introduction
************

The Gudhi library (Geometry Understanding in Higher Dimensions) is a generic
open source C++ library for Computational Topology and Topological Data
Analysis (TDA). The GUDHI library intends to help the development of new
algorithmic solutions in TDA and their transfer to applications. It provides
robust, efficient, flexible and easy to use implementations of
state-of-the-art algorithms and data structures.

The current release of the GUDHI library includes:

* Data structures to represent, construct and manipulate simplicial complexes.
* Algorithms to compute persistent homology and multi-field persistent homology.
* Simplication of simplicial complexes by edge contraction.

All data-structures are generic and several of their aspects can be
parameterized via template classes. We refer to [14] for a detailed
description of the design of the library.

Data structures
***************

Alpha complex
=============

.. image:: img/alpha_complex_representation.png
   :align: center

Alpha_complex is a simplicial complex constructed from the finite cells of a
Delaunay Triangulation.

The filtration value of each simplex is computed as the square of the
circumradius of the simplex if the circumsphere is empty (the simplex is then
said to be Gabriel), and as the minimum of the filtration values of the
codimension 1 cofaces that make it not Gabriel otherwise. All simplices that
have a filtration value strictly greater than a given alpha squared value are
not inserted into the complex.

This package requires having CGAL version 4.7 or higher (4.8.1 is advised for
better perfomances).

* :ref:`alphacomplex`


.. automodule:: gudhi
   :members:
   :undoc-members:
   :show-inheritance:

