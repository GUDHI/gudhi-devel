.. GUDHI documentation master file, created by
   sphinx-quickstart on Thu Jun 30 09:55:51 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GUDHI's documentation
#####################

.. image:: img/Gudhi_banner.png
   :align: center

Introduction
************

The Gudhi library (Geometry Understanding in Higher Dimensions) is a generic
open source C++ library for Computational Topology and Topological Data
Analysis (`TDA <https://en.wikipedia.org/wiki/Topological_data_analysis>`_).
The GUDHI library intends to help the development of new algorithmic solutions
in TDA and their transfer to applications. It provides robust, efficient,
flexible and easy to use implementations of state-of-the-art algorithms and
data structures.

The current release of the GUDHI library includes:

* Data structures to represent, construct and manipulate simplicial complexes.
* Algorithms to compute persistent homology and multi-field persistent homology.
* Simplication of simplicial complexes by edge contraction.

All data-structures are generic and several of their aspects can be
parameterized via template classes. We refer to :cite:`gudhilibrary_ICMS14`
for a detailed description of the design of the library.

Data structures
***************

Alpha complex
=============

.. include:: alpha_complex_sum.rst

Cubical complex
===============

.. include:: cubical_complex_sum.rst

Simplex tree
============

.. include:: simplex_tree_sum.rst

Witness complex
===============

.. include:: witness_complex_sum.rst


Toolbox
*******

Persistence cohomology
======================

.. include:: persistent_cohomology_sum.rst

Persistence graphical tools
===========================

.. include:: persistence_graphical_tools_sum.rst
