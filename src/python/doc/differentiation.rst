:orphan:

.. To get rid of WARNING: document isn't included in any toctree

======================
Differentiation manual
======================

.. include:: differentiation_sum.inc

In this module, we provide neural network models for computing persistent homology. In particular, we provide TensorFlow 2 models that allow to compute persistence diagrams from complexes available in the Gudhi library, including simplex trees, cubical complexes and Vietoris-Rips complexes. These models can be incorporated at each step of a given neural network architecture, and can be used in addition to `PersLay <https://github.com/MathieuCarriere/gudhi/blob/perslay/src/python/gudhi/representations/perslay.py>`_ to produce topological features.

TensorFlow models
-----------------
.. automodule:: gudhi.differentiation
   :members:
   :special-members:
   :show-inheritance:
