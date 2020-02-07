:orphan:

.. To get rid of WARNING: document isn't included in any toctree

==============
PersLay manual
==============

.. include:: perslay_sum.inc

This module, originally available at https://github.com/MathieuCarriere/perslay, aims at bridging the gap between persistence diagrams and deep learning, by providing the implementation in tensorflow of a neural network layer for persistence diagrams, called PersLay, defined in https://arxiv.org/abs/1904.09378. More specifically, it provides code for reproducing most of the persistence diagram representations in the literature, and optimize their parameters with a neural network.

Since tensorflow requires its inputs to be fixed-size matrices, persistence diagrams are padded  to a fixed number of points, with a third boolean coordinate for each points assessing whether the point is real or just a dummy point added after padding. A diagram is thus represented as a tensorflow tensor of shape (n,3), as can be obtained from :func:`~gudhi.SimplexTree.persistence_intervals_in_dimension` and :func:`~gudhi.representations.preprocessing.Padding` for instance. Points at infinity are represented as a tensorflow tensors of shape (n,2), storing only the birth times.

A small example is provided

.. only:: builder_html

    * :download:`perslay_visu.py <../example/perslay_visu.py>`

Readme
------

.. mdinclude:: perslay.md

Function
--------
.. autofunction:: gudhi.perslay.perslay_channel

Basic example
-------------
