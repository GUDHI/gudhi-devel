=================================================================  ===================================  ===================================
:Author: Siargey Kachanovich                                       :Introduced in: GUDHI 2.0.0          :Copyright: GPL v3
=================================================================  ===================================  ===================================

+-----------------------------------------------------------------+----------------------------------------------------------------------+
| .. image::                                                      | Witness complex :math:`Wit(W,L)` is a simplicial complex defined on  |
|      img/Witness_complex_representation.png                     | two sets of points in :math:`\mathbb{R}^D`:Wit(W,L)` is a simplicial |
|                                                                 | complex defined on two sets of points in :math:`\mathbb{R}^D`:       |
|                                                                 |                                                                      |
|                                                                 | * :math:`W` set of **witnesses** and                                 |
|                                                                 | * :math:`L \subseteq W` set of **landmarks**.                        |
|                                                                 |                                                                      |
|                                                                 | The simplices are based on landmarks and a simplex belongs to the    |
|                                                                 | witness complex if and only if it is witnessed, that is:             |
|                                                                 |                                                                      |
|                                                                 | :math:`\sigma \subset L` is witnessed if there exists a point        |
|                                                                 | :math:`w \in W` such that w is closer to the vertices of             |
|                                                                 | :math:`\sigma` than other points in :math:`L` and all of its faces   |
|                                                                 | are witnessed as well.                                               |
|                                                                 |                                                                      |
|                                                                 | The data structure is described in                                   |
|                                                                 | :cite:`boissonnatmariasimplextreealgorithmica`.                      |
+-----------------------------------------------------------------+----------------------------------------------------------------------+
| :doc:`witness_complex_user`                                     | :doc:`witness_complex_ref`                                           |
+-----------------------------------------------------------------+----------------------------------------------------------------------+
