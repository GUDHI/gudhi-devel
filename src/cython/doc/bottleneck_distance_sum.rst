=================================================================  ===================================  ===================================
:Author: François Godi                                             :Introduced in: GUDHI 2.0.0          :Copyright: GPL v3
:Requires: CGAL :math:`\geq` 4.8.0
=================================================================  ===================================  ===================================

+-----------------------------------------------------------------+----------------------------------------------------------------------+
| .. figure::                                                     | Bottleneck distance measures the similarity between two persistence  |
|      img/perturb_pd.png                                         | diagrams. It's the shortest distance b for which there exists a      |
|      :figclass: align-center                                    | perfect matching between the points of the two diagrams (+ all the   |
|                                                                 | diagonal points) such that any couple of matched points are at       |
|      Bottleneck distance is the length of                       | distance at most b.                                                  |
|      the longest edge                                           |                                                                      |
+-----------------------------------------------------------------+----------------------------------------------------------------------+
| :doc:`bottleneck_distance_user`                                 | :doc:`bottleneck_distance_ref`                                       |
+-----------------------------------------------------------------+----------------------------------------------------------------------+
