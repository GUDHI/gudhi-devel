.. table::
   :widths: 30 40 30

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | Bottleneck distance measures the similarity between two persistence  | :Author: François Godi                                           |
   |      ../../doc/Bottleneck_distance/perturb_pd.png               | diagrams. It's the shortest distance b for which there exists a      |                                                                  |
   |      :figclass: align-center                                    | perfect matching between the points of the two diagrams (+ all the   | :Introduced in: GUDHI 2.0.0                                      |
   |                                                                 | diagonal points) such that any couple of matched points are at       |                                                                  |
   |      Bottleneck distance is the length of                       | distance at most b, where the distance between points is the sup     | :Copyright: MIT (`GPL v3 </licensing/>`_)                        |
   |      the longest edge                                           | norm in :math:`\mathbb{R}^2`.                                        |                                                                  |
   |                                                                 |                                                                      | :Requires: `CGAL <installation.html#cgal>`__ :math:`\geq` 4.11.0 |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`bottleneck_distance_user`                               |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
