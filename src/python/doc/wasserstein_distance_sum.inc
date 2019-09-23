.. table::
   :widths: 30 50 20

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | The p-Wasserstein distance measures the similarity between two       | :Author: Theo Lacombe                                            |
   |      ../../doc/Bottleneck_distance/perturb_pd.png               | persistence diagrams. It's the minimum value c that can be achieve by|                                                                  |
   |      :figclass: align-center                                    | a perfect matching between the points of the two diagrams (+ all the | :Introduced in: GUDHI 2.0.0                                      |
   |                                                                 | diagonal points), where the value of a matching is defined as the    |                                                                  |
   |      Wasserstein distance is the p-th root of the sum of the    | p-th root of the sum of all edges lengths to the power p. Edges      | :Copyright: MIT (`GPL v3 </licensing/>`_)                        |
   |      edges lengths to the power p.                              | lengths are measured in norm q, for $1 \leq q \leq \infty$.          |                                                                  |
   |                                                                 |                                                                      | :Requires: `Python Optimal Transport (POT)`                      |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`wasserstein_distance_user`                              |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
