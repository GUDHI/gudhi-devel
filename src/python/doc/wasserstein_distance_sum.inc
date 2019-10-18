.. table::
   :widths: 30 50 20

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | The p-Wasserstein distance measures the similarity between two       | :Author: Theo Lacombe                                            |
   |      ../../doc/Bottleneck_distance/perturb_pd.png               | persistence diagrams. It's the minimum value c that can be achieved  |                                                                  |
   |      :figclass: align-center                                    | by a perfect matching between the points of the two diagrams (+ all  | :Introduced in: GUDHI 3.1.0                                      |
   |                                                                 | diagonal points), where the value of a matching is defined as the    |                                                                  |
   |      Wasserstein distance is the p-th root of the sum of the    | p-th root of the sum of all edge lengths to the power p. Edge lengths| :Copyright: MIT                                                  |
   |      edge lengths to the power p.                               | are measured in norm q, for :math:`1 \leq q \leq \infty`.            |                                                                  |
   |                                                                 |                                                                      | :Requires: Python Optimal Transport (POT) :math:`\geq` 0.5.1     |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`wasserstein_distance_user`                              |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
