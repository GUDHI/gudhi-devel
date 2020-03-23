.. table::
   :widths: 30 40 30

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | The q-Wasserstein distance measures the similarity between two       | :Author: Theo Lacombe                                            |
   |      ../../doc/Bottleneck_distance/perturb_pd.png               | persistence diagrams. It's the minimum value c that can be achieved  |                                                                  |
   |      :figclass: align-center                                    | by a perfect matching between the points of the two diagrams (+ all  | :Since: GUDHI 3.1.0                                              |
   |                                                                 | diagonal points), where the value of a matching is defined as the    |                                                                  |
   |      Wasserstein distance is the q-th root of the sum of the    | q-th root of the sum of all edge lengths to the power q. Edge lengths| :License: MIT                                                    |
   |      edge lengths to the power q.                               | are measured in norm p, for :math:`1 \leq p \leq \infty`.            |                                                                  |
   |                                                                 |                                                                      | :Requires: Python Optimal Transport (POT) :math:`\geq` 0.5.1     |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`wasserstein_distance_user`                              |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
