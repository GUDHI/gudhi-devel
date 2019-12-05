.. table::
   :widths: 30 50 20

   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | .. figure::                                                     | A Frechet mean (or barycenter) is a generalization of the arithmetic | :Author: Theo Lacombe                                            |
   |      ./img/barycenter.png                                       | mean in a non linear space such as the one of persistence diagrams.  |                                                                  |
   |      :figclass: align-center                                    | Given a set of persistence diagrams :math:`\mu_1 \dots \mu_n`, it is | :Introduced in: GUDHI 3.1.0                                      |
   |                                                                 | defined as a minimizer of the variance functional, that is of        |                                                                  |
   |      Illustration of Frechet mean between persistence           | :math:`\mu \mapsto \sum_{i=1}^n d_2(\mu,\mu_i)^2`.                   | :Copyright: MIT                                                  |
   |      diagrams.                                                  | where :math:`d_2` denotes the Wasserstein-2 distance between persis- |                                                                  |
   |                                                                 | tence diagrams.                                                      |                                                                  |
   |                                                                 | It is known to exist and is generically unique. However, an exact    |                                                                  |
   |                                                                 | computation is in general untractable. Current implementation avai-  | :Requires: Python Optimal Transport (POT) :math:`\geq` 0.5.1     |
   |                                                                 | -lable is based on [Turner et al, 2014], and uses an EM-scheme to    |                                                                  |
   |                                                                 | provide a local minimum of the variance functional (somewhat similar |                                                                  |
   |                                                                 | to the Lloyd algorithm to estimate a solution to the k-means         |                                                                  |
   |                                                                 | problem). The local minimum returned depends on the initialization of|                                                                  |
   |                                                                 | the barycenter.                                                      |                                                                  |
   |                                                                 | The combinatorial structure of the algorithm limits its              |                                                                  |
   |                                                                 | scaling on large scale problems (thousands of diagrams and of points |                                                                  |
   |                                                                 | per diagram).                                                        |                                                                  |
   +-----------------------------------------------------------------+----------------------------------------------------------------------+------------------------------------------------------------------+
   | * :doc:`barycenter_user`                                        |                                                                                                                                         |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
