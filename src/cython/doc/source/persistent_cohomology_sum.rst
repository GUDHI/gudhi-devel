=====================================  =====================================  =====================================
:Author: Clément Maria                 :Introduced in: GUDHI PYTHON 1.4.0     :Copyright: GPL v3
=====================================  =====================================  =====================================

+---------------------------------------------+----------------------------------------------------------------------+
| .. image::                                  | The theory of homology consists in attaching to a topological space  |
|      img/3DTorus_poch.png                   | a sequence of (homology) groups, capturing global topological        |
|                                             | features like connected components, holes, cavities, etc. Persistent |
|                                             | homology studies the evolution -- birth, life and death -- of these  |
|                                             | features when the topological space is changing. Consequently, the   |
|                                             | theory is essentially composed of three elements: topological spaces,|
|                                             | their homology groups and an evolution scheme.                       |
|                                             |                                                                      |
|                                             | Computation of persistent cohomology using the algorithm of          |
|                                             | :cite:`DBLP:journals/dcg/SilvaMV11` and                              |
|                                             | :cite:`DBLP:journals/corr/abs-1208-5018` and the Compressed          |
|                                             | Annotation Matrix implementation of                                  |
|                                             | :cite:`DBLP:conf/esa/BoissonnatDM13`.                                |
|                                             |                                                                      |
+---------------------------------------------+----------------------------------------------------------------------+
|  :doc:`persistent_cohomology_user`          | Please refer to each data structure that contains persistence        |
|                                             | feature for reference:                                               |
|                                             |                                                                      |
|                                             | * :doc:`alpha_complex_ref`                                           |
|                                             | * :doc:`cubical_complex_ref`                                         |
|                                             | * :doc:`simplex_tree_ref`                                            |
|                                             | * :doc:`witness_complex_ref`                                         |
+---------------------------------------------+----------------------------------------------------------------------+
