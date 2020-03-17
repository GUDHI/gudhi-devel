.. table::
   :widths: 30 40 30

   +-----------------------------------------------------------------+-----------------------------------------------------------------------+-----------------------------------------------+
   | .. figure::                                                     | The theory of homology consists in attaching to a topological space   | :Author: Clément Maria                        |
   |      ../../doc/Persistent_cohomology/3DTorus_poch.png           | a sequence of (homology) groups, capturing global topological         |                                               |
   |      :figclass: align-center                                    | features like connected components, holes, cavities, etc. Persistent  | :Introduced in: GUDHI 2.0.0                   |
   |                                                                 | homology studies the evolution -- birth, life and death -- of these   |                                               |
   |      Rips Persistent Cohomology on a 3D                         | features when the topological space is changing. Consequently, the    | :Copyright: MIT                               |
   |      Torus                                                      | theory is essentially composed of three elements: topological spaces, |                                               |
   |                                                                 | their homology groups and an evolution scheme.                        |                                               |
   |                                                                 |                                                                       |                                               |
   |                                                                 | Computation of persistent cohomology using the algorithm of           |                                               |
   |                                                                 | :cite:`DBLP:journals/dcg/SilvaMV11` and                               |                                               |
   |                                                                 | :cite:`DBLP:journals/corr/abs-1208-5018` and the Compressed           |                                               |
   |                                                                 | Annotation Matrix implementation of                                   |                                               |
   |                                                                 | :cite:`DBLP:conf/esa/BoissonnatDM13`.                                 |                                               |
   |                                                                 |                                                                       |                                               |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------+-----------------------------------------------+
   |  * :doc:`persistent_cohomology_user`                            | Please refer to each data structure that contains persistence                                                         |
   |                                                                 | feature for reference:                                                                                                |
   |                                                                 |                                                                                                                       |
   |                                                                 | * :doc:`simplex_tree_ref`                                                                                             |
   |                                                                 | * :doc:`cubical_complex_ref`                                                                                          |
   |                                                                 | * :doc:`periodic_cubical_complex_ref`                                                                                 |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
