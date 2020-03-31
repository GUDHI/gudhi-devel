GUDHI Python modules documentation
##################################

.. figure::
     ../../doc/common/Gudhi_banner.png
     :alt: Gudhi banner
     :figclass: align-center

Data structures for cell complexes
**********************************

Cubical complexes
=================

.. include:: cubical_complex_sum.inc

Simplicial complexes
====================

Simplex tree
------------

.. include:: simplex_tree_sum.inc

Filtrations and reconstructions
*******************************

Alpha complex
=============

.. include:: alpha_complex_sum.inc

Rips complex
============

.. include:: rips_complex_sum.inc

Witness complex
===============

.. include:: witness_complex_sum.inc

Cover complexes
===============

.. include:: nerve_gic_complex_sum.inc

Tangential complex
==================

.. include:: tangential_complex_sum.inc

Topological descriptors computation
***********************************

Persistence cohomology
======================

.. include:: persistent_cohomology_sum.inc

Topological descriptors tools
*****************************

Metrics
=======

.. table::
   :widths: 30 40 30

   +---------------------------------------------------+----------------------------------------------------------------------+--------------------------------------------------------------------+
   | .. figure::                                       | Bottleneck distance measures the similarity between two persistence  | :Author: François Godi                                             |
   |      ../../doc/Bottleneck_distance/perturb_pd.png | diagrams. It's the shortest distance b for which there exists a      |                                                                    |
   |      :figclass: align-center                      | perfect matching between the points of the two diagrams (+ all the   | :Since: GUDHI 2.0.0                                                |
   |                                                   | diagonal points), where the value of a matching is defined as the    |                                                                    |
   |                                                   | q-th root of the sum of all edge lengths to the power q. Edge        | :License: MIT (`GPL v3 </licensing/>`_)                            |
   |                                                   | lengths are measured in norm p, for :math:`1 \leq p \leq \infty`.    |                                                                    |
   |                                                   |                                                                      | :Requires: `CGAL <installation.html#cgal>`__ :math:`\geq` 4.11.0   |
   |                                                   +----------------------------------------------------------------------+--------------------------------------------------------------------+
   |                                                   | * :doc:`bottleneck_distance_user`                                                                                                         |
   |                                                   +----------------------------------------------------------------------+--------------------------------------------------------------------+
   |                                                   | The q-Wasserstein distance measures the similarity between two       | :Author: Theo Lacombe                                              |
   |                                                   | persistence diagrams. It's the minimum value c that can be achieved  |                                                                    |
   |                                                   | by a perfect matching between the points of the two diagrams (+ all  | :Since: GUDHI 3.1.0                                                |
   |                                                   | diagonal points), where the value of a matching is defined as the    |                                                                    |
   |                                                   | q-th root of the sum of all edge lengths to the power q. Edge        | :License: MIT                                                      |
   |                                                   | lengths are measured in norm p, for :math:`1 \leq p \leq \infty`.    |                                                                    |
   |                                                   |                                                                      | :Requires: Python Optimal Transport (POT) :math:`\geq` 0.5.1       |
   |                                                   +----------------------------------------------------------------------+--------------------------------------------------------------------+
   |                                                   | * :doc:`wasserstein_distance_user`                                                                                                        |
   |                                                   +----------------------------------------------------------------------+--------------------------------------------------------------------+
   |                                                   | This other q-Wasserstein distance comes from                         | :Author: Marc Glisse                                               |
   |                                                   | `Hera <https://bitbucket.org/grey_narn/hera/src/master/>`_           |                                                                    |
   |                                                   | (BSD-3-Clause) which is based on "Geometry Helps to Compare          | :Since: GUDHI 3.2.0                                                |
   |                                                   | Persistence Diagrams" :cite:`Kerber:2017:GHC:3047249.3064175`        |                                                                    |
   |                                                   | by Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov.              | :License: MIT (`GPL v3 </licensing/>`_)                            |
   |                                                   |                                                                      |                                                                    |
   |                                                   |                                                                      | :Uses: `Hera <https://bitbucket.org/grey_narn/hera/src/master/>`__ |
   |                                                   +----------------------------------------------------------------------+--------------------------------------------------------------------+
   |                                                   | * :doc:`wasserstein_distance_user`                                                                                                        |
   +---------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------+


Persistence representations
===========================

.. include:: representations_sum.inc

Persistence graphical tools
===========================

.. include:: persistence_graphical_tools_sum.inc

Point cloud utilities
*********************

.. include:: point_cloud_sum.inc

Bibliography
************

.. bibliography:: ../../biblio/bibliography.bib
   :filter: docnames
   :style: unsrt
