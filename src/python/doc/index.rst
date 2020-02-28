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
   :widths: 30 15 15 20 20

   +----------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | .. figure::                                        |                                                                                                                                                                       |
   |     ../../doc/Bottleneck_distance/perturb_pd.png   +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |     :figclass: align-center                        | :Author: François Godi | :Introduced in: GUDHI 2.0.0 | :Copyright: MIT (`GPL v3 </licensing/>`_) | :Requires: `CGAL <installation.html#cgal>`__ :math:`\geq` 4.11.0   |
   |                                                    +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |                                                    | Bottleneck distance measures the similarity between two persistence diagrams. It's the shortest distance b for which there exists a perfect matching between the      |
   |                                                    | points of the two diagrams (+ all the diagonal points) such that any couple of matched points are at distance at most b, where the distance between points is the     |
   |                                                    | sup norm in :math:`\mathbb{R}^2`.                                                                                                                                     |
   |                                                    |                                                                                                                                                                       |
   |                                                    | * :doc:`bottleneck_distance_user`                                                                                                                                     |
   |                                                    +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |                                                    | :Author: Theo Lacombe  | :Introduced in: GUDHI 3.1.0 | :Copyright: MIT                           | :Requires: Python Optimal Transport (POT) :math:`\geq` 0.5.1       |
   |                                                    +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |                                                    | The q-Wasserstein distance measures the similarity between two persistence diagrams. It's the minimum value c that can be achieved by a perfect matching between the  |
   |                                                    | points of the two diagrams (+ all diagonal points), where the value of a matching is defined as the q-th root of the sum of all edge lengths to the power q. Edge     |
   |                                                    | lengths are measured in norm p, for :math:`1 \leq p \leq \infty`.                                                                                                     |
   |                                                    |                                                                                                                                                                       |
   |                                                    | * :doc:`wasserstein_distance_user`                                                                                                                                    |
   |                                                    +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |                                                    | :Author: Marc Glisse   | :Introduced in: GUDHI 3.2.0 | :Copyright: MIT                           | :Uses: `Hera <https://bitbucket.org/grey_narn/hera/src/master/>`__ |
   |                                                    +------------------------+-----------------------------+-------------------------------------------+--------------------------------------------------------------------+
   |                                                    | This other q-Wasserstein distance comes from `Hera <https://bitbucket.org/grey_narn/hera/src/master/>`_ (BSD-3-Clause) which is based on "Geometry Helps to Compare   |
   |                                                    | Persistence Diagrams" :cite:`Kerber:2017:GHC:3047249.3064175` by Michael Kerber, Dmitriy Morozov, and Arnur Nigmetov.                                                 |
   |                                                    |                                                                                                                                                                       |
   |                                                    | * :doc:`wasserstein_distance_user`                                                                                                                                    |
   +----------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+


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
