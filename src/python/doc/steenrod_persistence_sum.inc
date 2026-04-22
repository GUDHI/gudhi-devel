.. table::
   :widths: 30 40 30

   +-----------------------------------------------------------------+-----------------------------------------------------------------------+-----------------------------------------------+
   | .. figure::                                                     | Persistence Steenrod barcodes extend the ordinary persistence         | :Author: Lupo, Medina-Mardones, Tauzin        |
   |      ../../doc/Steenrod_persistence/steenrod_barcode.png        | barcode using Steenrod cohomology operations Sq\ :sup:`k`             | :Since: GUDHI 3.13.0                          |
   |      :figclass: align-center                                    | over :math:`\mathbb{F}_2`. Two filtered complexes with                | :License: MIT                                 |
   |                                                                 | identical ordinary barcodes may be distinguished by their             |                                               |
   |                                                                 | Steenrod barcodes.                                                    |                                               |
   |                                                                 |                                                                       |                                               |
   |                                                                 | The computation reduces a persistence matrix while storing            |                                               |
   |                                                                 | cocycle representatives, then applies the cup-i product               |                                               |
   |                                                                 | formulas of :cite:`medina2023cupifast`.                               |                                               |
   |                                                                 |                                                                       |                                               |
   |                                                                 | The algorithm is described in                                         |                                               |
   |                                                                 | :cite:`lupo2022persistencesteenrod`.                                  |                                               |
   |                                                                 |                                                                       |                                               |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------+-----------------------------------------------+
   | * :doc:`steenrod_persistence_user`                              |  * :doc:`simplex_tree_ref` — ``SimplexTree.compute_steenrod_barcodes``                                                |
   +-----------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------+
