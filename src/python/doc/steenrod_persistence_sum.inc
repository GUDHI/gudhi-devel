.. table::
   :widths: 30 40 30

   +----------------------------------------------------------+-------------------------------------------------------------------------+---------------------------------------------------+
   | .. figure::                                              | Persistence Steenrod barcodes extend the ordinary persistence barcode   | :Author: Umberto Lupo, Anibal M.                  |
   |      ../../doc/common/Gudhi_banner.png                   | using Steenrod cohomology operations Sq\ :sup:`k` over                 | Medina-Mardones, Guillaume Tauzin                 |
   |      :figclass: align-center                             | :math:`\mathbb{F}_2`. Two filtered complexes with identical ordinary    | (port to Gudhi: Maximiliano Alvarez)              |
   |                                                          | barcodes may be distinguished by their Steenrod barcodes. The          |                                                   |
   |                                                          | computation uses a persistence matrix reduction that also stores        | :Since: GUDHI 3.13.0                              |
   |                                                          | cocycle representatives, then applies the explicit cup-i product        |                                                   |
   |                                                          | formulas of :cite:`medina2023cupifast`.                                 | :License: MIT                                     |
   |                                                          |                                                                         |                                                   |
   |                                                          | The algorithm is described in :cite:`lupo2022persistencesteenrod`.      |                                                   |
   |                                                          |                                                                         |                                                   |
   +----------------------------------------------------------+-------------------------------------------------------------------------+---------------------------------------------------+
   |                                                          | * :doc:`simplex_tree_ref` — ``SimplexTree.compute_steenrod_barcodes``  |                                                   |
   +----------------------------------------------------------+-------------------------------------------------------------------------+---------------------------------------------------+
