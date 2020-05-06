.. table::
   :widths: 30 40 30

   +----------------------------------------------------------------+------------------------------------------------------------------------+--------------------------------------------------------------------------------+
   | .. figure::                                                    | Alpha complex is a simplicial complex constructed from the finite      | :Author: Vincent Rouvreau                                                      |
   |      ../../doc/Alpha_complex/alpha_complex_representation.png  | cells of a Delaunay Triangulation.                                     |                                                                                |
   |      :alt: Alpha complex representation                        |                                                                        | :Since: GUDHI 2.0.0                                                            |
   |      :figclass: align-center                                   | The filtration value of each simplex is computed as the **square** of  |                                                                                |
   |                                                                | the circumradius of the simplex if the circumsphere is empty (the      | :License: MIT (`GPL v3 </licensing/>`_)                                        |
   |                                                                | simplex is then said to be Gabriel), and as the minimum of the         |                                                                                |
   |                                                                | filtration values of the codimension 1 cofaces that make it not        | :Requires: :ref:`Eigen` :math:`\geq` 3.1.0 and :ref:`CGAL` :math:`\geq` 4.11.0 |
   |                                                                | Gabriel otherwise.                                                     |                                                                                |
   |                                                                |                                                                        |                                                                                |
   |                                                                | For performances reasons, it is advised to use :ref:`CGAL`             |                                                                                |
   |                                                                | :math:`\geq` 5.0.0.                                                    |                                                                                |
   +----------------------------------------------------------------+------------------------------------------------------------------------+--------------------------------------------------------------------------------+
   | * :doc:`alpha_complex_user`                                    | * :doc:`alpha_complex_ref`                                                                                                                              |
   +----------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
