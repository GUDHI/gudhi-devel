.. table::
   :widths: 30 40 30

   +----------------------------------------------------------------+------------------------------------------------------------------------+------------------------------------------------------------------------------------+
   | .. figure::                                                    | A Tangential Delaunay complex is a simplicial complex designed to      | :Author:    Clément Jamin                                                          |
   |      ../../doc/Tangential_complex/tc_examples.png              | reconstruct a :math:`k`-dimensional manifold embedded in               |                                                                                    |
   |      :figclass: align-center                                   | :math:`d`-dimensional Euclidean space. The input is a point sample     | :Since:    GUDHI 2.0.0                                                             |
   |                                                                | coming from an unknown manifold. The running time depends only linearly|                                                                                    |
   |                                                                | on the extrinsic dimension :math:`d` and exponentially on the intrinsic| :License: MIT (`GPL v3 </licensing/>`_)                                            |
   |                                                                | dimension :math:`k`.                                                   |                                                                                    |
   |                                                                |                                                                        | :Requires: `Eigen <installation.html#eigen>`_ and `CGAL <installation.html#cgal>`_ |
   +----------------------------------------------------------------+------------------------------------------------------------------------+------------------------------------------------------------------------------------+
   | * :doc:`tangential_complex_user`                               | * :doc:`tangential_complex_ref`                                                                                                                             |
   +----------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------+
