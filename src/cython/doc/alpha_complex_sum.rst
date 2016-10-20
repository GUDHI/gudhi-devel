=====================================  =====================================  =====================================
:Author: Vincent Rouvreau              :Introduced in: GUDHI 1.3.0            :Copyright: GPL v3
=====================================  =====================================  =====================================
:Requires:                             CGAL &ge 4.7.0                         Eigen3
=====================================  =====================================  =====================================

+-------------------------------------------+----------------------------------------------------------------------+
| .. image::                                | Alpha_complex is a simplicial complex constructed from the finite    |
|      img/alpha_complex_representation.png | cells of a Delaunay Triangulation.                                   |
|                                           |                                                                      |
|                                           | The filtration value of each simplex is computed as the square of the|
|                                           | circumradius of the simplex if the circumsphere is empty (the simplex|
|                                           | is then said to be Gabriel), and as the minimum of the filtration    |
|                                           | values of the codimension 1 cofaces that make it not Gabriel         |
|                                           | otherwise. All simplices that have a filtration value strictly       |
|                                           | greater than a given alpha squared value are not inserted into the   |
|                                           | complex.                                                             |
|                                           |                                                                      |
|                                           | This package requires having CGAL version 4.7 or higher (4.8.1 is    |
|                                           | advised for better perfomances).                                     |
+-------------------------------------------+----------------------------------------------------------------------+
| :doc:`alpha_complex_user`                 | :doc:`alpha_complex_ref`                                             |
+-------------------------------------------+----------------------------------------------------------------------+
