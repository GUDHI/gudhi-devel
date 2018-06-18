/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALPHA_COMPLEX_3D_OPTIONS_H_
#define ALPHA_COMPLEX_3D_OPTIONS_H_


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>


namespace Gudhi {

namespace alpha_complex {

class Alpha_shapes_3d {
private:
  // Alpha_shape_3 templates type definitions
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel>;
  using Fb = CGAL::Alpha_shape_cell_base_3<Kernel>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb>;
  using Triangulation_3 = CGAL::Delaunay_triangulation_3<Kernel, Tds>;

public:
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3>;
  using Point_3 = Kernel::Point_3;

  static const bool exact = false;
  static const bool weighted = false;
  static const bool periodic = false;

};

class Exact_alpha_shapes_3d {
private:
  // Alpha_shape_3 templates type definitions
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Exact_tag = CGAL::Tag_true;
  using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel, CGAL::Default, Exact_tag>;
  using Fb = CGAL::Alpha_shape_cell_base_3<Kernel, CGAL::Default, Exact_tag>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb>;
  using Triangulation_3 = CGAL::Delaunay_triangulation_3<Kernel, Tds>;

public:
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3, Exact_tag>;
  using Point_3 = Kernel::Point_3;

  static const bool exact = true;
  static const bool weighted = false;
  static const bool periodic = false;
};

class Weighted_alpha_shapes_3d {
private:
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Rvb = CGAL::Regular_triangulation_vertex_base_3<Kernel>;
  using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel, Rvb>;
  using Rcb = CGAL::Regular_triangulation_cell_base_3<Kernel>;
  using Cb = CGAL::Alpha_shape_cell_base_3<Kernel, Rcb>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Triangulation_3 = CGAL::Regular_triangulation_3<Kernel, Tds>;


public:
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3>;
  using Point_3 = Triangulation_3::Bare_point;
  using Weighted_point_3 = Triangulation_3::Weighted_point;

  static const bool exact = false;
  static const bool weighted = true;
  static const bool periodic = false;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_3D_H_
