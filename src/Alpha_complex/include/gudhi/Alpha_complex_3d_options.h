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
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>


namespace Gudhi {

namespace alpha_complex {

class Alpha_shapes_3d {
private:
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Alpha_shape_vertex_base_3<Kernel>;
  using Fb = CGAL::Alpha_shape_cell_base_3<Kernel>;
  using Tds = CGAL::Triangulation_data_structure_3<Vb, Fb>;
  using Triangulation_3 = CGAL::Delaunay_triangulation_3<Kernel, Tds>;

public:
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Triangulation_3>;
  using Point_3 = Kernel::Point_3;

  static const bool weighted = false;
  static const bool periodic = false;

  template<class Filtration_value, class Alpha_value_iterator>
  static Filtration_value value_from_iterator(const Alpha_value_iterator avi){
    return /*std::sqrt*/ *avi;
  }
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

  static const bool weighted = false;
  static const bool periodic = false;
  static const bool exact = true;

  template<class Filtration_value, class Alpha_value_iterator>
  static Filtration_value value_from_iterator(const Alpha_value_iterator avi){
    return /*std::sqrt*/ CGAL::to_double(avi->exact());
  }
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

  static const bool weighted = true;
  static const bool periodic = false;
  static const bool exact = false;

  template<class Filtration_value, class Alpha_value_iterator>
  static Filtration_value value_from_iterator(const Alpha_value_iterator avi){
    return /*std::sqrt*/ *avi;
  }
};

class Periodic_alpha_shapes_3d {
private:
  // Traits
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using PK = CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>;
// Vertex type
  using DsVb = CGAL::Periodic_3_triangulation_ds_vertex_base_3<>;
  using Vb = CGAL::Triangulation_vertex_base_3<PK, DsVb>;
  using AsVb = CGAL::Alpha_shape_vertex_base_3<PK, Vb>;
// Cell type
  using DsCb = CGAL::Periodic_3_triangulation_ds_cell_base_3<>;
  using Cb = CGAL::Triangulation_cell_base_3<PK, DsCb>;
  using AsCb = CGAL::Alpha_shape_cell_base_3<PK, Cb>;
  using Tds = CGAL::Triangulation_data_structure_3<AsVb, AsCb>;

public:
  using Periodic_delaunay_triangulation_3 = CGAL::Periodic_3_Delaunay_triangulation_3<PK, Tds>;
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Periodic_delaunay_triangulation_3>;
  using Point_3 = PK::Point_3;
  using Iso_cuboid_3 = PK::Iso_cuboid_3;

  static const bool weighted = false;
  static const bool periodic = true;
  static const bool exact = false;

  template<class Filtration_value, class Alpha_value_iterator>
  static Filtration_value value_from_iterator(const Alpha_value_iterator avi){
    return /*std::sqrt*/ *avi;
  }
};

class Weighted_periodic_alpha_shapes_3d {
private:
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Periodic_kernel = CGAL::Periodic_3_regular_triangulation_traits_3<Kernel>;
  using DsVb = CGAL::Periodic_3_triangulation_ds_vertex_base_3<>;
  using Vb = CGAL::Regular_triangulation_vertex_base_3<Periodic_kernel, DsVb>;
  using AsVb = CGAL::Alpha_shape_vertex_base_3<Periodic_kernel, Vb>;
  using DsCb = CGAL::Periodic_3_triangulation_ds_cell_base_3<>;
  using Cb = CGAL::Regular_triangulation_cell_base_3<Periodic_kernel, DsCb>;
  using AsCb = CGAL::Alpha_shape_cell_base_3<Periodic_kernel, Cb>;
  using Tds = CGAL::Triangulation_data_structure_3<AsVb, AsCb>;

public:
  using Periodic_delaunay_triangulation_3 = CGAL::Periodic_3_regular_triangulation_3<Periodic_kernel, Tds>;
  using Alpha_shape_3 = CGAL::Alpha_shape_3<Periodic_delaunay_triangulation_3>;
  using Point_3 = Periodic_delaunay_triangulation_3::Bare_point;
  using Weighted_point_3 = Periodic_delaunay_triangulation_3::Weighted_point;
  using Iso_cuboid_3 = Periodic_kernel::Iso_cuboid_3;

  static const bool weighted = true;
  static const bool periodic = true;
  static const bool exact = false;

  template<class Filtration_value, class Alpha_value_iterator>
  static Filtration_value value_from_iterator(const Alpha_value_iterator avi){
    return /*std::sqrt*/ *avi;
  }
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_3D_OPTIONS_H_
