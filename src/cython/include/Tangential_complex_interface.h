/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef TANGENTIAL_COMPLEX_INTERFACE_H
#define	TANGENTIAL_COMPLEX_INTERFACE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Tangential_complex.h>
#include <gudhi/Points_off_io.h>
#include <CGAL/Epick_d.h>

#include "Simplex_tree_interface.h"

#include <vector>
#include <utility>  // std::pair
#include <iostream>

namespace Gudhi {

namespace tangential_complex {

class Tangential_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;
  typedef typename Simplex_tree<>::Simplex_handle Simplex_handle;
  typedef typename std::pair<Simplex_handle, bool> Insertion_result;
  using Simplex = std::vector<Vertex_handle>;
  using Filtered_complex = std::pair<Simplex, Filtration_value>;
  using Complex_tree = std::vector<Filtered_complex>;
  using TC = Tangential_complex<Dynamic_kernel, CGAL::Dynamic_dimension_tag, CGAL::Parallel_tag>;

 public:
  Tangential_complex_interface(std::vector<std::vector<double>>&points) {
    Dynamic_kernel k;
    unsigned intrisic_dim = 0;
    if (points.size() > 0)
      intrisic_dim = points[0].size() - 1;
    
    tangential_complex_ = new TC(points, intrisic_dim, k);
    tangential_complex_->compute_tangential_complex();
    num_inconsistencies_ = tangential_complex_->number_of_inconsistent_simplices();
  }

  Tangential_complex_interface(std::string off_file_name, bool from_file = true) {
    Gudhi::Points_off_reader<Point_d> off_reader(off_file_name);
    Dynamic_kernel k;
    unsigned intrisic_dim = 0;
    std::vector<Point_d> points = off_reader.get_point_cloud();
    if (points.size() > 0)
      intrisic_dim = points[0].size() - 1;

    tangential_complex_ = new TC(points, intrisic_dim, k);
    tangential_complex_->compute_tangential_complex();
    num_inconsistencies_ = tangential_complex_->number_of_inconsistent_simplices();
  }

  std::vector<double> get_point(unsigned vh) {
    std::vector<double> vd;
    if (vh < tangential_complex_->number_of_vertices()) {
      Point_d ph = tangential_complex_->get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); coord++)
        vd.push_back(*coord);
    }
    return vd;
  }

  unsigned number_of_vertices() {
    return tangential_complex_->number_of_vertices();
  }

  unsigned number_of_simplices() {
    return num_inconsistencies_.num_simplices;
  }

  unsigned number_of_inconsistent_simplices() {
    return num_inconsistencies_.num_inconsistent_simplices;
  }

  unsigned number_of_inconsistent_stars() {
    return num_inconsistencies_.num_inconsistent_stars;
  }

  void fix_inconsistencies_using_perturbation(double max_perturb, double time_limit) {
      tangential_complex_->fix_inconsistencies_using_perturbation(max_perturb, time_limit);
      num_inconsistencies_ = tangential_complex_->number_of_inconsistent_simplices();
  }

  void create_simplex_tree(Simplex_tree<>* simplex_tree) {
    int max_dim = tangential_complex_->create_complex<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>(*simplex_tree);
    // FIXME
    simplex_tree->set_dimension(max_dim);
    simplex_tree->initialize_filtration();
  }

 private:
  TC* tangential_complex_;
  TC::Num_inconsistencies num_inconsistencies_;
};

}  // namespace tangential_complex

} // namespace Gudhi

#endif  // TANGENTIAL_COMPLEX_INTERFACE_H

