/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_TANGENTIAL_COMPLEX_INTERFACE_H_
#define INCLUDE_TANGENTIAL_COMPLEX_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Tangential_complex.h>
#include <gudhi/Points_off_io.h>
#include <CGAL/Epick_d.h>

#include "Simplex_tree_interface.h"

#include <vector>
#include <utility>  // std::pair
#include <iostream>
#include <string>

namespace Gudhi {

namespace tangential_complex {

class Tangential_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;
  using TC = Tangential_complex<Dynamic_kernel, CGAL::Dynamic_dimension_tag, CGAL::Parallel_tag>;

 public:
  Tangential_complex_interface(int intrisic_dim, const std::vector<std::vector<double>>& points) {
    Dynamic_kernel k;

    tangential_complex_ = new TC(points, intrisic_dim, k);
  }

  Tangential_complex_interface(int intrisic_dim, const std::string& off_file_name, bool from_file = true) {
    Dynamic_kernel k;

    Gudhi::Points_off_reader<Point_d> off_reader(off_file_name);
    std::vector<Point_d> points = off_reader.get_point_cloud();

    tangential_complex_ = new TC(points, intrisic_dim, k);
  }

  ~Tangential_complex_interface() {
    delete tangential_complex_;
  }

  void compute_tangential_complex() {
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
    tangential_complex_->create_complex<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>(*simplex_tree);
    simplex_tree->initialize_filtration();
  }

  void set_max_squared_edge_length(double max_squared_edge_length) {
    tangential_complex_->set_max_squared_edge_length(max_squared_edge_length);
  }

private:
  TC* tangential_complex_;
  TC::Num_inconsistencies num_inconsistencies_;
};

}  // namespace tangential_complex

}  // namespace Gudhi

#endif  // INCLUDE_TANGENTIAL_COMPLEX_INTERFACE_H_
