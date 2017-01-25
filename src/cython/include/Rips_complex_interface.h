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

#ifndef RIPS_COMPLEX_INTERFACE_H
#define	RIPS_COMPLEX_INTERFACE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <string>

namespace Gudhi {

namespace rips_complex {

class Rips_complex_interface {
  using Point_d = std::vector<double>;

 public:
  Rips_complex_interface(std::vector<std::vector<double>>&points, double threshold) {
    rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(points, threshold,
                                                                                 Euclidean_distance());
  }

  Rips_complex_interface(std::string off_file_name, double threshold, bool from_file = true) {
    Gudhi::Points_off_reader<Point_d> off_reader(off_file_name);
    rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(off_reader.get_point_cloud(),
                                                                                 threshold, Euclidean_distance());
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, int dim_max) {
    rips_complex_->create_complex(*simplex_tree, dim_max);
    simplex_tree->initialize_filtration();
  }

 private:
  Rips_complex<Simplex_tree_interface<>::Filtration_value>* rips_complex_;
};

}  // namespace rips_complex

}  // namespace Gudhi

#endif  // RIPS_COMPLEX_INTERFACE_H
