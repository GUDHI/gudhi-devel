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

#ifndef INCLUDE_RIPS_COMPLEX_INTERFACE_H_
#define INCLUDE_RIPS_COMPLEX_INTERFACE_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <utility>  // std::pair
#include <string>

namespace Gudhi {

namespace rips_complex {

class Rips_complex_interface {
  using Point_d = std::vector<double>;
  using Distance_matrix = std::vector<std::vector<Simplex_tree_interface<>::Filtration_value>>;

 public:
  Rips_complex_interface(const std::vector<std::vector<double>>& values, double threshold, bool euclidean) {
    if (euclidean) {
      // Rips construction where values is a vector of points
      rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(values, threshold,
                                                                                   Gudhi::Euclidean_distance());
    } else {
      // Rips construction where values is a distance matrix
      rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(values, threshold);
    }
  }

  Rips_complex_interface(const std::string& file_name, double threshold, bool euclidean, bool from_file = true) {
    if (euclidean) {
      // Rips construction where file_name is an OFF file
      Gudhi::Points_off_reader<Point_d> off_reader(file_name);
      rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(off_reader.get_point_cloud(),
                                                                                   threshold,
                                                                                   Gudhi::Euclidean_distance());
    } else {
      // Rips construction where values is a distance matrix
      Distance_matrix distances =
          Gudhi::read_lower_triangular_matrix_from_csv_file<Simplex_tree_interface<>::Filtration_value>(file_name);
      rips_complex_ = new Rips_complex<Simplex_tree_interface<>::Filtration_value>(distances, threshold);
    }
  }

  ~Rips_complex_interface() {
    delete rips_complex_;
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

#endif  // INCLUDE_RIPS_COMPLEX_INTERFACE_H_
