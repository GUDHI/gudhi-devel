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

#ifndef ALPHA_COMPLEX_INTERFACE_H
#define	ALPHA_COMPLEX_INTERFACE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Alpha_complex.h>
#include <CGAL/Epick_d.h>

#include "Simplex_tree_interface.h"

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace alpha_complex {

class Alpha_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;

 public:
  Alpha_complex_interface(const std::vector<std::vector<double>>& points) {
    alpha_complex_ = new Alpha_complex<Dynamic_kernel>(points);
  }

  Alpha_complex_interface(const std::string& off_file_name, bool from_file = true) {
    alpha_complex_ = new Alpha_complex<Dynamic_kernel>(off_file_name);
  }

  std::vector<double> get_point(int vh) {
    std::vector<double> vd;
    try {
      Point_d ph = alpha_complex_->get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); coord++)
        vd.push_back(*coord);
    } catch (std::out_of_range outofrange) {
      // std::out_of_range is thrown in case not found. Other exceptions must be re-thrown
    }
    return vd;
  }

  void create_simplex_tree(Simplex_tree_interface<>* simplex_tree, double max_alpha_square) {
    alpha_complex_->create_complex(*simplex_tree, max_alpha_square);
    simplex_tree->initialize_filtration();
  }

 private:
  Alpha_complex<Dynamic_kernel>* alpha_complex_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_INTERFACE_H
