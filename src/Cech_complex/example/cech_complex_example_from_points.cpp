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

#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <string>
#include <vector>
#include <array>

int main() {
  // Type definitions
  using Point_cloud = std::vector<std::array<double, 2>>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;

  Point_cloud points;
  points.push_back({1.0, 1.0});
  points.push_back({7.0, 0.0});
  points.push_back({4.0, 6.0});
  points.push_back({9.0, 6.0});
  points.push_back({0.0, 14.0});
  points.push_back({2.0, 19.0});
  points.push_back({9.0, 17.0});

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from points
  // ----------------------------------------------------------------------------
  // 7.1 is a magic number to force one blocker, and one non-blocker
  Filtration_value threshold = 7.1;
  Cech_complex cech_complex_from_points(points, threshold, Gudhi::Euclidean_distance());

  Simplex_tree stree;
  cech_complex_from_points.create_complex(stree, 2);
  // ----------------------------------------------------------------------------
  // Display information about the one skeleton Cech complex
  // ----------------------------------------------------------------------------
  std::cout << "Cech complex is of dimension " << stree.dimension() <<
               " - " << stree.num_simplices() << " simplices - " <<
               stree.num_vertices() << " vertices." << std::endl;

  std::cout << "Iterator on Cech complex simplices in the filtration order, with [filtration value]:" <<
               std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  return 0;
}
