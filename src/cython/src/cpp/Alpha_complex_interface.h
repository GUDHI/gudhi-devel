/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA Saclay (France)
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

#include <gudhi/Alpha_complex.h>
#include <CGAL/Epick_d.h>

#include <vector>
#include <utility>  // std::pair
#include <iostream>

namespace Gudhi {

namespace alphacomplex {

class Alpha_complex_interface : public Alpha_complex< CGAL::Epick_d< CGAL::Dynamic_dimension_tag > > {
  using Alpha_complex = Alpha_complex< CGAL::Epick_d< CGAL::Dynamic_dimension_tag > > ;
  typedef typename Alpha_complex::Simplex_handle Simplex_handle;
  typedef typename std::pair<Simplex_handle, bool> Insertion_result;
  using Simplex = std::vector<Vertex_handle>;
  using Filtered_complex = std::pair<Simplex, Filtration_value>;
  using Complex_tree = std::vector<Filtered_complex>;
  using Point_d = Alpha_complex::Point_d;

 public:

  Alpha_complex_interface(std::vector<std::vector<double>>&points, double max_alpha_square)
  : Alpha_complex(points, max_alpha_square) {
  }

  bool find_simplex(const Simplex& vh) {
    return (Alpha_complex::find(vh) != Alpha_complex::null_simplex());
  }

  bool insert_simplex_and_subfaces(const Simplex& complex, Filtration_value filtration = 0) {
    Insertion_result result = Alpha_complex::insert_simplex_and_subfaces(complex, filtration);
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& complex) {
    return Alpha_complex::filtration(Alpha_complex::find(complex));
  }

  void remove_maximal_simplex(const Simplex& complex) {
    return Alpha_complex::remove_maximal_simplex(Alpha_complex::find(complex));
  }

  Complex_tree get_filtered_tree() {
    Complex_tree filtered_tree;
    for (auto f_simplex : Alpha_complex::filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : Alpha_complex::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtered_tree.push_back(std::make_pair(simplex, Alpha_complex::filtration(f_simplex)));
    }
    return filtered_tree;

  }

  Complex_tree get_skeleton_tree(int dimension) {
    Complex_tree skeleton_tree;
    for (auto f_simplex : Alpha_complex::skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : Alpha_complex::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeleton_tree.push_back(std::make_pair(simplex, Alpha_complex::filtration(f_simplex)));
    }
    return skeleton_tree;
  }

  Complex_tree get_star_tree(const Simplex& complex) {
    Complex_tree star_tree;
    for (auto f_simplex : Alpha_complex::star_simplex_range(Alpha_complex::find(complex))) {
      Simplex simplex;
      for (auto vertex : Alpha_complex::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      star_tree.push_back(std::make_pair(simplex, Alpha_complex::filtration(f_simplex)));
    }
    return star_tree;
  }

  Complex_tree get_coface_tree(const Simplex& complex, int dimension) {
    Complex_tree coface_tree;
    for (auto f_simplex : Alpha_complex::cofaces_simplex_range(Alpha_complex::find(complex), dimension)) {
      Simplex simplex;
      for (auto vertex : Alpha_complex::simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      coface_tree.push_back(std::make_pair(simplex, Alpha_complex::filtration(f_simplex)));
    }
    return coface_tree;
  }
  
  std::vector<double> get_point(int vh) {
    std::vector<double> vd;
    try {
      Point_d ph = Alpha_complex::get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); coord++)
        vd.push_back(*coord);
    } catch (std::out_of_range outofrange) {
      // std::out_of_range is thrown in case not found. Other exceptions must be re-thrown
    }
    return vd;
  }

};

}  // namespace alphacomplex

} // namespace Gudhi

#endif  // ALPHA_COMPLEX_INTERFACE_H

