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

#include "Persistent_cohomology_interface.h"
#include "Simplex_tree_interface.h"

#include <vector>
#include <utility>  // std::pair
#include <iostream>

namespace Gudhi {

namespace alpha_complex {

class Alpha_complex_interface {
  using Dynamic_kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag >;
  using Point_d = Dynamic_kernel::Point_d;
  typedef typename Simplex_tree<>::Simplex_handle Simplex_handle;
  typedef typename std::pair<Simplex_handle, bool> Insertion_result;
  using Simplex = std::vector<Vertex_handle>;
  using Filtered_complex = std::pair<Simplex, Filtration_value>;
  using Complex_tree = std::vector<Filtered_complex>;

  typedef typename Simplex_tree<>::Simplex_key Simplex_key;

 public:
  Alpha_complex_interface(std::vector<std::vector<double>>&points)
  : pcoh_(nullptr) {
    alpha_complex_ = new Alpha_complex<Dynamic_kernel>(points);
  }

  Alpha_complex_interface(std::string off_file_name, bool from_file = true)
  : pcoh_(nullptr) {
    alpha_complex_ = new Alpha_complex<Dynamic_kernel>(off_file_name);
  }

  bool find_simplex(const Simplex& vh) {
    return (simplex_tree_.find(vh) != simplex_tree_.null_simplex());
  }

  bool insert_simplex_and_subfaces(const Simplex& complex, Filtration_value filtration = 0) {
    Insertion_result result = simplex_tree_.insert_simplex_and_subfaces(complex, filtration);
    return (result.second);
  }

  Filtration_value simplex_filtration(const Simplex& complex) {
    return simplex_tree_.filtration(simplex_tree_.find(complex));
  }

  void remove_maximal_simplex(const Simplex& complex) {
    return simplex_tree_.remove_maximal_simplex(simplex_tree_.find(complex));
  }

  Complex_tree get_filtered_tree() {
    Complex_tree filtered_tree;
    for (auto f_simplex : simplex_tree_.filtration_simplex_range()) {
      Simplex simplex;
      for (auto vertex : simplex_tree_.simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      filtered_tree.push_back(std::make_pair(simplex, simplex_tree_.filtration(f_simplex)));
    }
    return filtered_tree;

  }

  Complex_tree get_skeleton_tree(int dimension) {
    Complex_tree skeleton_tree;
    for (auto f_simplex : simplex_tree_.skeleton_simplex_range(dimension)) {
      Simplex simplex;
      for (auto vertex : simplex_tree_.simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      skeleton_tree.push_back(std::make_pair(simplex, simplex_tree_.filtration(f_simplex)));
    }
    return skeleton_tree;
  }

  Complex_tree get_star_tree(const Simplex& complex) {
    Complex_tree star_tree;
    for (auto f_simplex : simplex_tree_.star_simplex_range(simplex_tree_.find(complex))) {
      Simplex simplex;
      for (auto vertex : simplex_tree_.simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      star_tree.push_back(std::make_pair(simplex, simplex_tree_.filtration(f_simplex)));
    }
    return star_tree;
  }

  Complex_tree get_coface_tree(const Simplex& complex, int dimension) {
    Complex_tree coface_tree;
    for (auto f_simplex : simplex_tree_.cofaces_simplex_range(simplex_tree_.find(complex), dimension)) {
      Simplex simplex;
      for (auto vertex : simplex_tree_.simplex_vertex_range(f_simplex)) {
        simplex.insert(simplex.begin(), vertex);
      }
      coface_tree.push_back(std::make_pair(simplex, simplex_tree_.filtration(f_simplex)));
    }
    return coface_tree;
  }

  // Specific to Witness complex because no inheritance
  Filtration_value filtration() const {
    return simplex_tree_.filtration();
  }

  void set_filtration(Filtration_value fil) {
    simplex_tree_.set_filtration(fil);
  }

  void initialize_filtration() {
    simplex_tree_.initialize_filtration();
  }

  size_t num_vertices() const {
    return simplex_tree_.num_vertices();
  }

  size_t num_simplices() {
    return simplex_tree_.num_simplices();
  }

  int dimension() const {
    return simplex_tree_.dimension();
  }

  void set_dimension(int dimension) {
    simplex_tree_.set_dimension(dimension);
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

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence(int homology_coeff_field, double min_persistence) {
    if (pcoh_ != nullptr) {
      delete pcoh_;
    }
    pcoh_ = new Persistent_cohomology_interface<Simplex_tree<>>(&simplex_tree_);
    return pcoh_->get_persistence(homology_coeff_field, min_persistence);
  }
  
  std::vector<int> get_betti_numbers() const {
    if (pcoh_ != nullptr) {
      return pcoh_->betti_numbers();
    }
    std::vector<int> betti_numbers;
    return betti_numbers;
  }

  std::vector<int> get_persistent_betti_numbers(Filtration_value from, Filtration_value to) const {
    if (pcoh_ != nullptr) {
      return pcoh_->persistent_betti_numbers(from, to);
    }
    std::vector<int> persistent_betti_numbers;
    return persistent_betti_numbers;
  }
  
  void create_simplex_tree(Simplex_tree_interface<>& simplex_tree, double max_alpha_square) {
    alpha_complex_->create_complex(simplex_tree, max_alpha_square);
    simplex_tree.initialize_filtration();
  }

 private:
  Simplex_tree<> simplex_tree_;
  Persistent_cohomology_interface<Simplex_tree<>>* pcoh_;
  Alpha_complex<Dynamic_kernel>* alpha_complex_;
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_INTERFACE_H
