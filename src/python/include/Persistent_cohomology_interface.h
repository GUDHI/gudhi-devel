/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
#define INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_

#include <gudhi/Persistent_cohomology.h>

#include <vector>
#include <utility>  // for std::pair
#include <algorithm>  // for sort
#include <unordered_map>

namespace Gudhi {

template<class FilteredComplex>
class Persistent_cohomology_interface : public
persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> {
 private:
   typedef persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> Base;
  /*
   * Compare two intervals by dimension, then by length.
   */
  struct cmp_intervals_by_dim_then_length {
    explicit cmp_intervals_by_dim_then_length(FilteredComplex * sc)
        : sc_(sc) { }

    template<typename Persistent_interval>
    bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
      if (sc_->dimension(get < 0 > (p1)) == sc_->dimension(get < 0 > (p2)))
        return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
                > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
      else
        return (sc_->dimension(get < 0 > (p1)) > sc_->dimension(get < 0 > (p2)));
    }
    FilteredComplex* sc_;
  };

 public:
  Persistent_cohomology_interface(FilteredComplex* stptr, bool persistence_dim_max=false)
      : Base(*stptr, persistence_dim_max),
        stptr_(stptr) { }

  // TODO: move to the constructors?
  void compute_persistence(int homology_coeff_field, double min_persistence) {
    Base::init_coefficients(homology_coeff_field);
    Base::compute_persistent_cohomology(min_persistence);
  }

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence() {
    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp(stptr_);
    auto persistent_pairs = Base::get_persistent_pairs();
    std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
    std::vector<std::pair<int, std::pair<double, double>>> persistence;
    for (auto pair : persistent_pairs) {
      persistence.push_back(std::make_pair(stptr_->dimension(get<0>(pair)),
                                           std::make_pair(stptr_->filtration(get<0>(pair)),
                                                          stptr_->filtration(get<1>(pair)))));
    }
    return persistence;
  }

  std::vector<std::vector<int>> cofaces_of_cubical_persistence_pairs() {

    // Warning: this function is meant to be used with CubicalComplex only!!

    auto pairs = persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::get_persistent_pairs();

    // Gather all top-dimensional cells and store their simplex handles
    std::vector<int> max_splx; for (auto splx : stptr_->top_dimensional_cells_range()){
      max_splx.push_back(splx); 
    }
    // Sort these simplex handles and compute the ordering function
    // This function allows to go directly from the simplex handle to the position of the corresponding top-dimensional cell in the input data 
    std::unordered_map<int, int> order;  
    //std::sort(max_splx.begin(), max_splx.end()); 
    for (unsigned int i = 0; i < max_splx.size(); i++)  order.emplace(max_splx[i], i);

    std::vector<std::vector<int>> persistence_pairs;
    for (auto pair : pairs) {
      int h = stptr_->dimension(get<0>(pair));
      // Recursively get the top-dimensional cell / coface associated to the persistence generator
      int face0 = stptr_->get_top_dimensional_coface_of_a_cell(get<0>(pair));
      // Retrieve the index of the corresponding top-dimensional cell in the input data
      int splx0 = order[face0];

      int splx1 = -1;
      if (isfinite(stptr_->filtration(get<1>(pair)))){
      // Recursively get the top-dimensional cell / coface associated to the persistence generator
      int face1 = stptr_->get_top_dimensional_coface_of_a_cell(get<1>(pair));
      // Retrieve the index of the corresponding top-dimensional cell in the input data
      splx1 = order[face1];
      }
      std::vector<int> vect{ h, splx0, splx1};
      persistence_pairs.push_back(vect);
    }
    return persistence_pairs;
  }

  std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs() {
    auto pairs = Base::get_persistent_pairs();

    std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs;
    persistence_pairs.reserve(pairs.size());
    for (auto pair : pairs) {
      std::vector<int> birth;
      if (get<0>(pair) != stptr_->null_simplex()) {
        for (auto vertex : stptr_->simplex_vertex_range(get<0>(pair))) {
          birth.push_back(vertex);
        }
      }

      std::vector<int> death;
      if (get<1>(pair) != stptr_->null_simplex()) {
        for (auto vertex : stptr_->simplex_vertex_range(get<1>(pair))) {
          death.push_back(vertex);
        }
      }

      persistence_pairs.push_back(std::make_pair(birth, death));
    }
    return persistence_pairs;
  }

 private:
  // A copy
  FilteredComplex* stptr_;
};

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
