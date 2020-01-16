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

namespace Gudhi {

template<class FilteredComplex>
class Persistent_cohomology_interface : public
persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> {
 private:
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
  Persistent_cohomology_interface(FilteredComplex* stptr)
      : persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp>(*stptr),
      stptr_(stptr) { }

  Persistent_cohomology_interface(FilteredComplex* stptr, bool persistence_dim_max)
      : persistent_cohomology::Persistent_cohomology<FilteredComplex,
          persistent_cohomology::Field_Zp>(*stptr, persistence_dim_max),
        stptr_(stptr) { }

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence(int homology_coeff_field,
                                                                         double min_persistence) {
    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::init_coefficients(homology_coeff_field);
    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::compute_persistent_cohomology(min_persistence);

    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp(stptr_);
    auto persistent_pairs = persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::get_persistent_pairs();
    std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);

    std::vector<std::pair<int, std::pair<double, double>>> persistence;
    for (auto pair : persistent_pairs) {
      persistence.push_back(std::make_pair(stptr_->dimension(get<0>(pair)),
                                           std::make_pair(stptr_->filtration(get<0>(pair)),
                                                          stptr_->filtration(get<1>(pair)))));
    }
    return persistence;
  }

  std::vector<std::pair<int, std::pair<std::pair<double, std::vector<int>>, std::pair<double,std::vector<int>>>>> get_persistence_generators(int homology_coeff_field,
                                                                                                                                             double min_persistence) {
    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::init_coefficients(homology_coeff_field);
    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::compute_persistent_cohomology(min_persistence);

    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp(stptr_);
    auto persistent_pairs = persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::get_persistent_pairs();
    std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);

    std::vector<std::pair<int, std::pair<std::pair<double, std::vector<int>>, std::pair<double,std::vector<int>>>>> persistence;
    for (auto pair : persistent_pairs) {
      std::vector<int> splx0, splx1;
      for (auto vertex : stptr_->simplex_vertex_range(get<0>(pair))){splx0.push_back(vertex);}
      if (isfinite(stptr_->filtration(get<1>(pair)))){ for (auto vertex : stptr_->simplex_vertex_range(get<1>(pair))){splx1.push_back(vertex);}}
      persistence.push_back(std::make_pair(stptr_->dimension(get<0>(pair)), std::make_pair(std::make_pair(stptr_->filtration(get<0>(pair)), splx0), std::make_pair(stptr_->filtration(get<1>(pair)), splx1))));
    }
    return persistence;
  }

  void top_dimensional_cofaces(std::vector<int> & cofaces, int splx){
    if (stptr_->dimension(stptr_->simplex(splx)) == stptr_->dimension()){cofaces.push_back(stptr_->simplex(splx));}
    else{  for (auto v : stptr_->coboundary_simplex_range(stptr_->simplex(splx))){top_dimensional_cofaces(cofaces, stptr_->key(v));}  }
  }

  std::vector<std::pair<int, std::pair<std::pair<double, int>, std::pair<double, int>>>> get_persistence_cubical_generators(int homology_coeff_field,
                                                                                                                            double min_persistence) {

    // Gather all top-dimensional cells and store their simplex handles
    std::vector<int> max_splx; for (auto splx : stptr_->filtration_simplex_range()){ if (stptr_->dimension(splx) == stptr_->dimension())  max_splx.push_back(splx); }
    // Sort these simplex handles and compute the ordering function
    // This function allows to go directly from the simplex handle to the position of the corresponding top-dimensional cell in the input data 
    std::map<int, int> order; std::sort(max_splx.begin(), max_splx.end()); for (int i = 0; i < max_splx.size(); i++)  order.insert(std::make_pair(max_splx[i], i));

    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::init_coefficients(homology_coeff_field);
    persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::compute_persistent_cohomology(min_persistence);

    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp(stptr_);
    auto persistent_pairs = persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::get_persistent_pairs();
    std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);

    std::vector<std::pair<int, std::pair<std::pair<double, int>, std::pair<double, int>>>> persistence;
    for (auto pair : persistent_pairs) {

      int splx0, splx1;

      double f0 = stptr_->filtration(get<0>(pair));
      // Recursively get the top-dimensional cells / cofaces associated to the persistence generator
      std::vector<int> faces0; top_dimensional_cofaces(faces0, stptr_->key(get<0>(pair)));
      // Find the top-dimensional cell / coface with the same filtration value
      int cf; for (int i = 0; i < faces0.size(); i++){ if (stptr_->filtration(faces0[i]) == f0){cf = i; break;}}
      // Retrieve the index of the corresponding top-dimensional cell in the input data
      splx0 = order[faces0[cf]];

      if (isfinite(stptr_->filtration(get<1>(pair)))){
      double f1 = stptr_->filtration(get<1>(pair));
      // Recursively get the top-dimensional cells / cofaces associated to the persistence generator
      std::vector<int> faces1; top_dimensional_cofaces(faces1, stptr_->key(get<1>(pair)));
      // Find the top-dimensional cell / coface with the same filtration value
      int cf; for (int i = 0; i < faces0.size(); i++){ if (stptr_->filtration(faces0[i]) == f0){cf = i; break;}}
      // Retrieve the index of the corresponding top-dimensional cell in the input data
      splx1 = order[faces1[cf]];
      }

      persistence.push_back(std::make_pair(stptr_->dimension(get<0>(pair)), std::make_pair(std::make_pair(stptr_->filtration(get<0>(pair)), splx0), std::make_pair(stptr_->filtration(get<1>(pair)), splx1))));
    }
    return persistence;
  }

  std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs() {
    auto pairs = persistent_cohomology::Persistent_cohomology<FilteredComplex,
      persistent_cohomology::Field_Zp>::get_persistent_pairs();

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
