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
   typedef persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> Base;
  /*
   * Compare two intervals by dimension, then by length.
   */
  struct cmp_intervals_by_dim_then_length {
    template<typename Persistent_interval>
    bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
      if (std::get<0>(p1) == std::get<0>(p2)) {
        auto& i1 = std::get<1>(p1);
        auto& i2 = std::get<1>(p2);
        return std::get<1>(i1) - std::get<0>(i1) > std::get<1>(i2) - std::get<0>(i2);
      }
      else
        return (std::get<0>(p1) > std::get<0>(p2));
        // Why does this sort by decreasing dimension?
    }
  };

 public:
  Persistent_cohomology_interface(FilteredComplex* stptr)
      : Base(*stptr),
      stptr_(stptr) { }

  Persistent_cohomology_interface(FilteredComplex* stptr, bool persistence_dim_max)
      : Base(*stptr, persistence_dim_max),
        stptr_(stptr) { }

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence(int homology_coeff_field,
                                                                         double min_persistence) {
    Base::init_coefficients(homology_coeff_field);
    Base::compute_persistent_cohomology(min_persistence);

    auto const& persistent_pairs = Base::get_persistent_pairs();
    std::vector<std::pair<int, std::pair<double, double>>> persistence;
    persistence.reserve(persistent_pairs.size());
    for (auto pair : persistent_pairs) {
      persistence.emplace_back(stptr_->dimension(get<0>(pair)),
                               std::make_pair(stptr_->filtration(get<0>(pair)),
                                              stptr_->filtration(get<1>(pair))));
    }
    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp;
    std::sort(std::begin(persistence), std::end(persistence), cmp);
    return persistence;
  }

  std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs() {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs;
    auto const& pairs = Base::get_persistent_pairs();
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
        death.reserve(birth.size()+1);
        for (auto vertex : stptr_->simplex_vertex_range(get<1>(pair))) {
          death.push_back(vertex);
        }
      }

      persistence_pairs.emplace_back(std::move(birth), std::move(death));
    }
    return persistence_pairs;
  }

  // TODO: (possibly at the python level)
  // - an option to return only some of those vectors?
  typedef std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> Generators;

  Generators lower_star_generators(double min_persistence) {
    Generators out;
    // diags[i] should be interpreted as vector<array<int,2>>
    auto& diags = out.first;
    // diagsinf[i] should be interpreted as vector<int>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      if(stptr_->filtration(t) - stptr_->filtration(s) <= min_persistence)
        continue;
      int dim = stptr_->dimension(s);
      auto v = stptr_->vertex_with_same_filtration(s);
      if(t == stptr_->null_simplex()) {
        while(diagsinf.size() < dim+1) diagsinf.emplace_back();
        diagsinf[dim].push_back(v);
      } else {
        while(diags.size() < dim+1) diags.emplace_back();
        auto w = stptr_->vertex_with_same_filtration(t);
        diags[dim].push_back(v);
        diags[dim].push_back(w);
      }
    }
    return out;
  }

  // An alternative, to avoid those different sizes, would be to "pad" vertex generator v as (v, v) or (v, -1). When using it as index, this corresponds to adding the vertex filtration values either on the diagonal of the distance matrix, or as an extra row or column.
  // We could also merge the vectors for different dimensions into a single one, with an extra column for the dimension (converted to type double).
  Generators flag_generators(double min_persistence) {
    Generators out;
    // diags[0] should be interpreted as vector<array<int,3>> and other diags[i] as vector<array<int,4>>
    auto& diags = out.first;
    // diagsinf[0] should be interpreted as vector<int> and other diagsinf[i] as vector<array<int,2>>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      if(stptr_->filtration(t) - stptr_->filtration(s) <= min_persistence)
        continue;
      int dim = stptr_->dimension(s);
      bool infinite = t == stptr_->null_simplex();
      if(infinite) {
        if(dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if(diagsinf.size()==0)diagsinf.emplace_back();
          diagsinf[0].push_back(v);
        } else {
          auto e = stptr_->edge_with_same_filtration(s);
          auto&& e_vertices = stptr_->simplex_vertex_range(e);
          auto i = std::begin(e_vertices);
          auto v1 = *i;
          auto v2 = *++i;
          GUDHI_CHECK(++i==std::end(e_vertices), "must be an edge");
          while(diagsinf.size() < dim+1) diagsinf.emplace_back();
          diagsinf[dim].push_back(v1);
          diagsinf[dim].push_back(v2);
        }
      } else {
        auto et = stptr_->edge_with_same_filtration(t);
        auto&& et_vertices = stptr_->simplex_vertex_range(et);
        auto it = std::begin(et_vertices);
        auto w1 = *it;
        auto w2 = *++it;
        GUDHI_CHECK(++it==std::end(et_vertices), "must be an edge");
        if(dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if(diags.size()==0)diags.emplace_back();
          diags[0].push_back(v);
          diags[0].push_back(w1);
          diags[0].push_back(w2);
        } else {
          auto es = stptr_->edge_with_same_filtration(s);
          auto&& es_vertices = stptr_->simplex_vertex_range(es);
          auto is = std::begin(es_vertices);
          auto v1 = *is;
          auto v2 = *++is;
          GUDHI_CHECK(++is==std::end(es_vertices), "must be an edge");
          while(diags.size() < dim+1) diags.emplace_back();
          diags[dim].push_back(v1);
          diags[dim].push_back(v2);
          diags[dim].push_back(w1);
          diags[dim].push_back(w2);
        }
      }
    }
    return out;
  }

 private:
  // A copy
  FilteredComplex* stptr_;
};

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
