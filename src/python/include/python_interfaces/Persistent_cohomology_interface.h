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

#include <cstdlib>
#include <vector>
#include <utility>    // for std::pair
#include <algorithm>  // for sort
#include <unordered_map>

#include <nanobind/nanobind.h>

#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Simplex_tree.h>  // for Extended_simplex_type

namespace Gudhi {

template <class FilteredComplex>
class Persistent_cohomology_interface
    : public persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp>
{
 private:
  typedef persistent_cohomology::Persistent_cohomology<FilteredComplex, persistent_cohomology::Field_Zp> Base;

  /*
   * Compare two intervals by dimension, then by length.
   */
  struct cmp_intervals_by_dim_then_length {
    template <typename Persistent_interval>
    bool operator()(const Persistent_interval& p1, const Persistent_interval& p2)
    {
      if (std::get<0>(p1) == std::get<0>(p2)) {
        auto& i1 = std::get<1>(p1);
        auto& i2 = std::get<1>(p2);
        return std::get<1>(i1) - std::get<0>(i1) > std::get<1>(i2) - std::get<0>(i2);
      } else
        return (std::get<0>(p1) > std::get<0>(p2));
      // Why does this sort by decreasing dimension?
    }
  };

 public:
  Persistent_cohomology_interface(FilteredComplex& stptr, bool persistence_dim_max = false)
      : Base(stptr, persistence_dim_max), stptr_(&stptr)
  {}

  // TODO: move to the constructors?
  void compute_persistence(int homology_coeff_field, double min_persistence)
  {
    Base::init_coefficients(homology_coeff_field);
    Base::compute_persistent_cohomology(min_persistence);
  }

  std::vector<std::pair<int, std::pair<double, double>>> get_persistence()
  {
    std::vector<std::pair<int, std::pair<double, double>>> persistence;
    auto const& persistent_pairs = Base::get_persistent_pairs();
    persistence.reserve(persistent_pairs.size());
    for (auto pair : persistent_pairs) {
      persistence.emplace_back(stptr_->dimension(get<0>(pair)),
                               std::make_pair(stptr_->filtration(get<0>(pair)), stptr_->filtration(get<1>(pair))));
    }
    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp;
    std::sort(std::begin(persistence), std::end(persistence), cmp);
    return persistence;
  }

  // This function computes the top-dimensional cofaces associated to the positive and negative
  // simplices of a cubical complex. The output format is a vector of vectors of three integers,
  // which are [homological dimension, index of top-dimensional coface of positive simplex,
  // index of top-dimensional coface of negative simplex]. If the topological feature is essential,
  // then the index of top-dimensional coface of negative simplex is arbitrarily set to -1.
  std::vector<std::vector<int>> cofaces_of_cubical_persistence_pairs()
  {
    // Warning: this function is meant to be used with CubicalComplex only!!

    auto&& pairs = Base::get_persistent_pairs();

    // Compute the ordering function of the top-dimensional cells simplex handles
    // This function allows to go directly from the simplex handle to the position of the corresponding top-dimensional
    // cell in the input data
    std::unordered_map<std::size_t, int> order;
    unsigned idx = 0;
    for (auto splx : stptr_->top_dimensional_cells_range()) {
      order.emplace(splx, idx);
      idx++;
    }

    std::vector<std::vector<int>> persistence_pairs;
    {
      // mostly just in case the coface methods or others end up using parallel processing one day.
      nanobind::gil_scoped_release release;

      for (auto pair : pairs) {
        int h = static_cast<int>(stptr_->dimension(get<0>(pair)));
        // Recursively get the top-dimensional cell / coface associated to the persistence generator
        std::size_t face0 = stptr_->get_top_dimensional_coface_of_a_cell(get<0>(pair));
        // Retrieve the index of the corresponding top-dimensional cell in the input data
        int splx0 = order[face0];

        int splx1 = -1;
        if (get<1>(pair) != stptr_->null_simplex()) {
          // Recursively get the top-dimensional cell / coface associated to the persistence generator
          std::size_t face1 = stptr_->get_top_dimensional_coface_of_a_cell(get<1>(pair));
          // Retrieve the index of the corresponding top-dimensional cell in the input data
          splx1 = order[face1];
        }
        persistence_pairs.push_back({h, splx0, splx1});
      }
    }
    return persistence_pairs;
  }

  // This function computes the vertices associated to the positive and negative
  // simplices of a cubical complex. The output format is a vector of vectors of three integers,
  // which are [homological dimension, index of vertex of positive simplex,
  // index of vertex of negative simplex]. If the topological feature is essential,
  // then the index of vertex of negative simplex is arbitrarily set to -1.
  std::vector<std::vector<int>> vertices_of_cubical_persistence_pairs()
  {
    // Warning: this function is meant to be used with CubicalComplex only!!
    auto&& pairs = Base::get_persistent_pairs();

    // Compute the ordering function of the vertices simplex handles
    // This function allows to go directly from the simplex handle to the position of the corresponding vertex in the
    // input data
    std::unordered_map<std::size_t, int> order;
    unsigned idx = 0;
    for (auto splx : stptr_->vertices_range()) {
      order.emplace(splx, idx);
      idx++;
    }

    std::vector<std::vector<int>> persistence_pairs;
    {
      // mostly just in case the vertex methods or others end up using parallel processing one day.
      nanobind::gil_scoped_release release;

      for (auto pair : pairs) {
        int h = static_cast<int>(stptr_->dimension(get<0>(pair)));
        // Recursively get the vertex associated to the persistence generator
        std::size_t face0 = stptr_->get_vertex_of_a_cell(get<0>(pair));
        // Retrieve the index of the corresponding vertex in the input data
        int splx0 = order[face0];

        int splx1 = -1;
        if (get<1>(pair) != stptr_->null_simplex()) {
          // Recursively get the vertex associated to the persistence generator
          std::size_t face1 = stptr_->get_vertex_of_a_cell(get<1>(pair));
          // Retrieve the index of the corresponding vertex in the input data
          splx1 = order[face1];
        }
        persistence_pairs.push_back({h, splx0, splx1});
      }
    }
    return persistence_pairs;
  }

  std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs()
  {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> persistence_pairs;
    auto const& pairs = Base::get_persistent_pairs();
    persistence_pairs.reserve(pairs.size());
    std::vector<int> birth;
    std::vector<int> death;
    for (auto pair : pairs) {
      birth.clear();
      if (get<0>(pair) != stptr_->null_simplex()) {
        for (auto vertex : stptr_->simplex_vertex_range(get<0>(pair))) {
          birth.push_back(vertex);
        }
      }

      death.clear();
      if (get<1>(pair) != stptr_->null_simplex()) {
        death.reserve(birth.size() + 1);
        for (auto vertex : stptr_->simplex_vertex_range(get<1>(pair))) {
          death.push_back(vertex);
        }
      }

      persistence_pairs.emplace_back(birth, death);
    }
    return persistence_pairs;
  }

  // TODO: (possibly at the python level)
  // - an option to return only some of those vectors?
  typedef std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> Generators;

  Generators lower_star_generators()
  {
    Generators out;
    // diags[i] should be interpreted as vector<array<int,2>>
    auto& diags = out.first;
    // diagsinf[i] should be interpreted as vector<int>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      int dim = stptr_->dimension(s);
      auto v = stptr_->vertex_with_same_filtration(s);
      if (t == stptr_->null_simplex()) {
        while (static_cast<int>(diagsinf.size()) < dim + 1) diagsinf.emplace_back();
        diagsinf[dim].push_back(v);
      } else {
        while (static_cast<int>(diags.size()) < dim + 1) diags.emplace_back();
        auto w = stptr_->vertex_with_same_filtration(t);
        auto& d = diags[dim];
        d.insert(d.end(), {v, w});
      }
    }
    return out;
  }

  // An alternative, to avoid those different sizes, would be to "pad" vertex generator v as (v, v) or (v, -1). When
  // using it as index, this corresponds to adding the vertex filtration values either on the diagonal of the distance
  // matrix, or as an extra row or column. We could also merge the vectors for different dimensions into a single one,
  // with an extra column for the dimension (converted to type double).
  Generators flag_generators()
  {
    Generators out;
    // diags[0] should be interpreted as vector<array<int,3>> and other diags[i] as vector<array<int,4>>
    auto& diags = out.first;
    // diagsinf[0] should be interpreted as vector<int> and other diagsinf[i] as vector<array<int,2>>
    auto& diagsinf = out.second;
    for (auto pair : Base::get_persistent_pairs()) {
      auto s = std::get<0>(pair);
      auto t = std::get<1>(pair);
      int dim = stptr_->dimension(s);
      bool infinite = t == stptr_->null_simplex();
      if (infinite) {
        if (dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if (diagsinf.size() == 0) diagsinf.emplace_back();
          diagsinf[0].push_back(v);
        } else {
          auto e = stptr_->edge_with_same_filtration(s);
          auto&& e_vertices = stptr_->simplex_vertex_range(e);
          auto i = std::begin(e_vertices);
          auto v1 = *i;
          auto v2 = *++i;
          GUDHI_CHECK(++i == std::end(e_vertices), "must be an edge");
          while (diagsinf.size() < dim + 1) diagsinf.emplace_back();
          auto& d = diagsinf[dim];
          d.insert(d.end(), {v1, v2});
        }
      } else {
        auto et = stptr_->edge_with_same_filtration(t);
        auto&& et_vertices = stptr_->simplex_vertex_range(et);
        auto it = std::begin(et_vertices);
        auto w1 = *it;
        auto w2 = *++it;
        GUDHI_CHECK(++it == std::end(et_vertices), "must be an edge");
        if (dim == 0) {
          auto v = *std::begin(stptr_->simplex_vertex_range(s));
          if (diags.size() == 0) diags.emplace_back();
          auto& d = diags[0];
          d.insert(d.end(), {v, w1, w2});
        } else {
          auto es = stptr_->edge_with_same_filtration(s);
          auto&& es_vertices = stptr_->simplex_vertex_range(es);
          auto is = std::begin(es_vertices);
          auto v1 = *is;
          auto v2 = *++is;
          GUDHI_CHECK(++is == std::end(es_vertices), "must be an edge");
          while (diags.size() < dim + 1) diags.emplace_back();
          auto& d = diags[dim];
          d.insert(d.end(), {v1, v2, w1, w2});
        }
      }
    }
    return out;
  }

  using Filtration_value = typename FilteredComplex::Filtration_value;
  using Birth_death = std::pair<Filtration_value, Filtration_value>;
  using Persistence_subdiagrams = std::vector<std::vector<std::pair<int, Birth_death>>>;

  Persistence_subdiagrams compute_extended_persistence_subdiagrams(Filtration_value min_persistence)
  {
    Persistence_subdiagrams pers_subs(4);
    auto const& persistent_pairs = Base::get_persistent_pairs();
    for (auto pair : persistent_pairs) {
      std::pair<Filtration_value, Extended_simplex_type> px =
          stptr_->decode_extended_filtration(stptr_->filtration(get<0>(pair)), stptr_->efd);
      std::pair<Filtration_value, Extended_simplex_type> py =
          stptr_->decode_extended_filtration(stptr_->filtration(get<1>(pair)), stptr_->efd);
      std::pair<int, Birth_death> pd_point =
          std::make_pair(stptr_->dimension(get<0>(pair)), std::make_pair(px.first, py.first));
      if (std::abs(px.first - py.first) > min_persistence) {
        // Ordinary
        if (px.second == Extended_simplex_type::UP && py.second == Extended_simplex_type::UP) {
          pers_subs[0].push_back(pd_point);
        }
        // Relative
        else if (px.second == Extended_simplex_type::DOWN && py.second == Extended_simplex_type::DOWN) {
          pers_subs[1].push_back(pd_point);
        } else {
          // Extended+
          if (px.first < py.first) {
            pers_subs[2].push_back(pd_point);
          }
          // Extended-
          else {
            pers_subs[3].push_back(pd_point);
          }
        }
      }
    }
    return pers_subs;
  }

 private:
  // A copy
  FilteredComplex* stptr_;
};

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENT_COHOMOLOGY_INTERFACE_H_
