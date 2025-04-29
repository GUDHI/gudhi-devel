/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: Reorganization.
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Multi_parameter_filtered_complex.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex class.
 */

#ifndef MP_FILTERED_COMPLEX_H_INCLUDED
#define MP_FILTERED_COMPLEX_H_INCLUDED

#include <algorithm>
#include <cstdint>  //std::uint32_t
#include <numeric>
#include <ostream>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_parameter_filtration.h>  //for lex order

namespace Gudhi {
namespace multi_persistence {

// TODO: better name
template <class Multi_filtration_value>
class Multi_parameter_filtered_complex
{
 public:
  using Index = std::uint32_t;
  using Filtration_value = Multi_filtration_value;
  using T = typename Filtration_value::value_type;
  using Filtration_value_container = std::vector<Filtration_value>;
  using Boundary = std::vector<Index>;
  using Boundary_container = std::vector<Boundary>;
  using Dimension = int;
  using Dimension_container = std::vector<Dimension>;

  Multi_parameter_filtered_complex()
      : boundaries_(), dimensions_(), filtrationValues_(), maxDimension_(0), isOrderedByDimension_(true)
  {}

  // assumes boundary Idxs corresponds to container Idxs
  Multi_parameter_filtered_complex(const Boundary_container& boundaries,
                                   const Dimension_container& dimensions,
                                   const Filtration_value_container& filtrationValues)
      : boundaries_(boundaries), dimensions_(dimensions), filtrationValues_(filtrationValues), maxDimension_(0)
  {
    _initialize_dimension_utils();
  }

  // assumes boundary Idxs corresponds to container Idxs
  Multi_parameter_filtered_complex(Boundary_container&& boundaries,
                                   Dimension_container&& dimensions,
                                   Filtration_value_container&& filtrationValues)
      : boundaries_(std::move(boundaries)),
        dimensions_(std::move(dimensions)),
        filtrationValues_(std::move(filtrationValues)),
        maxDimension_(0)
  {
    _initialize_dimension_utils();
  }

  Index get_number_of_cycle_generators() const { return boundaries_.size(); }

  Index get_number_of_parameters() const
  {
    if (filtrationValues_.empty()) return 0;
    filtrationValues_[0].num_parameters();
  }

  bool is_ordered_by_dimension() const { return isOrderedByDimension_; }

  const Filtration_value_container& get_filtration_values() const { return filtrationValues_; }

  Filtration_value_container& get_filtration_values() { return filtrationValues_; }

  const Dimension_container& get_dimensions() const { return dimensions_; }

  Dimension_container& get_dimensions() { return dimensions_; }

  const Boundary_container& get_boundaries() const { return boundaries_; }

  Boundary_container& get_boundaries() { return boundaries_; }

  Dimension get_dimension(Index i) const { return dimensions_[i]; }

  Dimension get_max_dimension() const { return maxDimension_; }

  void sort_by_dimension_co_lexicographically()
  {
    sort([&](Index i, Index j) -> bool {
      if (dimensions_[i] == dimensions_[j]) {
        return multi_filtration::is_strict_less_than_lexicographically<true>(filtrationValues_[i],
                                                                             filtrationValues_[j]);
      }
      return dimensions_[i] < dimensions_[j];
    });

    isOrderedByDimension_ = true;
  }

  template <typename Comp>
  void sort(Comp&& comparaison)
  {
    std::vector<Index> perm(boundaries_.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::vector<Index> pos = perm;
    std::sort(perm.begin(), perm.end(), comparaison);
    std::vector<Index> inv(boundaries_.size());
    for (Index i = 0; i < perm.size(); ++i) inv[perm[i]] = i;

    for (Index curr = 0; curr < perm.size(); ++curr) {
      Index p = perm[curr];
      Index i = pos[p];
      if (i != curr) {
        GUDHI_CHECK(curr < i, "Something is wrong");
        std::swap(boundaries_[curr], boundaries_[i]);
        for (Index& b : boundaries_[curr]) b = inv[b];
        std::swap(dimensions_[curr], dimensions_[i]);
        swap(filtrationValues_[curr], filtrationValues_[i]);
        std::swap(pos[curr], pos[i]);
      }
    }
  }

  // warning: shuffles order if not ordered by dimension
  int prune_above_dimension(int maxDim)
  {
    if (!isOrderedByDimension_) sort_by_dimension_co_lexicographically();
    Index i = 0;
    while (i < dimensions_.size() && dimensions_[i] < maxDim + 1) ++i;
    boundaries_.resize(i);
    dimensions_.resize(i);
    filtrationValues_.resize(i);
    maxDimension_ = dimensions_.empty() ? -1 : dimensions_.back();
  }

  void coarsen_on_grid(const std::vector<std::vector<T> >& grid, bool coordinate = true)
  {
    for (auto gen = 0u; gen < filtrationValues_.size(); ++gen) {
      filtrationValues_[gen].project_onto_grid(grid, coordinate);
    }
  }

  friend Multi_parameter_filtered_complex build_permuted_complex(const Multi_parameter_filtered_complex& complex,
                                                                 const std::vector<Index>& permutation)
  {
    if (permutation.size() != complex.get_number_of_cycle_generators())
      throw std::invalid_argument("Invalid permutation size.");

    std::vector<Index> inv(permutation.size());
    for (Index i = 0; i < permutation.size(); ++i) inv[permutation[i]] = i;

    Boundary_container newBoundaries;
    newBoundaries.reserve(permutation.size());
    Dimension_container newDimensions;
    newDimensions.reserve(permutation.size());
    Filtration_value_container newFiltrationValues;
    newBoundaries.reserve(permutation.size());

    for (Index i : permutation) {
      Boundary boundary(complex.boundaries_[i]);
      for (Index& b : boundary) b = inv[b];
      newBoundaries.emplace_back(std::move(boundary));
      newDimensions.push_back(complex.dimensions_[i]);
      newFiltrationValues.emplace_back(complex.filtrationValues_[i]);
    }

    return Multi_parameter_filtered_complex(
        std::move(newBoundaries), std::move(newDimensions), std::move(newFiltrationValues));
  }

  friend std::pair<Multi_parameter_filtered_complex, std::vector<Index> > build_permuted_complex(
      const Multi_parameter_filtered_complex& complex)
  {
    std::vector<Index> perm(complex.get_number_of_cycle_generators());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](Index i, Index j) -> bool {
      if (complex.dimensions_[i] == complex.dimensions_[j]) {
        return multi_filtration::is_strict_less_than_lexicographically<true>(complex.filtrationValues_[i],
                                                                             complex.filtrationValues_[j]);
      }
      return complex.dimensions_[i] < complex.dimensions_[j];
    });
    return build_permuted_complex(complex, perm);
  }

  friend auto build_complex_coarsen_on_grid(const Multi_parameter_filtered_complex& complex,
                                            const std::vector<std::vector<T> >& grid)
  {
    using Return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    using Return_complex = Multi_parameter_filtered_complex<Return_filtration_value>;
    typename Return_complex::Filtration_value_container coords(complex.get_number_of_cycle_generators());
    for (std::size_t gen = 0u; gen < coords.size(); ++gen) {
      coords[gen] = multi_filtration::compute_coordinates_in_grid<std::int32_t>(complex.filtrationValues_[gen], grid);
    }
    return Return_complex(complex.boundaries_, complex.dimensions_, coords);
  }

  friend std::ostream& operator<<(std::ostream& stream, const Multi_parameter_filtered_complex& complex)
  {
    stream << "Boundary:\n";
    stream << "{\n";
    for (Index i = 0; i < complex.boundaries_.size(); ++i) {
      const auto& boundary = complex.boundaries_[i];
      stream << i << ": {";
      for (auto b : boundary) stream << b << ", ";
      if (!boundary.empty()) stream << "\b" << "\b ";
      stream << "},\n";
    }
    stream << "}\n";

    stream << "Dimensions: (max " << complex.get_max_dimension() << ")\n";
    stream << "{";
    for (auto d : complex.dimensions_) stream << d << ", ";
    if (!complex.dimensions_.empty()) {
      stream << "\b" << "\b";
    }
    stream << "}\n";

    stream << "Filtration values:\n";
    stream << "{";
    for (auto f : complex.filtrationValues_) stream << f << "\n";
    stream << "}\n";

    return stream;
  }

 private:
  Boundary_container boundaries_;
  Dimension_container dimensions_;
  Filtration_value_container filtrationValues_;
  Dimension maxDimension_;
  bool isOrderedByDimension_;

  void _initialize_dimension_utils()
  {
    isOrderedByDimension_ = true;
    for (Index i = 0; i < dimensions_.size() - 1; ++i) {
      maxDimension_ = std::max(dimensions_[i], maxDimension_);
      if (dimensions_[i] > dimensions_[i + 1]) isOrderedByDimension_ = false;
    }
    maxDimension_ = std::max(dimensions_.back(), maxDimension_);
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_FILTERED_COMPLEX_H_INCLUDED
