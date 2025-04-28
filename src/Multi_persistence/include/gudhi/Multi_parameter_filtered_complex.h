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

#include <cstdint>  //std::uint32_t
#include <ostream>
#include <vector>

// #include <gudhi/Debug_utils.h>

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
  using Boundary_container = std::vector<std::vector<Index> >;
  using Dimension = int;

  // has to guarantee sorted by dimension.
  Multi_parameter_filtered_complex() {}

  Index get_number_of_cycle_generators() const;
  Index get_number_of_parameters() const;

  bool filtration_is_one_critical() const;

  const Filtration_value_container& get_filtration_values() const;
  Filtration_value_container& get_filtration_values();

  const std::vector<Dimension>& get_dimensions() const;
  // std::vector<Dimension>& get_dimensions();

  const Boundary_container& get_boundaries() const;
  // Boundary_container& get_boundaries();

  Dimension get_dimension(Index i) const;
  Dimension get_max_dimension() const;

  int prune_above_dimension(int maxDim);

  void coarsen_on_grid(const std::vector<std::vector<T>> &grid, bool coordinate = true);

  // Warning: be carefull that the permutation preserves ordered by dim
  friend Multi_parameter_filtered_complex build_permuted_complex(const Multi_parameter_filtered_complex& complex,
                                                                 const std::vector<Index>& permutation);

  friend std::pair<Multi_parameter_filtered_complex, std::vector<Index> > build_permuted_complex(
      const Multi_parameter_filtered_complex& complex);

  friend auto build_complex_coarsen_on_grid(const Multi_parameter_filtered_complex& complex,
                                            const std::vector<std::vector<T> >& grid)
  {
    // using return_filtration_value = decltype(std::declval<Filtration_value>().template as_type<std::int32_t>());
    // using return_complex = Multi_parameter_filtered_complex<return_filtration_value>;
    // typename return_complex::Filtration_value_container coords(slicer.get_number_of_cycle_generators());
    // for (std::size_t gen = 0u; gen < coords.size(); ++gen) {
    //   coords[gen] = compute_coordinates_in_grid<std::int32_t>(generator_filtration_values[gen], grid);
    // }
    // Multi_parameter_filtered_complex<return_filtration_value> outComplex();
  }

  friend std::ostream& operator<<(std::ostream& stream, const Multi_parameter_filtered_complex& complex);
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_FILTERED_COMPLEX_H_INCLUDED
