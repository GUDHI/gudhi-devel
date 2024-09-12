/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which
 * is released under MIT. See file LICENSE or go to
 * https://gudhi.inria.fr/licensing/ for full license details. Author(s): David
 * Loiseaux, Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#ifndef SIMPLEX_TREE_MULTI_H_
#define SIMPLEX_TREE_MULTI_H_

#include <algorithm>
#include <gudhi/Simplex_tree.h>

namespace Gudhi::multi_persistence {

/** Model of SimplexTreeOptions, with a multiparameter filtration.
 * \ingroup multi_persistence
 * */
template <typename Filtration>
struct Simplex_tree_options_multidimensional_filtration {
public:
  typedef linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  using Filtration_value = Filtration;
  typedef typename Filtration::value_type value_type;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
  static const bool is_multi_parameter = true;
};

/**
 * \brief Turns a 1-parameter simplextree into a multiparameter simplextree,
 * and keeps the 1-filtration in the 1st axis.
 * Default values can be specified.
 * \ingroup multi_persistence
 * \tparam simplextree_std A non-multi simplextree
 * \tparam simplextree_multi A multi simplextree
 * \param st Simplextree to copy
 * \param st_multi Multiparameter simplextree container to fill.
 * \param default_values If given, this vector is assume to be of size `num_parameters-1` and contains the default
 * values of axes `1` to `num_parameters`.
 * */
template <class simplextree_std, class simplextree_multi>
void multify(simplextree_std &st, simplextree_multi &st_multi,
             const int num_parameters,
             const typename simplextree_multi::Options::Filtration_value
                 &default_values = {}) {
  typename simplextree_multi::Options::Filtration_value f(num_parameters);
  static_assert(
      !simplextree_std::Options::is_multi_parameter &&
          simplextree_multi::Options::is_multi_parameter,
      "Can only convert non-multiparameter to multiparameter simplextree.");
  unsigned int num_default_values;
  if constexpr (simplextree_multi::Options::Filtration_value::
                    is_multi_critical) {
    num_default_values = default_values[0].size();
  } else {
    num_default_values = default_values.size();
  }
  for (auto i = 0u; i < std::min(num_default_values,
                                 static_cast<unsigned int>(num_parameters - 1));
       i++)
    if constexpr (simplextree_multi::Options::Filtration_value::
                      is_multi_critical) {
      f[0][i + 1] = default_values[0][i];
    } else {
      f[i + 1] = default_values[i];
    }

  std::vector<int> simplex;
  simplex.reserve(st.dimension() + 1);
  for (auto &simplex_handle : st.complex_simplex_range()) {
    simplex.clear();
    for (auto vertex : st.simplex_vertex_range(simplex_handle))
      simplex.push_back(vertex);

    if (num_parameters > 0) {
      if constexpr (simplextree_multi::Options::Filtration_value::
                        is_multi_critical) {
        f[0][0] = st.filtration(simplex_handle);
      } else {
        f[0] = st.filtration(simplex_handle);
      }
    }
    st_multi.insert_simplex(simplex, f);
  }
  st_multi.set_number_of_parameters(num_parameters);
}

/**
 * \brief Turns a multiparameter-parameter simplextree into a 1-parameter
 * simplextree.
 * \ingroup multi_persistence
 * \tparam simplextree_std A non-multi simplextree
 * \tparam simplextree_multi A multi simplextree
 * \param st Simplextree to fill.
 * \param st_multi Multiparameter simplextree to convert into a 1 parameter simplex tree.
 * \param dimension The filtration parameter to put into the 1 parameter simplextree.
 * */
template <class simplextree_std, class simplextree_multi>
void flatten(simplextree_std &st, simplextree_multi &st_multi,
             const int dimension = 0) {
  static_assert(
      !simplextree_std::Options::is_multi_parameter &&
          simplextree_multi::Options::is_multi_parameter,
      "Can only convert multiparameter to non-multiparameter simplextree.");
  for (const auto &simplex_handle : st_multi.complex_simplex_range()) {
    std::vector<int> simplex;
    typename simplextree_multi::Options::value_type f;
    for (auto vertex : st_multi.simplex_vertex_range(simplex_handle))
      simplex.push_back(vertex);
    if constexpr (simplextree_multi::Filtration_value::is_multi_critical) {
      f = dimension >= 0 ? st_multi.filtration(simplex_handle)[0][dimension]
                         : 0;
    } else {
      f = dimension >= 0 ? st_multi.filtration(simplex_handle)[dimension] : 0;
    }
    st.insert_simplex(simplex, f);
  }
}

} // namespace Gudhi::multi_persistence

#endif // SIMPLEX_TREE_MULTI_H_
