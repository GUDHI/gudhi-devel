/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/04 Hannah Schreiber: simplifications with new simplex tree constructors + name changes
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file multi_simplex_tree_helpers.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Simplex_tree_options_multidimensional_filtration struct,
 * as well as the helper methods @ref Gudhi::multi_persistence::make_multi_dimensional,
 * @ref Gudhi::multi_persistence::make_one_dimensional, @ref Gudhi::multi_persistence::build_simplex_tree_from_complex
 * and @ref Gudhi::multi_persistence::fill_axis_with_lowerstar.
 */

#ifndef MP_MULTI_SIMPLEX_TREE_HELPERS_H_
#define MP_MULTI_SIMPLEX_TREE_HELPERS_H_

#include <cstddef>
#include <type_traits>
#include <stdexcept>

#include <gudhi/Debug_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/simplex_tree_options.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Multi_persistence/Line.h>
#include <gudhi/Multi_persistence/utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @ingroup multi_persistence
 *
 * @brief Model of @ref SimplexTreeOptions. Same as @ref Gudhi::Simplex_tree_options_default but with a custom
 * filtration value type.
 *
 * @tparam MultiFiltrationValue Has to respect the @ref FiltrationValue concept.
 */
template <typename MultiFiltrationValue>
struct Simplex_tree_options_multidimensional_filtration : Simplex_tree_options_default {
  using Filtration_value = MultiFiltrationValue;
};

/**
 * @ingroup multi_persistence
 *
 * @brief Constructs a multi-dimensional simplex tree from the given one-dimensional simplex tree.
 *
 * All simplices are copied from the one-dimensional simplex tree \f$ st \f$ to the multi-dimensional simplex tree
 * \f$ st_multi \f$. To begin, all filtration values of \f$ st_multi \f$ are initialized to the given default value.
 * Then, all filtration values of \f$ st \f$ are projected onto \f$ st_multi \f$ at the given dimension index.
 *
 * @tparam MultiDimSimplexTreeOptions Options for the multi-dimensional simplex tree. Should follow the
 * @ref SimplexTreeOptions concept. It has to define a @ref FiltrationValue with the additional methods:
 * `num_parameters()` which returns the number of parameters, `num_generators()` which returns the number of generators
 * and `operator(g, p)` which return a (modifiable) reference to the \f$ p^{th} \f$ element of the \f$ g^{th} \f$
 * generator. It should also define a type `value_type` with the type of an element in the filtration value.
 * @tparam OneDimSimplexTree Type of the one-dimensional @ref Gudhi::Simplex_tree. The `Filtration_value` type has to
 * be convertible to the `Filtration_value::value_type` of `MultiDimSimplexTreeOptions`.
 * @param st Simplex tree to project.
 * @param default_value Default value of the multi-dimensional filtration values. Has therefore to contain at least
 * one generator and as many parameters than the final tree should have. One of the elements of the first generator
 * will take the value of the projected value, so make sure to initialize the default value such that the `operator()`
 * makes the change of value possible.
 * @param dimension Dimension index to which the filtration values should be projected onto.
 */
template <class MultiDimSimplexTreeOptions, class OneDimSimplexTree>
Simplex_tree<MultiDimSimplexTreeOptions> make_multi_dimensional(
    const OneDimSimplexTree &st, const typename MultiDimSimplexTreeOptions::Filtration_value &default_value,
    const std::size_t dimension = 0) {
  using OneDimF = typename OneDimSimplexTree::Options::Filtration_value;
  using MultiDimF = typename MultiDimSimplexTreeOptions::Filtration_value;

  static_assert(std::is_convertible_v<OneDimF, typename MultiDimF::value_type>,
                "A filtration value of the one dimensional tree should be convertible to an element of a filtration "
                "value of the multi dimensional simplex tree.");

  auto num_param = default_value.num_parameters();

  GUDHI_CHECK(dimension < num_param,
              "Given dimension is too high, it has to be smaller than the number of parameters.");
  GUDHI_CHECK(default_value.num_generators() > 0,
              "The default value for the filtration values should contain at least one generator.");

  auto translate = [&](const OneDimF &f) -> MultiDimF {
    auto res = default_value;
    res(0, dimension) = f;
    return res;
  };

  Simplex_tree<MultiDimSimplexTreeOptions> multi_st(st, translate);
  multi_st.set_num_parameters(num_param);

  return multi_st;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Constructs a one-dimensional simplex tree from the given multi-dimensional simplex tree.
 *
 * All simplices are copied from the multi-dimensional simplex tree \f$ st \f$ to the one-dimensional simplex tree
 * \f$ st_one \f$. All filtration values of \f$ st_one \f$ are initialized with the element value at given dimension
 * index of the first generator of the corresponding multi-dimensional filtration value in \f$ st \f$.
 *
 * @tparam OneDimSimplexTreeOptions Options for the one-dimensional simplex tree. Should follow the
 * @ref SimplexTreeOptions concept.
 * @tparam MultiDimSimplexTree Type of the multi-dimensional @ref Gudhi::Simplex_tree. It has to define a
 * @ref FiltrationValue with the additional methods: `num_parameters()` which returns the number of parameters,
 * `num_generators()` which returns the number of generators and `operator(g, p)` which returns the value of the
 * \f$ p^{th} \f$ element of the \f$ g^{th} \f$ generator. It should also define a type `value_type` with the type of
 * an element in the filtration value, which has to be convertible to `Filtration_value` of `OneDimSimplexTreeOptions`.
 * @param st Simplex tree to project.
 * @param dimension Dimension index in the first generator to project.
 */
template <class OneDimSimplexTreeOptions, class MultiDimSimplexTree>
Simplex_tree<OneDimSimplexTreeOptions> make_one_dimensional(const MultiDimSimplexTree &st,
                                                            const std::size_t dimension = 0) {
  using OneDimF = typename OneDimSimplexTreeOptions::Filtration_value;
  using MultiDimF = typename MultiDimSimplexTree::Options::Filtration_value;

  static_assert(std::is_convertible_v<typename MultiDimF::value_type, OneDimF>,
                "An element of a filtration value of the multi dimensional tree should be convertible to a filtration "
                "value of the one dimensional simplex tree.");

  auto translate = [dimension](const MultiDimF &f) -> OneDimF {
    GUDHI_CHECK(dimension < f.num_parameters(),
                "Given dimension is too high, it has to be smaller than the number of parameters.");
    GUDHI_CHECK(f.num_generators() > 0, "A filtration value of the multi tree should contain at least one generator.");
    return f(0, dimension);
  };

  Simplex_tree<OneDimSimplexTreeOptions> one_st(st, translate);
  one_st.set_num_parameters(1);

  return one_st;
}

// TODO: unit test
/**
 * @brief Constructs a one-dimensional simplex tree from the given multi-dimensional simplex tree.
 *
 * All simplices are copied from the multi-dimensional simplex tree \f$ st \f$ to the one-dimensional simplex tree
 * \f$ st_one \f$. All filtration values of \f$ st_one \f$ are initialized with the element value at the intersection
 * point of the given line and the cone spanned by the corresponding multi-dimensional filtration value in \f$ st \f$.
 * 
 * @tparam OneDimSimplexTreeOptions Options for the one-dimensional simplex tree. Should follow the
 * @ref SimplexTreeOptions concept.
 * @tparam MultiDimSimplexTree Type of the multi-dimensional @ref Gudhi::Simplex_tree. It has to define a
 * @ref FiltrationValue with the additional methods: `num_parameters()` which returns the number of parameters,
 * `num_generators()` which returns the number of generators and `operator(g, p)` which returns the value of the
 * \f$ p^{th} \f$ element of the \f$ g^{th} \f$ generator. It should also define a type `value_type` with the type of
 * an element in the filtration value, which has to be convertible to `Filtration_value` of `OneDimSimplexTreeOptions`.
 * @tparam U Template argument of the @ref Gudhi::multi_persistence::Line class.
 * @param st Simplex tree to project.
 * @param line Line with positive slope into which to project the multi parameter filtration onto to obtain
 * a 1-parameter filtration.
 * @param dimension Coordinate of the point resulting from the projection into the line to store as filtration value.
 */
template <class OneDimSimplexTreeOptions, class MultiDimSimplexTree,
          typename U = typename MultiDimSimplexTree::Filtration_value::value_type>
Simplex_tree<OneDimSimplexTreeOptions> make_one_dimensional(const MultiDimSimplexTree &st, const Line<U> line,
                                                            const std::size_t dimension = 0) {
  using OneDimF = typename OneDimSimplexTreeOptions::Filtration_value;
  using MultiDimF = typename MultiDimSimplexTree::Options::Filtration_value;

  static_assert(std::is_convertible_v<typename MultiDimF::value_type, OneDimF>,
                "An element of a filtration value of the multi dimensional tree should be convertible to a filtration "
                "value of the one dimensional simplex tree.");

  auto translate = [dimension, &line](const MultiDimF &f) -> OneDimF {
    GUDHI_CHECK(dimension < f.num_parameters(),
                "Given dimension is too high, it has to be smaller than the number of parameters.");
    GUDHI_CHECK(f.num_generators() > 0, "A filtration value of the multi tree should contain at least one generator.");
    return line[line.compute_forward_intersection(f)][dimension];
  };

  Simplex_tree<OneDimSimplexTreeOptions> one_st(st, translate);
  one_st.set_num_parameters(1);

  return one_st;
}

namespace detail {

template <class SimplexTreeOptions, typename Index, class MultiFiltrationValue, typename I, typename D, class F>
void _insert_simplices(const std::vector<std::vector<Index>> &simplices, int numParam,
                       const Multi_parameter_filtered_complex<MultiFiltrationValue, I, D> &cpx, F &&insert_simplex) {
  using ST = Simplex_tree<SimplexTreeOptions>;
  for (std::size_t i = 0; i < simplices.size(); ++i) {
    auto &f_i = cpx.get_filtration_values()[i];

    if constexpr (std::is_same_v<MultiFiltrationValue, typename ST::Filtration_value>) {
      if (numParam >= 0 && static_cast<std::size_t>(numParam) != f_i.num_parameters()) {
        auto f = f_i.copy(numParam, f_i.num_generators());
        insert_simplex(simplices[i], f);
      } else {
        insert_simplex(simplices[i], f_i);
      }
    } else {
      if (numParam >= 0 && static_cast<std::size_t>(numParam) != f_i.num_parameters()) {
        auto f = f_i.copy(numParam, f_i.num_generators()).template as_type<typename ST::Filtration_value>();
        insert_simplex(simplices[i], f);
      } else {
        auto f = f_i.template as_type<typename ST::Filtration_value>();
        insert_simplex(simplices[i], f);
      }
    }
  }
}

}  // namespace detail

// TODO: unit test
/**
 * @brief Constructs a multi-parameter filtered simplex tree from the given complex.
 * 
 * @tparam SimplexTreeOptions Options of the @ref Gudhi::Simplex_tree to construct.
 * @tparam MultiFiltrationValue First template argument of
 * @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex. Must be convertible into
 * @ref Gudhi::Simplex_tree::Filtration_value of the resulting simplex tree with a `as_type` method.
 * @tparam I Second template argument of @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex.
 * @tparam D Third template argument of @ref Gudhi::multi_persistence::Multi_parameter_filtered_complex.
 * @param cpx Complex to translate.
 * @param maxDim Maximal dimension to include in the translation. If negative, all dimensions are kept. Default: -1.
 * @param numParam New number of parameters if it should change. If negative, it remains the same. If the new value
 * is smaller than the old one, the generators are shortened from the end to fit. If the new value is greater,
 * the generators are extended from the end with -inf if `Co` and with +inf otherwise. Default: -1.
 */
template <class SimplexTreeOptions, class MultiFiltrationValue, typename I, typename D>
inline Simplex_tree<SimplexTreeOptions> build_simplex_tree_from_complex(
    const Multi_parameter_filtered_complex<MultiFiltrationValue, I, D> &cpx, int maxDim = -1, int numParam = -1) {
  // TODO: is_multi_filtration will discriminate all pre-made multi filtration classes, but not any user made
  // class following the MultiFiltrationValue concept (as it was more thought for inner use). The tests should be
  // re-thought or this one just removed.
  static_assert(multi_filtration::detail::RangeTraits<MultiFiltrationValue>::is_multi_filtration,
                "Target filtration value type has to correspond to the MultiFiltrationValue concept.");

  using ST = Simplex_tree<SimplexTreeOptions>;
  using Index = typename ST::Vertex_handle;

  const auto numberOfSimplices = cpx.get_number_of_cycle_generators();
  ST st;
  st.set_num_parameters(cpx.num_parameters());

  if (numberOfSimplices == 0) return st;

  if (cpx.is_ordered_by_dimension()) {
    std::vector<std::vector<Index>> simplices =
        get_vertices_from_ordered_boundaries<Index>(cpx.get_boundaries(), cpx.get_dimensions(), maxDim);
    detail::_insert_simplices<SimplexTreeOptions>(
        simplices, numParam, cpx,
        [&st](const std::vector<Index> &s, const typename ST::Filtration_value &f) { st.insert_simplex(s, f); });
  } else {
    std::vector<std::vector<Index>> simplices =
        get_vertices_from_boundaries<Index>(cpx.get_boundaries(), cpx.get_dimensions(), maxDim);
    detail::_insert_simplices<SimplexTreeOptions>(
        simplices, numParam, cpx, [&st](const std::vector<Index> &s, const typename ST::Filtration_value &f) {
          // if out of scope of maxDim, the simplex is empty
          if (!s.empty()) st.insert_simplex_and_subfaces(ST::Filtration_maintenance::IGNORE_VALIDITY, s, f);
        });
  }

  return st;
}

// TODO: unit test
/**
 * @brief Fills the values at given parameter of the first generator of all filtration values in the given simplex tree
 * with a lower star filtration generated from the given vertex filtration values.
 *
 * @tparam MultiDimSimplexTree Type of the multi-dimensional @ref Gudhi::Simplex_tree. It has to define a
 * @ref FiltrationValue with the additional methods: `num_parameters()` which returns the number of parameters,
 * and `operator(g, p)` which returns the value of the \f$ p^{th} \f$ element of the \f$ g^{th} \f$ generator.
 * It should also define a type `value_type` with the type of an element in the filtration value.
 * @tparam RandomAccessRange Random access range of value convertible into
 * `MultiDimSimplexTree::Filtration_value::value_type`.
 * @param st Simplex tree to modify.
 * @param vertexFiltration Initial vertex filtration values for the 1D lower star filtration.
 * @param axis Parameter to fill with the lower star filtration.
 */
template <class MultiDimSimplexTree, class RandomAccessRange>
void fill_axis_with_lowerstar(MultiDimSimplexTree &st, const RandomAccessRange &vertexFiltration, std::size_t axis) {
  using T = typename MultiDimSimplexTree::Filtration_value::value_type;
  for (auto sh : st.complex_simplex_range()) {
    auto &current_birth = st.get_filtration_value(sh);
    T maxValue = Gudhi::multi_filtration::detail::MF_T_m_inf<T>;
    for (auto vertex : st.simplex_vertex_range(sh)) {
      GUDHI_CHECK(static_cast<std::size_t>(vertex) < vertexFiltration.size(),
                  std::invalid_argument("Vertex filtration values does not have a value for every vertex."));
      GUDHI_CHECK(!Gudhi::multi_filtration::detail::_is_nan(vertexFiltration[vertex]),
                  std::invalid_argument("Filtration value should not be NaN."));
      maxValue = std::max(maxValue, static_cast<T>(vertexFiltration[vertex]));
    }
    GUDHI_CHECK(axis < current_birth.num_parameters(), std::invalid_argument("Axis is not a valid parameter index."));
    current_birth(0, axis) = maxValue;
  }
}

// TODO: unit test
/**
 * @brief Fills the values at given parameter of the first generator of all filtration values in the given simplex tree
 * with a lower star filtration generated from the given distance matrix.
 * 
 * @tparam MultiDimSimplexTree Type of the multi-dimensional @ref Gudhi::Simplex_tree. It has to define a
 * @ref FiltrationValue with the additional methods: `num_parameters()` which returns the number of parameters,
 * and `operator(g, p)` which returns the value of the \f$ p^{th} \f$ element of the \f$ g^{th} \f$ generator.
 * It should also define a type `value_type` with the type of an element in the filtration value.
 * @tparam RandomAccess2DRange 2-dimensional random access range of value convertible into
 * `MultiDimSimplexTree::Filtration_value::value_type`.
 * @param st Simplex tree to modify.
 * @param distanceMatrix Distance matrix. Must take the vertex handles of the simplex tree as indices (e.g. if the
 * simplex tree stores vertex 0, 1 and 3, there should be a distance entry at indices 0, 1 and 3. The index 2 will
 * be skipped) and be symmetric.
 * @param vertexValue Value for the vertices.
 * @param axis Parameter to fill with the lower star filtration.
 */
template <class MultiDimSimplexTree, class RandomAccess2DRange>
void fill_axis_with_distance_matrix(MultiDimSimplexTree &st, const RandomAccess2DRange &distanceMatrix,
                                    typename MultiDimSimplexTree::Filtration_value::value_type vertexValue,
                                    std::size_t axis) {
  using T = typename MultiDimSimplexTree::Filtration_value::value_type;
  for (auto sh : st.complex_simplex_range()) {
    auto &current_birth = st.get_filtration_value(sh);
    T maxValue = vertexValue;
    for (auto v1 : st.simplex_vertex_range(sh)) {
      for (auto v2 : st.simplex_vertex_range(sh)) {
        if (v1 < v2) {
          GUDHI_CHECK(static_cast<std::size_t>(v1) < distanceMatrix.size(),
                      std::invalid_argument("Distance matrix values does not have a column for every vertex."));
          GUDHI_CHECK(static_cast<std::size_t>(v2) < distanceMatrix[v1].size(),
                      std::invalid_argument("Distance matrix values does not have a row for every vertex."));
          auto val = distanceMatrix[v1][v2];
          GUDHI_CHECK(!Gudhi::multi_filtration::detail::_is_nan(val),
                      std::invalid_argument("Filtration value should not be NaN."));
          maxValue = std::max(maxValue, static_cast<T>(val));
        }
      }
    }
    GUDHI_CHECK(axis < current_birth.num_parameters(), std::invalid_argument("Axis is not a valid parameter index."));
    current_birth(0, axis) = maxValue;
  }
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MULTI_SIMPLEX_TREE_HELPERS_H_
