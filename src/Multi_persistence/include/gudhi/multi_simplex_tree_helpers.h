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
 * as well as the two helper methods @ref Gudhi::multi_persistence::make_multi_dimensional and
 * @ref Gudhi::multi_persistence::make_one_dimensional.
 */

#ifndef MP_MULTI_SIMPLEX_TREE_HELPERS_H_
#define MP_MULTI_SIMPLEX_TREE_HELPERS_H_

#include <cstddef>
#include <type_traits>

#include <gudhi/Debug_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Simplex_tree/simplex_tree_options.h>

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
    const OneDimSimplexTree &st,
    const typename MultiDimSimplexTreeOptions::Filtration_value &default_value,
    const std::size_t dimension = 0)
{
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
                                                            const std::size_t dimension = 0)
{
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

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_MULTI_SIMPLEX_TREE_HELPERS_H_
