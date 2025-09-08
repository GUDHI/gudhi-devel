/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file MultiFiltrationValue.h
 * @author Hannah Schreiber
 * @brief Contains the concept for multi-parameter filtration values.
 */

#include <initializer_list>
#include <vector>

namespace Gudhi {
namespace multi_persistence {

/**
 * @brief Concept for multi-parameter filtration values templated in various classes of the module:
 * @ref Multi_parameter_filtered_complex, @ref Projective_cover_kernel, @ref Slicer...
 * The concept is realized by @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration and
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 */
class MultiFiltrationValue
{
 public:
  using value_type = unspecified;           /**< Numerical value type. */
  using Underlying_container = unspecified; /**< Underlying container for values. */
  using size_type = unspecified;            /**< Size type. */
  using reference = unspecified;            /**< Reference type. */
  using const_reference = unspecified;      /**< Const reference type. */

  /**
   * @brief Default constructor. Builds filtration value with one generator and given number of parameters (if not
   * fixed by the class type).
   *
   * @param number_of_parameters If negative, takes the default value instead. Can be marked as `[[maybe_unused]]`
   * if the number of parameters is fixed for the class. Default value: 2.
   */
  MultiFiltrationValue(int number_of_parameters = 2);

  /**
   * @brief Builds filtration value with one generator and given number of parameters (if not fixed by the class type).
   * All values are initialized at the given value.
   *
   * @param number_of_parameters If negative, is set to 2 instead. Can be marked as `[[maybe_unused]]` if the number
   * of parameters is fixed for the class.
   * @param value Initialization value for every value in the generator.
   */
  MultiFiltrationValue(int number_of_parameters, value_type value);

  /**
   * @brief Builds filtration value with one generator that is initialized with the given range. The number of
   * parameters are therefore deduced from the length of the range. If the number of parameters is fixed for the
   * class type and the length is longer, then any exceeding value should just be ignored. If the length is shorter,
   * the behaviour is allowed to be undefined.
   *
   * @tparam ValueRange Range of types convertible to `value_type`. Can require a begin() and end() method.
   * @param range Values of the generator.
   */
  template <class ValueRange = std::initializer_list<value_type> >
  MultiFiltrationValue(const ValueRange &range);

  /**
   * @brief Builds filtration value with given number of parameters (can be ignored if the number of parameters
   * is fixed for the class) and values from the given range. Lets \f$ p \f$ be the number of parameters. The \f$ p \f$
   * first elements of the range have to correspond to the first generator, the \f$ p \f$ next elements to the second
   * generator and so on... So the length of the range has to be a multiple of \f$ p \f$ and the number of generators
   * will be \f$ length / p \f$. The range is represented by two iterators.
   *
   * @tparam Iterator Iterator type that can require satisfying the requirements of standard LegacyInputIterator and
   * dereferenced elements have to be convertible to `value_type`.
   * @param it_begin Iterator pointing to the start of the range.
   * @param it_end Iterator pointing to the end of the range.
   * @param number_of_parameters Negative values are associated to 0.
   */
  template <class Iterator>
  MultiFiltrationValue(Iterator it_begin, Iterator it_end, int number_of_parameters);

  /**
   * @brief Builds filtration value with given number of parameters (can be ignored if the number of parameters
   * is fixed for the class) and values from the given range. Let \f$ p \f$ be the number of parameters. The \f$ p \f$
   * first elements of the range have to correspond to the first generator, the \f$ p \f$ next elements to the second
   * generator and so on... So the length of the range has to be a multiple of \f$ p \f$ and the number of generators
   * will be \f$ length / p \f$. The range is represented by @ref MultiFiltrationValue::Underlying_container "" and
   * **moved** into the underlying container of the class.
   *
   * @param generators Values to move.
   * @param number_of_parameters Negative values are associated to 0.
   */
  MultiFiltrationValue(Underlying_container &&generators, int number_of_parameters);

  /**
   * @brief Copy constructor.
   */
  MultiFiltrationValue(const MultiFiltrationValue &other);

  /**
   * @brief Assign operator.
   */
  MultiFiltrationValue &operator=(const MultiFiltrationValue &other);

  /**
   * @brief Swap operator.
   */
  friend void swap(MultiFiltrationValue &f1, MultiFiltrationValue &f2) noexcept;

  /**
   * @brief Returns reference to value of parameter `p` of generator `g`.
   */
  reference operator()(size_type g, size_type p);

  /**
   * @brief Returns const reference to value of parameter `p` of generator `g`.
   */
  const_reference operator()(size_type g, size_type p) const;

  /**
   * @brief Returns a copy with entries casted into the type given as template parameter.
   *
   * @tparam U New type for the entries.
   * @return Copy with new entry type.
   */
  template <typename U>
  MultiFiltrationValue as_type() const;

  /**
   * @brief Returns the number of parameters in the filtration value.
   */
  size_type num_parameters() const;

  /**
   * @brief Returns the number of generators in the filtration value, i.e. the criticality of the element.
   */
  size_type num_generators() const;

  /**
   * @brief Returns a filtration value with given number of parameters equal to infinity (on all parameters).
   * The parameter @p number_of_parameters can be marked as `[[maybe_unused]]`.
   */
  static MultiFiltrationValue inf(int number_of_parameters);

  /**
   * @brief Returns a filtration value with given number of parameters equal to minus infinity (on all parameters).
   * The parameter @p number_of_parameters can be marked as `[[maybe_unused]]`.
   */
  static MultiFiltrationValue minus_inf(int number_of_parameters);

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically strictly less than the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator< "", this order
   * has to be total.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse>
  friend bool is_strict_less_than_lexicographically(const MultiFiltrationValue &a, const MultiFiltrationValue &b);

  /**
   * @brief Returns `true` if and only if the first argument is lexicographically less than or equal to the second
   * argument. The "words" considered for the lexicographical order are all the generators concatenated together
   * in order of generator index and then in order of parameter index. Different from @ref operator<= "", this order
   * has to be total.
   *
   * @tparam inverse If true, the parameter index and generator index order is inverted.
   */
  template <bool inverse>
  friend bool is_less_or_equal_than_lexicographically(const MultiFiltrationValue &a, const MultiFiltrationValue &b);

  /**
   * @brief Returns `true` if and only if the cones generated by @p b are strictly contained in the
   * cones generated by @p a (the cones can be positive or negative, but it has to be consistent within the same
   * class type). Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that this order does not have to be total (and in most cases will not).
   */
  friend bool operator<(const MultiFiltrationValue &a, const MultiFiltrationValue &b);

  /**
   * @brief Returns `true` if and only if the cones generated by @p a are strictly contained in the
   * cones generated by @p b (the cones can be positive or negative, but it has to be consistent within the same
   * class type). Both @p a and @p b  have to have the same number of parameters.
   *
   * Note that this order does not have to be total (and in most cases will not).
   */
  friend bool operator<=(const MultiFiltrationValue &a, const MultiFiltrationValue &b);

  /**
   * @brief Returns `true` if and only if the number of parameters and generators are the same and for each
   * \f$ i,j \f$, \f$ a(i,j) \f$ is equal to \f$ b(i,j) \f$.
   */
  friend bool operator==(const MultiFiltrationValue &a, const MultiFiltrationValue &b);

  /**
   * @brief Returns a filtration value such that an entry at index \f$ i,j \f$ is equal to \f$ -f(i,j) \f$.
   *
   * Used conventions:
   * - \f$ -NaN = NaN \f$.
   */
  friend MultiFiltrationValue operator-(const MultiFiltrationValue &f);

  /**
   * @brief Sets each generator of the filtration value to the least common upper bound between it and the given value.
   *
   * More formally, it pushes the current generator to the cone \f$ \{ y \in \mathbb R^n : y \ge x \} \f$
   * originating in \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \ge this \} \cap \{ y \in \mathbb R^n : y \ge x \} \f$.
   *
   * @tparam GeneratorRange Range of elements convertible to `value_type`. Can require a begin(), end() and size()
   * method.
   * @param x Range towards to push. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type> >
  bool push_to_least_common_upper_bound(const GeneratorRange &x, bool exclude_infinite_values = false);

  /**
   * @brief Sets each generator of the filtration value to the greatest common lower bound between it and the given
   * value.
   *
   * More formally, it pulls the current generator to the cone \f$ \{ y \in \mathbb R^n : y \le x \} \f$
   * originating in \f$ x \f$. The resulting value corresponds to the intersection of both
   * cones: \f$ \mathrm{this} = \min \{ y \in \mathbb R^n : y \le this \} \cap \{ y \in \mathbb R^n : y \le x \} \f$.
   *
   * @tparam GeneratorRange Range of elements convertible to `value_type`. Can require a begin(), end() and size()
   * method.
   * @param x Range towards to pull. Has to have as many elements than @ref num_parameters().
   * @param exclude_infinite_values If true, values at infinity or minus infinity are not affected.
   * @return true If the filtration value was actually modified.
   * @return false Otherwise.
   */
  template <class GeneratorRange = std::initializer_list<value_type> >
  bool pull_to_greatest_common_lower_bound(const GeneratorRange &x, bool exclude_infinite_values = false);

  /**
   * @brief Projects the filtration value into the given grid. If @p coordinate is false, the entries are set to
   * the nearest upper bound value with the same parameter in the grid. Otherwise, the entries are set to the indices
   * of those nearest upper bound values.
   * The grid has to be represented as a vector of ordered ranges of values convertible into `value_type`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OneDimArray A range of values convertible into `value_type` ordered by increasing value. Can require
   * a begin, end and operator[] method.
   * @param grid Vector of @p OneDimArray with size at least number of filtration parameters.
   * @param coordinate If true, the values are set to the coordinates of the projection in the grid. If false,
   * the values are set to the values at the coordinates of the projection.
   */
  template <typename OneDimArray>
  void project_onto_grid(const std::vector<OneDimArray> &grid, bool coordinate = true);

  /**
   * @brief Returns a generator with the minimal values of all parameters in any generator of the given filtration
   * value. That is, the greatest lower bound of all generators.
   */
  friend MultiFiltrationValue factorize_below(const MultiFiltrationValue &f);

  /**
   * @brief Returns a generator with the maximal values of all parameters in any generator of the given filtration
   * value. That is, the least upper bound of all generators.
   */
  friend MultiFiltrationValue factorize_above(const MultiFiltrationValue &f);

  /**
   * @brief Computes the coordinates in the given grid, corresponding to the nearest upper bounds of the entries
   * in the given filtration value.
   * The grid has to be represented as a vector of vectors of ordered values convertible into `OutValue`. An index
   * \f$ i \f$ of the vector corresponds to the same parameter as the index \f$ i \f$ in a generator of the filtration
   * value. The ranges correspond to the possible values of the parameters, ordered by increasing value, forming
   * therefore all together a 2D grid.
   *
   * @tparam OutValue Signed arithmetic type.
   * @tparam U Type which is convertible into `OutValue`.
   * @param f Filtration value to project.
   * @param grid Vector of vectors to project into.
   * @return Filtration value \f$ out \f$ whose entry correspond to the indices of the projected values. That is,
   * the projection of \f$ f(g,p) \f$ is \f$ grid[p][out(g,p)] \f$.
   */
  template <typename OutValue, typename U>
  friend MultiFiltrationValue compute_coordinates_in_grid(MultiFiltrationValue f,
                                                          const std::vector<std::vector<U> > &grid);

  /**
   * @brief Infinity value of an entry of the filtration value.
   */
  constexpr static const value_type T_inf = unspecified;
};

}  // namespace multi_persistence
}  // namespace Gudhi
