/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Nested_array_filtration.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_filtration::Nested_array_filtration class.
 */

#ifndef MF_NESTED_ARRAY_FILTRATION_H_
#define MF_NESTED_ARRAY_FILTRATION_H_

#include <cstddef>      //std::size_t
#include <limits>       //std::numerical_limits
#include <type_traits>  //std::is_arithmetic
#include <algorithm>    //std::sort
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/serialization_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

/**
 * @class Nested_array_filtration Nested_array_filtration.h gudhi/Multi_filtration/Nested_array_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value in a vector of vectors.
 * Implements the concept @ref StoragePolicy.
 *
 * @details For more documentation of the public interface, see @ref StoragePolicy.
 *
 * `std::numeric_limits<Multi_parameter_filtration_value<Nested_array_filtration>>` will behave such that:
 * - `::has_infinity` returns `true`,
 * - `::has_quiet_NaN` returns `true`,
 * - `::infinity(int)` returns `Nested_array_filtration::inf(size_type)`,
 * - `::minus_infinity(int)` returns `Nested_array_filtration::minus_inf(size_type)`,
 * - `::max(int)` returns a @ref Nested_array_filtration with one generator with all parameters at
 * std::numeric_limits<T>::max()`,
 * - `::quiet_NaN(int)` returns `Nested_array_filtration::nan(size_type)`.
 *
 * @tparam T Arithmetic type of an entry for one parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 */
template <typename T>
class Nested_array_filtration {
 private:
  using Generator = std::vector<T>;

 public:
  using value_type = T;                                       /**< `T` */
  using Underlying_container = std::vector<Generator>;        /**< std::vector<std::vector<T> > */
  using size_type = typename Underlying_container::size_type; /**< std::size_t */
  using reference = value_type &;                             /**< @ref value_type & */
  using const_reference = value_type;                         /**< @ref value_type */
  using iterator =
      typename Underlying_container::value_type::iterator;    /**< LegacyRandomAccessIterator Iterator type. */
  using const_iterator =
      typename Underlying_container::value_type::const_iterator; /**< LegacyRandomAccessIterator Const iterator type. */
  template <typename U>
  using As_type = Nested_array_filtration<U>;                 /**< Nested_array_filtration\<U\> */

  template <typename U = T>
  constexpr static const U T_inf = detail::MF_T_inf<U>;

  template <typename U = T>
  constexpr static const U T_m_inf = detail::MF_T_m_inf<U>;

  constexpr static const bool has_an_implicit_axis = false;
  constexpr static const bool has_lexicographical_storage = true;
  constexpr static const bool has_minimal_set_representation = true;
  template <bool Co>
  constexpr static const bool has_infinity = true;
  template <bool Co>
  constexpr static const bool has_minus_infinity = true;
  constexpr static const bool has_quiet_NaN = true;

  Nested_array_filtration() = default;

  Nested_array_filtration(size_type numberOfParameters, T value)
      : numberOfParameters_(numberOfParameters), generators_(1, Generator(numberOfParameters, value)) {}

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Nested_array_filtration(Iterator itBegin, Iterator itEnd) : generators_(1, Generator(itBegin, itEnd)) {
    numberOfParameters_ = generators_[0].size();
  }

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Nested_array_filtration(Iterator itBegin, Iterator itEnd, size_type numberOfParameters)
      : numberOfParameters_(numberOfParameters), generators_() {
    if (numberOfParameters_ == 0) return;

    // Will discard any value at the end which does not fit into a complete generator.
    const size_type num_gen = std::distance(itBegin, itEnd) / numberOfParameters_;

    Iterator it = itBegin;
    for (size_type i = 0; i < num_gen; ++i) {
      Iterator endIt = it;
      std::advance(endIt, numberOfParameters);
      generators_.emplace_back(it, endIt);
      it = endIt;
    }
  }

  Nested_array_filtration(const Underlying_container &generators, size_type numberOfParameters)
      : numberOfParameters_(numberOfParameters), generators_(generators) {
    GUDHI_CHECK(numberOfParameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));
  }

  Nested_array_filtration(Underlying_container &&generators, size_type numberOfParameters)
      : numberOfParameters_(numberOfParameters), generators_(std::move(generators)) {
    GUDHI_CHECK(numberOfParameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));
  }

  Nested_array_filtration(const Nested_array_filtration &other) = default;

  Nested_array_filtration(Nested_array_filtration &&other) noexcept = default;

  ~Nested_array_filtration() = default;

  Nested_array_filtration &operator=(const Nested_array_filtration &other) = default;

  Nested_array_filtration &operator=(Nested_array_filtration &&other) noexcept = default;

  friend void swap(Nested_array_filtration &f1, Nested_array_filtration &f2) noexcept {
    std::swap(f1.numberOfParameters_, f2.numberOfParameters_);
    f1.generators_.swap(f2.generators_);
  }

  template <typename U = T>
  As_type<U> copy(size_type numberOfParameters, size_type numberOfGenerators, U defaultValue) const {
    typename As_type<U>::Underlying_container newContainer(numberOfGenerators,
                                                           Generator(numberOfParameters, defaultValue));
    for (size_type g = 0; g < std::min(numberOfGenerators, num_generators()); ++g) {
      for (size_type p = 0; p < std::min(numberOfParameters, num_parameters()); ++p) {
        newContainer[g][p] = operator()(g, p);
      }
    }
    return As_type<U>(std::move(newContainer), numberOfParameters);
  }

  Underlying_container &get_underlying_container() { return generators_; }

  const Underlying_container &get_underlying_container() const { return generators_; }

  reference operator()(size_type g, size_type p) {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("First index out of bound."));
    GUDHI_CHECK(p < numberOfParameters_, std::out_of_range("Second index out of bound."));
    GUDHI_CHECK(!generators_[g].empty(), std::logic_error("Operator() cannot be called on NaN value."));
    if (generators_[g].size() < numberOfParameters_) generators_[g].resize(numberOfParameters_, generators_[g][0]);
    return generators_[g][p];
  }

  const_reference operator()(size_type g, size_type p) const {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("First index out of bound."));
    GUDHI_CHECK(p < numberOfParameters_, std::out_of_range("Second index out of bound."));
    GUDHI_CHECK(!generators_[g].empty(), std::logic_error("Operator() cannot be called on NaN value."));
    if (generators_[g].size() < numberOfParameters_) return generators_[g][0];
    return generators_[g][p];
  }

  iterator begin(size_type g) {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Generator index out of bound."));
    return generators_[g].begin();
  }

  iterator end(size_type g) {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Generator index out of bound."));
    return generators_[g].end();
  }

  const_iterator begin(size_type g) const {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Generator index out of bound."));
    return generators_[g].begin();
  }

  const_iterator end(size_type g) const {
    GUDHI_CHECK(g < generators_.size(), std::out_of_range("Generator index out of bound."));
    return generators_[g].end();
  }

  void reserve(size_type number_of_generators) { generators_.reserve(number_of_generators); }

  size_type num_parameters() const { return numberOfParameters_; }

  size_type num_generators() const { return generators_.size(); }

  size_type num_entries() const { return generators_.size() * numberOfParameters_; }

  template <bool Co>
  static Nested_array_filtration inf(size_type numberOfParameters) {
    if (numberOfParameters == 0) return Nested_array_filtration();
    Underlying_container out(1, Generator(1, T_inf<>));
    return Nested_array_filtration(std::move(out), numberOfParameters);
  }

  template <bool Co>
  static Nested_array_filtration minus_inf(size_type numberOfParameters) {
    if (numberOfParameters == 0) return Nested_array_filtration();
    Underlying_container out(1, Generator(1, T_m_inf<>));
    return Nested_array_filtration(std::move(out), numberOfParameters);
  }

  static Nested_array_filtration nan(size_type numberOfParameters) {
    if (numberOfParameters == 0) return Nested_array_filtration();
    Underlying_container out(1);
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      out[0].push_back(std::numeric_limits<T>::quiet_NaN());
    }
    return Nested_array_filtration(std::move(out), numberOfParameters);
  }

  template <bool Co>
  [[nodiscard]] bool is_plus_inf() const {
    for (const Generator &g : generators_) {
      if (!_is_plus_inf(g)) return false;
    }
    return !generators_.empty();
  }

  template <bool Co>
  [[nodiscard]] bool is_minus_inf() const {
    for (const Generator &g : generators_) {
      if (!_is_minus_inf(g)) return false;
    }
    return !generators_.empty();
  }

  [[nodiscard]] bool is_nan() const {
    if (generators_.empty()) return false;
    for (const Generator &g : generators_) {
      if (!_is_nan(g)) return false;
    }
    return true;
  }

  template <bool Co>
  [[nodiscard]] bool is_finite() const {
    bool isInf = true, isMinusInf = true, isNan = true;
    for (const Generator &g : generators_) {
      if (g.empty()) return false;
      for (const T &v : g) {
        if (v != T_inf<>) isInf = false;
        if (v != T_m_inf<>) isMinusInf = false;
        if (!detail::_is_nan(v)) isNan = false;
        if (!isInf && !isMinusInf && !isNan) return true;
      }
    }
    return false;
  }

  void set_num_generators(size_type g, value_type defaultValue) {
    generators_.resize(g, Generator(numberOfParameters_, defaultValue));
  }

  template <class Iterator>
  void emplace_back(Iterator startIt, Iterator endIt) {
    // TODO: add infinity test here to resize?
    generators_.emplace_back(startIt, endIt);
  }

  void sort() {
    // orders such that -inf < 'finite values'  < inf < NaN.
    std::sort(generators_.begin(), generators_.end(), [this](const Generator &g1, const Generator &g2) {
      bool g1IsNan = _is_nan(g1);
      if (g1IsNan || _is_nan(g2)) return !g1IsNan;
      if (_is_plus_inf(g1)) return false;
      if (_is_plus_inf(g2)) return true;
      if (_is_minus_inf(g2)) return false;
      if (_is_minus_inf(g1)) return true;

      // same size if reaching here
      for (size_type i = 0; i < g1.size(); ++i) {
        if (detail::_is_nan(g1[i]) && detail::_is_nan(g2[i])) continue;
        if (detail::_is_nan(g2[i])) return true;
        if (detail::_is_nan(g1[i])) return false;
        if (g1[i] != g2[i]) return g1[i] < g2[i];
      }
      return false;
    });
  }

  void swap_generators(size_type g1, size_type g2) {
    GUDHI_CHECK(g1 < num_generators() && g2 < num_generators(),
                std::out_of_range("Generator indices to swap are out of bound."));

    generators_[g1].swap(generators_[g2]);
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Nested_array_filtration &f) {
    stream << "( k = " << f.num_generators() << " ) ( p = " << f.num_parameters() << " ) [ ";
    for (const Generator &g : f.generators_) {
      stream << "[";
      for (value_type v : g) {
        stream << v << ", ";
      }
      stream << "]; ";
    }
    stream << " ]";

    return stream;
  }

  friend char *serialize_value_to_char_buffer(const Nested_array_filtration &value, char *start) {
    char *curr = start;
    curr = serialize_value_to_char_buffer(value.numberOfParameters_, curr);
    curr = serialize_value_to_char_buffer(value.generators_, curr);
    return curr;
  }

  friend const char *deserialize_value_from_char_buffer(Nested_array_filtration &value, const char *start) {
    const char *curr = start;
    curr = deserialize_value_from_char_buffer(value.numberOfParameters_, curr);
    curr = deserialize_value_from_char_buffer(value.generators_, curr);
    return curr;
  }

  friend std::size_t get_serialization_size_of(const Nested_array_filtration &value) {
    std::size_t size = get_serialization_size_of(value.numberOfParameters_);
    size += get_serialization_size_of(value.generators_);
    return size;
  }

 private:
  size_type numberOfParameters_;    /**< Number of parameters. */
  Underlying_container generators_; /**< Container of the filtration value elements. */

  [[nodiscard]] bool _is_plus_inf(const Generator &g) const {
    for (const T &v : g) {
      if (v != T_inf<>) return false;
    }
    return !g.empty();
  }

  [[nodiscard]] bool _is_minus_inf(const Generator &g) const {
    for (const T &v : g) {
      if (v != T_m_inf<>) return false;
    }
    return !g.empty();
  }

  [[nodiscard]] bool _is_nan(const Generator &g) const {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      for (const T &v : g) {
        if (!std::isnan(v)) return false;
      }
      return true;
    } else {
      return g.empty();
    }
  }
};

}  // namespace Gudhi::multi_filtration

#endif  // MF_NESTED_ARRAY_FILTRATION_H_
