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
 * @file Flat_array_filtration.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_filtration::Flat_array_filtration class.
 */

#ifndef MF_FLAT_ARRAY_FILTRATION_H_
#define MF_FLAT_ARRAY_FILTRATION_H_

#include <algorithm>    //std::lower_bound
#include <cmath>        //std::isnan, std::min, std::abs
#include <cstddef>      //std::size_t
#include <cstring>      //memcpy
#include <ostream>      //std::ostream
#include <limits>       //std::numerical_limits
#include <type_traits>  //std::is_arithmetic
#include <utility>      //std::swap, std::move
#include <numeric>      //std::iota
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/serialization_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>

namespace Gudhi::multi_filtration {

/**
 * @class Flat_array_filtration Flat_array_filtration.h gudhi/Multi_filtration/Flat_array_filtration.h
 * @ingroup multi_filtration
 *
 * @brief Class encoding the different generators, i.e., apparition times, of a \f$ k \f$-critical
 * \f$\mathbb R^n\f$-filtration value in a single vector.
 * Implements the concept @ref StoragePolicy.
 *
 * @details For more documentation of the public interface, see @ref StoragePolicy.
 *
 * `std::numeric_limits<Multi_parameter_filtration_value<Flat_array_filtration>, Co>` will behave such that:
 * - `::has_infinity` returns `true`,
 * - `::has_quiet_NaN` returns `std::numeric_limits<T>::has_quiet_NaN`,
 * - `::infinity(int)` returns `Flat_array_filtration::inf(size_type)`,
 * - `::minus_infinity(int)` returns `Flat_array_filtration::minus_inf(size_type)`,
 * - `::max(int)` throws if `Co` is true and otherwise returns a @ref Flat_array_filtration with one generator with
 * all parameters at std::numeric_limits<T>::max()`,
 * - `::quiet_NaN(int)` returns `Flat_array_filtration::nan(size_type)` if `std::numeric_limits<T>::has_quiet_NaN`
 * and throws otherwise.
 * 
 * @tparam T Arithmetic type of an entry for one parameter of a filtration value. Has to be **signed** and
 * to implement `std::isnan(T)`, `std::numeric_limits<T>::has_quiet_NaN`, `std::numeric_limits<T>::quiet_NaN()`,
 * `std::numeric_limits<T>::has_infinity`, `std::numeric_limits<T>::infinity()` and `std::numeric_limits<T>::max()`.
 * If `std::numeric_limits<T>::has_infinity` returns `false`, a call to `std::numeric_limits<T>::infinity()`
 * can simply throw. Examples are the native types `double`, `float` and `int`.
 */
template <typename T>
class Flat_array_filtration {
 private:
  using view_extents = extents<std::size_t, Gudhi::dynamic_extent, Gudhi::dynamic_extent>;
  using Viewer = Gudhi::Simple_mdspan<T, view_extents>;

 public:
  using value_type = T;                                       /**< `T` */
  using Underlying_container = std::vector<T>;                /**< std::vector<T> */
  using size_type = typename Underlying_container::size_type; /**< std::size_t */
  using reference = value_type &;                             /**< @ref value_type & */
  using const_reference = const value_type &;                 /**< const @ref value_type & */
  using iterator = typename Underlying_container::iterator;   /**< LegacyRandomAccessIterator Iterator type. */
  using const_iterator =
      typename Underlying_container::const_iterator;          /**< LegacyRandomAccessIterator Const iterator type. */
  template <typename U>
  using As_type = Flat_array_filtration<U>;                   /**< Flat_array_filtration\<U\> */

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
  /**
   * @brief std::numeric_limits\<T\>::has_quiet_NaN
   */
  constexpr static const bool has_quiet_NaN = std::numeric_limits<T>::has_quiet_NaN;

  Flat_array_filtration() : generators_(), generator_view_(generators_.data(), 0, 0) {}

  Flat_array_filtration(size_type numberOfParameters, T value)
      : generators_(numberOfParameters, value),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size()) {}

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Flat_array_filtration(Iterator itBegin, Iterator itEnd)
      : generators_(itBegin, itEnd),
        generator_view_(generators_.data(), generators_.empty() ? 0 : 1, generators_.size()) {}

  template <class Iterator, class = std::enable_if_t<!std::is_arithmetic_v<Iterator>>>
  Flat_array_filtration(Iterator itBegin, Iterator itEnd, size_type numberOfParameters)
      : generators_(itBegin, itEnd),
        generator_view_(generators_.data(), numberOfParameters == 0 ? 0 : generators_.size() / numberOfParameters,
                        numberOfParameters) {
    GUDHI_CHECK(numberOfParameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));
  }

  Flat_array_filtration(const Underlying_container &generators, size_type numberOfParameters)
      : generators_(generators),
        generator_view_(generators_.data(), numberOfParameters == 0 ? 0 : generators_.size() / numberOfParameters,
                        numberOfParameters) {
    GUDHI_CHECK(numberOfParameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));
  }

  Flat_array_filtration(Underlying_container &&generators, size_type numberOfParameters)
      : generators_(std::move(generators)),
        generator_view_(generators_.data(), numberOfParameters == 0 ? 0 : generators_.size() / numberOfParameters,
                        numberOfParameters) {
    GUDHI_CHECK(numberOfParameters > 0 || generators_.empty(),
                std::invalid_argument("Number of parameters cannot be 0 if the container is not empty."));
  }

  Flat_array_filtration(const Flat_array_filtration &other)
      : generators_(other.generators_),
        generator_view_(generators_.data(), other.num_generators(), other.num_parameters()) {}

  Flat_array_filtration(Flat_array_filtration &&other) noexcept
      : generators_(std::move(other.generators_)),
        generator_view_(generators_.data(), other.num_generators(), other.num_parameters()) {}

  ~Flat_array_filtration() = default;

  Flat_array_filtration &operator=(const Flat_array_filtration &other) {
    generators_ = other.generators_;
    generator_view_ = Viewer(generators_.data(), other.num_generators(), other.num_parameters());
    return *this;
  }

  Flat_array_filtration &operator=(Flat_array_filtration &&other) noexcept {
    generators_ = std::move(other.generators_);
    generator_view_ = Viewer(generators_.data(), other.num_generators(), other.num_parameters());
    return *this;
  }

  friend void swap(Flat_array_filtration &f1, Flat_array_filtration &f2) noexcept {
    f1.generators_.swap(f2.generators_);
    swap(f1.generator_view_, f2.generator_view_);
  }

  template <typename U = T>
  As_type<U> copy(size_type numberOfParameters, size_type numberOfGenerators, U defaultValue) const {
    typename As_type<U>::Underlying_container newContainer(numberOfParameters * numberOfGenerators, defaultValue);
    size_type i = 0;
    for (size_type g = 0; g < numberOfGenerators; ++g) {
      for (size_type p = 0; p < numberOfParameters; ++p) {
        if (g < num_generators() && p < num_parameters()) newContainer[i] = operator()(g, p);
        ++i;
      }
    }
    return As_type<U>(std::move(newContainer), numberOfParameters);
  }

  Underlying_container &get_underlying_container() { return generators_; }

  const Underlying_container &get_underlying_container() const { return generators_; }

  reference operator()(size_type g, size_type p) { return generator_view_(g, p); }

  const_reference operator()(size_type g, size_type p) const { return generator_view_(g, p); }

  iterator begin(size_type g) {
    GUDHI_CHECK(g < num_generators(), std::out_of_range("Generator index out of bound"));
    return generators_.begin() + (g * num_parameters());
  }

  iterator end(size_type g) {
    GUDHI_CHECK(g < num_generators(), std::out_of_range("Generator index out of bound"));
    return generators_.begin() + ((g + 1) * num_parameters());
  }

  const_iterator begin(size_type g) const {
    GUDHI_CHECK(g < num_generators(), std::out_of_range("Generator index out of bound"));
    return generators_.begin() + (g * num_parameters());
  }

  const_iterator end(size_type g) const {
    GUDHI_CHECK(g < num_generators(), std::out_of_range("Generator index out of bound"));
    return generators_.begin() + ((g + 1) * num_parameters());
  }

  void reserve(size_type number_of_generators) {
    generators_.reserve(num_parameters() * number_of_generators);
    generator_view_.update_data(generators_.data());
  }

  size_type num_parameters() const { return generator_view_.extent(1); }

  size_type num_generators() const { return generator_view_.extent(0); }

  size_type num_entries() const { return generators_.size(); }

  template <bool Co>
  static Flat_array_filtration inf(size_type numberOfParameters) {
    return Flat_array_filtration(numberOfParameters, T_inf<>);
  }

  template <bool Co>
  static Flat_array_filtration minus_inf(size_type numberOfParameters) {
    return Flat_array_filtration(numberOfParameters, T_m_inf<>);
  }

  /**
   * @throw When @ref has_quiet_NaN is false.
   */
  static Flat_array_filtration nan(size_type numberOfParameters) {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      return Flat_array_filtration(numberOfParameters, std::numeric_limits<T>::quiet_NaN());
    } else {
      throw std::logic_error("No NaN value exists.");
    }
  }

  template <bool Co>
  [[nodiscard]] bool is_plus_inf() const {
    for (const T &v : generators_) {
      if (v != T_inf<>) return false;
    }
    return !generators_.empty();
  }

  template <bool Co>
  [[nodiscard]] bool is_minus_inf() const {
    for (const T &v : generators_) {
      if (v != T_m_inf<>) return false;
    }
    return !generators_.empty();
  }

  [[nodiscard]] bool is_nan() const {
    if constexpr (std::numeric_limits<T>::has_quiet_NaN) {
      for (const auto &v : generators_) {
        if (!std::isnan(v)) return false;
      }
      return true;
    } else {
      return false;
    }
  }

  template <bool Co>
  [[nodiscard]] bool is_finite() const {
    if (generators_.empty()) return false;

    bool isInf = true, isMinusInf = true, isNan = true;
    for (const auto &v : generators_) {
      if (v != T_inf<>) isInf = false;
      if (v != T_m_inf<>) isMinusInf = false;
      if (!detail::_is_nan(v)) isNan = false;
      if (!isInf && !isMinusInf && !isNan) return true;
    }
    return false;
  }

  void set_num_generators(size_type g, value_type defaultValue) {
    generators_.resize(g * num_parameters(), defaultValue);
    generator_view_.update_data(generators_.data());  // in case it was relocated
    generator_view_.update_extent(0, g);
  }

  template <class Iterator>
  void emplace_back(Iterator startIt, Iterator endIt) {
    generators_.insert(generators_.end(), startIt, endIt);
    generator_view_.update_data(generators_.data());  // in case it was relocated
    generator_view_.update_extent(0, num_generators() + 1);
  }

  void sort() {
    auto numGen = num_generators();
    std::vector<size_type> newToOld(numGen);
    std::iota(newToOld.begin(), newToOld.end(), 0);
    // orders such that -inf < 'finite values'  < inf < NaN.
    std::sort(newToOld.begin(), newToOld.end(), [this](size_type g1, size_type g2) {
      for (auto it1 = begin(g1), it2 = begin(g2); it1 != end(g1); ++it1, ++it2) {
        value_type v1 = *it1;
        value_type v2 = *it2;
        if (detail::_is_nan(v1) && detail::_is_nan(v2)) continue;
        if (detail::_is_nan(v2)) return true;
        if (detail::_is_nan(v1)) return false;
        if (v1 != v2) return v1 < v2;
      }
      return false;
    });
    std::vector<size_type> oldToNew(numGen);
    for (size_type i = 0; i < numGen; ++i) oldToNew[newToOld[i]] = i;

    for (size_type i = 0; i < num_generators(); ++i) {
      size_type j = oldToNew[i];
      if (j != i) {
        swap_generators(i, j);
        std::swap(oldToNew[i], oldToNew[j]);
      }
    }
  }

  void swap_generators(size_type g1, size_type g2) {
    GUDHI_CHECK(g1 < num_generators() && g2 < num_generators(),
                std::out_of_range("Generator indices to swap are out of bound."));

    if (g1 == g2) return;

    for (auto it1 = begin(g1), it2 = begin(g2); it1 != end(g1); ++it1, ++it2) {
      if (*it1 != *it2) std::swap(*it1, *it2);
    }
  }

  friend std::ostream &operator<<(std::ostream &stream, const Flat_array_filtration &f) {
    const size_type num_gen = f.num_generators();
    const size_type num_param = f.num_parameters();

    stream << "( k = " << num_gen << " ) ( p = " << num_param << " ) [ ";
    for (size_type i = 0; i < f.generators_.size(); ++i) {
      stream << f.generators_[i] << " ";
    }
    stream << "]";

    return stream;
  }

  friend char *serialize_value_to_char_buffer(const Flat_array_filtration &value, char *start) {
    char *curr = start;
    curr = serialize_value_to_char_buffer(value.num_generators(), curr);
    curr = serialize_value_to_char_buffer(value.num_parameters(), curr);
    curr = serialize_value_to_char_buffer(value.generators_, curr);
    return curr;
  }

  friend const char *deserialize_value_from_char_buffer(Flat_array_filtration &value, const char *start) {
    const char *curr = start;
    size_type numParam, numGen;
    curr = deserialize_value_from_char_buffer(numGen, curr);
    curr = deserialize_value_from_char_buffer(numParam, curr);
    curr = deserialize_value_from_char_buffer(value.generators_, curr);
    value.generator_view_ = Viewer(value.generators_.data(), numGen, numParam);
    return curr;
  }

  friend std::size_t get_serialization_size_of(const Flat_array_filtration &value) {
    std::size_t size = sizeof(size_type) * 2;
    size += get_serialization_size_of(value.generators_);
    return size;
  }

 private:
  Underlying_container generators_; /**< Container of the filtration value elements. */
  Viewer generator_view_;           /**< Matrix view of the container. Has to be created after generators_. */
};

}  // namespace Gudhi::multi_filtration

#endif  // MF_FLAT_ARRAY_FILTRATION_H_
