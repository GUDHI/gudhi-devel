/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2026/02 Hannah Schreiber: reorganization + small optimizations + documentation
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Summand.h
 * @author David Loiseaux
 * @brief Contains the @ref Gudhi::multi_persistence::Summand class.
 */

#ifndef MP_SUMMAND_H_INCLUDED
#define MP_SUMMAND_H_INCLUDED

#include <cstddef>    //std::size_t
#include <ostream>    //std::ostream
#include <stdexcept>  //std::invalid_argument, std::runtime_error
#include <utility>
#include <vector>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/parallel_for.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Dynamic_multi_parameter_filtration.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Line.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Summand Summand.h gudhi/Multi_persistence/Summand.h
 * @ingroup multi_persistence
 *
 * @brief Summand of a multi-parameter persistence module, represented by a range of "birth corners" and a range of
 * "death corners".
 *
 * @tparam T Type of a parameter in the filtration value.
 * @tparam D Dimension type. Default: int.
 */
template <typename T, typename D = int>
class Summand {
 public:
  using value_type = T;                         /**< Filtration value parameter type. */
  using Dimension = D;                          /**< Dimension type. */
  /**
   * @brief Birth corners range type.
   */
  using Births = Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<value_type, false, false>;
  /**
   * @brief Death corners range type.
   */
  using Deaths = Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<value_type, true, false>;
  using Index = typename Births::size_type;     /**< Range index type. */

  static constexpr T T_inf = Births::T_inf;     /**< Parameter infinite value. */
  static constexpr T T_m_inf = Births::T_m_inf; /**< Parameter minus infinite value. */

  /**
   * @brief Default constructor. Constructs an empty summand.
   *
   * @param numberOfParameters Number of parameters of the filtration values. Default: 1.
   */
  Summand(int numberOfParameters = 1)
      : birthCorners_(numberOfParameters, T_inf),
        deathCorners_(numberOfParameters, T_m_inf),
        dimension_(get_null_value<Dimension>()) {}

  /**
   * @brief Constructs a summand from the given corners.
   *
   * @param birthCorners Birth corners.
   * @param deathCorners Death corners.
   * @param dimension Dimension of the summand.
   */
  Summand(const Births &birthCorners, const Deaths &deathCorners, Dimension dimension)
      : birthCorners_(birthCorners), deathCorners_(deathCorners), dimension_(dimension) {}

  /**
   * @brief Constructs a summand from the given corners. The corners are represented by ranges of coordinates such that
   * the \f$ p \f$ first elements of the range corresponds to the first corner, the \f$ p \f$ next elements to the
   * second corner and so on... Where \f$ p \f$ is the number of parameters.
   *
   * @tparam ValueRange Range of type convertible into `T` with a begin() and end() method.
   * @param birthCorners Birth corners.
   * @param deathCorners Death corners.
   * @param numberOfParameters Number of parameters.
   * @param dimension Dimension of the summand.
   */
  template <class ValueRange = std::initializer_list<value_type>>
  Summand(const ValueRange &birthCorners, const ValueRange &deathCorners, int numberOfParameters, Dimension dimension)
      : birthCorners_(birthCorners.begin(), birthCorners.end(), numberOfParameters),
        deathCorners_(deathCorners.begin(), deathCorners.end(), numberOfParameters),
        dimension_(dimension) {}

  /**
   * @brief Returns the dimension of the summand.
   */
  Dimension get_dimension() const { return dimension_; }

  /**
   * @brief Sets the dimension of summand to the given dimension.
   */
  void set_dimension(Dimension dimension) { dimension_ = dimension; }

  /**
   * @brief Returns the number of parameters.
   */
  [[nodiscard]] int get_number_of_parameters() const { return birthCorners_.num_parameters(); }

  /**
   * @brief Returns the number of birth corners (positive generators) in the summand.
   */
  [[nodiscard]] int get_number_of_birth_corners() const { return birthCorners_.num_generators(); }

  /**
   * @brief Returns the number of death corners (negative generators) in the summand.
   */
  [[nodiscard]] int get_number_of_death_corners() const { return deathCorners_.num_generators(); }

  /**
   * @brief Returns `true` if and only if the given filtration value is contained in summand.
   *
   * @tparam MultiFiltrationValue Either @ref Gudhi::multi_filtration::Multi_parameter_filtration,
   * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration or
   * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
   * @param x Multi-parameter filtration value.
   */
  template <class MultiFiltrationValue>
  bool contains(const MultiFiltrationValue &x) const {
    auto xPos = Gudhi::multi_filtration::as_type<Births>(x);
    auto xNeg = Gudhi::multi_filtration::as_type<Deaths>(x);
    return birthCorners_ <= xPos && deathCorners_ <= xNeg;
  }

  /**
   * @brief Returns the minimal box containing the summand.
   */
  Box<value_type> compute_bounds() const {
    if (birthCorners_.num_generators() == 0) return {};

    auto numParam = birthCorners_.num_parameters();
    typename Box<value_type>::Point_t m(numParam, T_inf);
    typename Box<value_type>::Point_t M(numParam, T_m_inf);

    for (Index g = 0; g < birthCorners_.num_generators(); ++g) {
      for (Index p = 0; p < numParam; ++p) {
        m[p] = std::min(m[p], birthCorners_(g, p));
      }
    }
    for (Index g = 0; g < deathCorners_.num_generators(); ++g) {
      for (Index p = 0; p < numParam; ++p) {
        auto corner_i = deathCorners_(g, p);
        if (corner_i != T_inf) M[p] = std::max(M[p], corner_i);
      }
    }

    return Box<value_type>(m, M);
  }

  /**
   * @brief Returns const reference to the range of birth corners.
   */
  const Births &get_upset() const { return birthCorners_; }

  /**
   * @brief Returns const reference to the range of death corners.
   */
  const Deaths &get_downset() const { return deathCorners_; }

  /**
   * @brief Flatten the range of birth corners as vector of `T`, such that the \f$ p \f$ first elements of the range
   * corresponds to the first corner, the \f$ p \f$ next elements to the second corner and so on... Where \f$ p \f$ is
   * the number of parameters.
   */
  std::vector<value_type> compute_flat_upset() const { return _compute_flat_set(birthCorners_); }

  /**
   * @brief Flatten the range of death corners as vector of `T`, such that the \f$ p \f$ first elements of the range
   * corresponds to the first corner, the \f$ p \f$ next elements to the second corner and so on... Where \f$ p \f$ is
   * the number of parameters.
   */
  std::vector<value_type> compute_flat_downset() const { return _compute_flat_set(deathCorners_); }

  /**
   * @brief Identifies/merges all birth corners whose maximal coordinate difference is smaller than the
   * given threshold. The new corner takes all minimal coordinates of the two, covering therefore both.
   *
   * Note that the merge is done in ordered pairs, so if \f$ b_1 \f$ is close enough to \f$ b_2 \f$ and \f$ b_2 \f$
   * close enough to \f$ b_3 \f$, once \f$ b_1 \f$ and \f$ b_2 \f$ merge, the new point is eventually not merged with
   * \f$ b_3 \f$.
   *
   * @param precision Distance threshold.
   */
  void identify_births(value_type precision) {
    if (!birthCorners_.is_finite()) return;

    const value_type error = 0.99;

    for (Index i = 0; i < birthCorners_.num_generators(); i++) {
      for (Index j = i + 1; j < birthCorners_.num_generators(); j++) {
        if (_d_inf(birthCorners_[i], birthCorners_[j]) < error * precision) {  // for machine error ?
          _factorize_min(birthCorners_[i], birthCorners_[j]);
          birthCorners_[j] = Births::Generator::inf();
          i++;
        }
      }
    }
    birthCorners_.simplify();
  }

  /**
   * @brief Identifies/merges all death corners whose maximal coordinate difference is smaller than the
   * given threshold. The new corner takes all maximal coordinates of the two, covering therefore both.
   *
   * Note that the merge is done in ordered pairs, so if \f$ b_1 \f$ is close enough to \f$ b_2 \f$ and \f$ b_2 \f$
   * close enough to \f$ b_3 \f$, once \f$ b_1 \f$ and \f$ b_2 \f$ merge, the new point is eventually not merged with
   * \f$ b_3 \f$.
   *
   * @param precision Distance threshold.
   */
  void identify_deaths(value_type precision) {
    if (!deathCorners_.is_finite()) return;

    const value_type error = 0.99;

    for (Index i = 0; i < deathCorners_.num_generators(); i++) {
      for (Index j = i + 1; j < deathCorners_.num_generators(); j++) {
        if (_d_inf(deathCorners_[i], deathCorners_[j]) < error * precision) {  // for machine error ?
          _factorize_max(deathCorners_[i], deathCorners_[j]);
          deathCorners_[j] = Deaths::Generator::minus_inf();
          i++;
        }
      }
    }
    deathCorners_.simplify();
  }

  /**
   * @brief Computes the two end points of the intersection between the summand and the given line.
   * Returns \f$ {inf, inf} \f$ if the intersection is empty.
   *
   * As the line has to have a positive slope, it can represent a 1-parameter sub-filtration and the intersection
   * is a bar in its barcode.
   *
   * @param l Line to intersect with. It has to have the same number of coordinates/parameters then the summand.
   * @return The two endpoints as an `std::array<T,2>`.
   */
  std::array<value_type, 2> get_bar(const Line<value_type> &l) const {
    value_type pushedBirth = T_inf;
    value_type pushedDeath = T_m_inf;

    for (const auto &birth : birthCorners_) {
      value_type pb = l.compute_forward_intersection(birth.begin(), birth.end());
      pushedBirth = std::min(pb, pushedBirth);
    }
    for (const auto &death : deathCorners_) {
      value_type pd = l.compute_backward_intersection(death.begin(), death.end());
      pushedDeath = std::max(pd, pushedDeath);
    }

    // !(<=) in case there is a NaN ?
    if (!(pushedBirth <= pushedDeath)) {
      return {T_inf, T_inf};
    }
    return {pushedBirth, pushedDeath};
  }

  /**
   * @brief Adds the given birth and death corners to the summand.
   *
   * Note that if a new corner is overshadowed by an already existing corner, the new corner will be ignored and will
   * not appear in the corner range.
   *
   * @tparam BirthGeneratorRange Range of elements convertible to `T`. Must have a begin(), end() method and the
   * iterator type should satisfy the requirements of the standard `LegacyForwardIterator`. Note that
   * @ref Births "Births::Generator" does fullfill those requirements.
   * @tparam DeathGeneratorRange  Range of elements convertible to `T`. Must have a begin(), end() method and the
   * iterator type should satisfy the requirements of the standard `LegacyForwardIterator`. Note that
   * @ref Deaths "Deaths::Generator" does fullfill those requirements.
   * @param birth Birth corner.
   * @param death Death corner.
   */
  template <class BirthGeneratorRange = std::initializer_list<T>, class DeathGeneratorRange = std::initializer_list<T>>
  void add_bar(const BirthGeneratorRange &birth, const DeathGeneratorRange &death) {
    birthCorners_.add_generator(birth);
    deathCorners_.add_generator(death);
  }

  /**
   * @brief Rescales the summand's corners for each parameter, that is, for each corner in the summand,
   * the \f$ p^{th} \f$ coordinate is multiplied by `rescaleFactors[p]`.
   *
   * @tparam RandomAccessValueRange Range with a size() and operator[] method.
   * @param rescaleFactors Rescale factors. There must be at least as many elements than parameters in the summand.
   */
  template <class RandomAccessValueRange>
  void rescale(const RandomAccessValueRange &rescaleFactors) {
    _transform(rescaleFactors, [](value_type &cornerVal, value_type fact) -> value_type { return cornerVal *= fact; });
  }

  /**
   * @brief Translates each parameter of the summand, that is, for each corner in the summand,
   * the \f$ p^{th} \f$ coordinate is summed with `translation[p]`.
   *
   * @tparam RandomAccessValueRange Range with a size() and operator[] method.
   * @param translation Translation factors. There must be at least as many elements than parameters in the summand.
   */
  template <class RandomAccessValueRange>
  void translate(const RandomAccessValueRange &translation) {
    _transform(translation, [](value_type &cornerVal, value_type fact) -> value_type { return cornerVal += fact; });
  }

  /**
   * @brief Snaps the summand's corner coordinates to their closest integer and interprets them as indices in a grid
   * to replaces them with the values at those indices in the given grid. For example, if \f$ c \f$ is a birth or
   * death corner in the summand, its new value at parameter \f$ p \f$ will be \f$ c[p] = grid[p][snap(c[p])] \f$.
   * If \f$ snap(c[p]) \f$ is out of bound, the value is replaced with infinity.
   *
   * @tparam GridRange Range with size() and operator[] method. The operator[] method must return a type with the
   * same methods and a value type convertible to `T`.
   * @param grid 2-dimensional range with first axis at least as long as number of parameters.
   */
  template <class GridRange>
  void evaluate_in_grid(const GridRange &grid) {
    if (birthCorners_.num_generators() == 0) return;

    GUDHI_CHECK(grid.size() >= birthCorners_.num_parameters(),
                std::invalid_argument(
                    "The size of the grid should correspond to the number of parameters in the filtration value."));

    auto snap = [](value_type x) -> std::size_t {
      value_type a = std::floor(x);
      value_type b = std::ceil(x);
      if (x - a < b - x) return a;
      return b;
    };
    auto evaluate_generator = [&](auto &corners) {
      return [&corners, &grid, &snap](std::size_t g) {
        for (Index p = 0; p < corners.num_parameters(); ++p) {
          GUDHI_CHECK(corners(g, p) >= 0, std::runtime_error("Values in the corners have to be positive."));
          if (corners(g, p) != T_inf) {
            std::size_t snapped = snap(corners(g, p));
            corners(g, p) = (snapped >= grid[p].size() ? T_inf : grid[p][snapped]);
          }
        }
      };
    };

#ifdef GUDHI_USE_TBB
    tbb::parallel_for(Index(0), birthCorners_.num_generators(), evaluate_generator(birthCorners_));
    if (deathCorners_.num_generators() == 0) return;
    tbb::parallel_for(Index(0), deathCorners_.num_generators(), evaluate_generator(deathCorners_));
#else
    for (Index g = 0; g < birthCorners_.num_generators(); ++g) {
      evaluate_generator(birthCorners_)(g);
    }
    for (Index g = 0; g < deathCorners_.num_generators(); ++g) {
      evaluate_generator(deathCorners_)(g);
    }
#endif
  }

  /**
   * @brief Returns the default value of the Summand type `Y`. E.g., `get_null_value<Summand::Dimension>()` corresponds
   * to the value of the dimension of the Summand if not explicitly set.
   */
  template <typename Y>
  static constexpr Y get_null_value() {
    return -1;
  }

  /**
   * @brief Equality operator.
   */
  friend bool operator==(const Summand &a, const Summand &b) {
    return a.dimension_ == b.dimension_ && a.birthCorners_ == b.birthCorners_ && a.deathCorners_ == b.deathCorners_;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream &operator<<(std::ostream &stream, const Summand &s) {
    stream << "Summand:\n";
    stream << "Dimension: " << s.dimension_ << "\n";
    stream << "Birth corners:\n";
    stream << s.birthCorners_ << "\n";
    stream << "Death corners:\n";
    stream << s.deathCorners_ << "\n";

    return stream;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Summand &sum1, Summand &sum2) noexcept {
    swap(sum1.birthCorners_, sum2.birthCorners_);
    swap(sum1.deathCorners_, sum2.deathCorners_);
    std::swap(sum1.dimension_, sum2.dimension_);
  }

 private:
  Births birthCorners_; /**< Birth corner range. */
  Deaths deathCorners_; /**< Death corner range. */
  Dimension dimension_; /**< Dimension. */

  template <class RandomAccessValueRange, class F>
  void _transform(const RandomAccessValueRange &factors, F &&operate) {
    if (birthCorners_.num_generators() == 0) return;

    auto dimension = birthCorners_.num_parameters();

    GUDHI_CHECK(dimension <= factors.size(), std::invalid_argument("Not enough factors in input."));

    // TODO: parallelize?
    for (unsigned int g = 0; g < birthCorners_.num_generators(); ++g) {
      for (unsigned int p = 0; p < dimension; ++p) {
        std::forward<F>(operate)(birthCorners_(g, p), factors[p]);
      }
    }
    for (unsigned int g = 0; g < deathCorners_.num_generators(); ++g) {
      for (unsigned int p = 0; p < dimension; ++p) {
        std::forward<F>(operate)(deathCorners_(g, p), factors[p]);
      }
    }
  }

  template <class Corners>
  static std::vector<value_type> _compute_flat_set(const Corners &corners) {
    std::vector<value_type> res(corners.num_generators() * corners.num_parameters());
    Simple_mdspan view(res.data(), corners.num_generators(), corners.num_parameters());
    for (Index g = 0; g < corners.num_generators(); ++g) {
      for (Index p = 0; p < corners.num_parameters(); ++p) {
        view(g, p) = corners(g, p);
      }
    }
    return res;
  }

  // TODO: better name?
  template <class Generator>
  static value_type _d_inf(const Generator &a, const Generator &b) {
    if (a.size() == 0 || b.size() == 0 || a.size() != b.size()) return T_inf;

    value_type d = std::abs(a[0] - b[0]);
    for (Index i = 1; i < a.size(); i++) d = std::max(d, std::abs(a[i] - b[i]));

    return d;
  }

  template <class Generator>
  static void _factorize_min(Generator &a, const Generator &b) {
    for (Index i = 0; i < std::min(b.size(), a.size()); i++) a[i] = std::min(a[i], b[i]);
  }

  template <class Generator>
  static void _factorize_max(Generator &a, const Generator &b) {
    for (Index i = 0; i < std::min(b.size(), a.size()); i++) a[i] = std::max(a[i], b[i]);
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SUMMAND_H_INCLUDED
