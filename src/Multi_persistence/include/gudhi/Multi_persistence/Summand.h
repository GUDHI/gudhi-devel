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

#include <cstddef>    // std::size_t
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
 * @brief
 *
 * @tparam T
 */
template <typename T, typename D = int>
class Summand {
 public:
  using value_type = T;
  using Births = Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<value_type, false, false>;
  using Deaths = Gudhi::multi_filtration::Dynamic_multi_parameter_filtration<value_type, true, false>;
  using Dimension = D;
  using Index = typename Births::size_type;

  static constexpr T T_inf = Births::T_inf;
  static constexpr T T_m_inf = Births::T_m_inf;

  Summand(int numberOfParameters = 1)
      : birthCorners_(numberOfParameters, T_inf),
        deathCorners_(numberOfParameters, T_m_inf),
        dimension_(get_null_value<Dimension>()) {}

  Summand(const Births &birthCorners, const Deaths &deathCorners, Dimension dimension)
      : birthCorners_(birthCorners), deathCorners_(deathCorners), dimension_(dimension) {}

  // Builds filtration value with given number of parameters and values from the given range. Lets \f$ p \f$
  // be the number of parameters. The \f$ p \f$ first elements of the range have to correspond to the first generator,
  // the \f$ p \f$ next elements to the second generator and so on... So the length of the range has to be a multiple
  // of \f$ p \f$ and the number of generators will be \f$ length / p \f$. The range is represented by two iterators.
  template <class ValueRange = std::initializer_list<value_type>>
  Summand(const ValueRange &birthCorners, const ValueRange &deathCorners, int numberOfParameters, Dimension dimension)
      : birthCorners_(birthCorners.begin(), birthCorners.end(), numberOfParameters),
        deathCorners_(deathCorners.begin(), deathCorners.end(), numberOfParameters),
        dimension_(dimension) {}

  Dimension get_dimension() const { return dimension_; }

  void set_dimension(Dimension dimension) { dimension_ = dimension; }

  [[nodiscard]] int get_number_of_parameters() const { return birthCorners_.num_parameters(); }

  [[nodiscard]] int get_number_of_birth_corners() const { return birthCorners_.num_generators(); }

  [[nodiscard]] int get_number_of_death_corners() const { return deathCorners_.num_generators(); }

  template <class MultiFiltrationValue>
  bool contains(const MultiFiltrationValue &x) const {
    auto xPos = Gudhi::multi_filtration::as_type<Births>(x);
    auto xNeg = Gudhi::multi_filtration::as_type<Deaths>(x);
    return birthCorners_ <= xPos && deathCorners_ <= xNeg;
  }

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

  const Births &get_upset() const { return birthCorners_; }

  const Deaths &get_downset() const { return deathCorners_; }

  std::vector<value_type> compute_birth_list() const { return _compute_list(birthCorners_); }

  std::vector<value_type> compute_death_list() const { return _compute_list(deathCorners_); }

  void complete_birth(value_type precision) {
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

  void complete_death(value_type precision) {
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

  // TODO: generalize as a "GeneratorRange" like for `add_generator` ?
  // Not usefull for multipers right now, just for C++ interface purposes
  void add_bar(const typename Births::Generator &birth, const typename Deaths::Generator &death) {
    birthCorners_.add_generator(birth);
    deathCorners_.add_generator(death);
  }

  template <class RandomAccessValueRange>
  void rescale(const RandomAccessValueRange &rescaleFactors) {
    _transform(rescaleFactors, [](value_type &cornerVal, value_type fact) -> value_type { return cornerVal *= fact; });
  }

  template <class RandomAccessValueRange>
  void translate(const RandomAccessValueRange &translation) {
    _transform(translation, [](value_type &cornerVal, value_type fact) -> value_type { return cornerVal += fact; });
  }

  template <class GridRange>
  void evaluate_in_grid(const GridRange &grid) {
    if (birthCorners_.num_generators() == 0) return;

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

  template <typename Y>
  static constexpr Y get_null_value() {
    return -1;
  }

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

  friend void swap(Summand &sum1, Summand &sum2) noexcept {
    swap(sum1.birthCorners_, sum2.birthCorners_);
    swap(sum1.deathCorners_, sum2.deathCorners_);
    std::swap(sum1.dimension_, sum2.dimension_);
  }

 private:
  Births birthCorners_;
  Deaths deathCorners_;
  Dimension dimension_;

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
  static std::vector<value_type> _compute_list(const Corners &corners) {
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
