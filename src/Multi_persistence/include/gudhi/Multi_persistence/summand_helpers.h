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
 * @file summand_helpers.h
 * @author David Loiseaux
 * @brief Contains the helper methods @ref Gudhi::multi_persistence::.
 */

#ifndef MP_SUMMAND_HELPERS_H_
#define MP_SUMMAND_HELPERS_H_

#include <algorithm>  //std::max
#include <cstddef>    //std::size_t
#include <stdexcept>  //std::invalid_argument

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Summand.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @private
 */
template <class Summand, class RandomAccessValueRange, class Corners, typename F>
inline typename Summand::value_type _compute_distance_to_front(const RandomAccessValueRange &x, const Corners &corners,
                                                               bool negative, F &&diff) {
  typename Summand::value_type distance = Summand::T_inf;

  for (typename Summand::Index g = 0; g < corners.num_generators(); ++g) {
    typename Summand::value_type tempDist = negative ? Summand::T_m_inf : 0;
    for (typename Summand::Index p = 0; p < corners.num_parameters(); ++p) {
      auto d = std::forward<F>(diff)(corners(g, p), x[p]);
      if (tempDist < d) tempDist = d;
    }
    if (distance > tempDist) distance = tempDist;
  }

  return distance;
}

/**
 * @brief Computes the distance from the given point `x` to the given Summand `sum`.
 * TODO: proper definition of the distance.
 * 
 * @tparam T First template argument of @ref Summand.
 * @tparam D Second template argument of @ref Summand.
 * @tparam RandomAccessValueRange Range of `T` with a size() and operator[] method.
 * @param sum Summand.
 * @param x Point.
 * @param negative If true, the distance is allowed to be signed.
 */
template <typename T, typename D, class RandomAccessValueRange>
inline T compute_summand_distance_to(const Summand<T, D> &sum, const RandomAccessValueRange &x, bool negative) {
  GUDHI_CHECK(x.size() >= sum.get_number_of_parameters(),
              std::invalid_argument("The given point does not have enough coordinates compared to the given Summand."));

  T lowerDist = _compute_distance_to_front<Summand<T, D> >(x, sum.get_upset(), negative,
                                                           [](T cornerVal, T xVal) -> T { return cornerVal - xVal; });
  T upperDist = _compute_distance_to_front<Summand<T, D> >(x, sum.get_downset(), negative,
                                                           [](T cornerVal, T xVal) -> T { return xVal - cornerVal; });
  return std::max(lowerDist, upperDist);
}

/**
 * @private
 */
template <typename T, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline T _get_summand_max_diagonal(const RandomAccessValueRange1 &birth, const RandomAccessValueRange2 &death,
                                   const Box<T> &box = {}) {
  // assumes birth and death to be never NaN
  GUDHI_CHECK(birth.size() == 1 || death.size() == 1 || birth.size() == death.size(),
              std::invalid_argument("Inputs must be of the same size !"));

  bool useThreshold = !box.is_trivial();

  GUDHI_CHECK((birth.size() == 1 && death.size() == 1) || !useThreshold || birth.size() == box.get_dimension() ||
                  death.size() == box.get_dimension(),
              std::invalid_argument("Inputs must be of the same size !"));

  auto get_val = [](const auto &r, std::size_t i) -> T {
    if (i < r.size()) return r[i];
    // never used if r.size() == 0
    return r[0];
  };

  T diag = Summand<T>::T_inf;
  if (useThreshold) {
    for (std::size_t i = 0; i < birth.size(); ++i) {
      T max_i = box.get_upper_corner()[i];
      T min_i = box.get_lower_corner()[i];
      T t_death = std::min(get_val(death, i), max_i);
      T t_birth = std::max(get_val(birth, i), min_i);
      diag = std::min(diag, t_death - t_birth);
    }
  } else {
    for (std::size_t i = 0; i < birth.size(); i++) {
      diag = std::min(diag, get_val(death, i) - get_val(birth, i));
    }
  }

  return diag;
}

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam D 
 * @param sum 
 * @param box 
 * @return T 
 */
template <typename T, typename D>
T compute_summand_interleaving(const Summand<T, D> &sum, const Box<T> &box) {
  T interleaving = 0;
  for (const auto &birth : sum.get_upset()) {
    for (const auto &death : sum.get_downset()) {
      interleaving = std::max(interleaving, _get_summand_max_diagonal<T>(birth, death, box));
    }
  }
  return interleaving;
}

/**
 * @private
 */
template <class BirthGenerator, class DeathGenerator, typename T>
inline T _summand_rectangle_volume(const BirthGenerator &birth, const DeathGenerator &death, const Box<T> &box) {
  // NaN?
  if (birth.size() == 0 || death.size() == 0) return 0;

  auto get_val = [](const auto &r, std::size_t i) -> T {
    if (i < r.size()) return r[i];
    // never used if r.size() == 0
    return r[0];
  };

  T volume = std::min(death[0], box.get_upper_corner()[0]) - std::max(birth[0], box.get_lower_corner()[0]);
  for (std::size_t i = 1; i < birth.size(); i++) {
    T t_death = std::min(get_val(death, i), box.get_upper_corner()[i]);
    T t_birth = std::max(get_val(birth, i), box.get_lower_corner()[i]);
    volume = volume * (t_death - t_birth);
  }
  return volume;
}

template <typename T, typename D, class RandomAccessValueRange>
inline T compute_summand_local_weight(const Summand<T, D> &sum, const RandomAccessValueRange &x, T delta) {
  using P = typename Box<T>::Point_t;

  GUDHI_CHECK(x.size() == sum.get_number_of_parameters(),
              std::invalid_argument("Input range does not have the right size."));

  bool rectangle = delta <= 0;

  // box on which to compute the local weight
  P mini(x.size());
  P maxi(x.size());
  for (typename Summand<T, D>::Index i = 0; i < x.size(); i++) {
    mini[i] = rectangle ? x[i] + delta : x[i] - delta;
    maxi[i] = rectangle ? x[i] - delta : x[i] + delta;
  }
  Box<T> threshold(rectangle ? maxi : mini, rectangle ? mini : maxi);

  T localWeight = 0;

  if (rectangle) {
    // local weight is the volume of the largest rectangle in the restricted
    for (const auto &birth : sum.get_upset()) {
      for (const auto &death : sum.get_downset()) {
        localWeight = std::max(localWeight, _summand_rectangle_volume(birth, death, threshold));
      }
    }
    return localWeight / std::pow(2 * std::abs(delta), x.size());
  }

  // local weight is interleaving to 0 of module restricted to the square
  localWeight = compute_summand_interleaving(sum, threshold);
  return localWeight / (2 * std::abs(delta));
}

template <typename T, typename D, class RandomAccessValueRange>
inline T compute_summand_landscape_value(const Summand<T, D> &sum, const RandomAccessValueRange &x) {
  T landscapeValue = 0;
  for (const auto &birth : sum.get_upset()) {
    for (const auto &death : sum.get_downset()) {
      T value = std::min(_get_summand_max_diagonal<T>(birth, x), _get_summand_max_diagonal<T>(x, death));
      landscapeValue = std::max(landscapeValue, value);
    }
  }
  return landscapeValue;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SUMMAND_HELPERS_H_
