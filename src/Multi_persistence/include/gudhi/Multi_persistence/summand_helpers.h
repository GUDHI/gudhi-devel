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
#include <type_traits>

#include <gudhi/Debug_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_persistence/Box.h>
#include <gudhi/Multi_persistence/Summand.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @ingroup multi_persistence
 *
 * @private
 */
template <typename Distance, class Summand, class RandomAccessValueRange, class Corners, typename F>
inline auto _compute_distance_to_front(const RandomAccessValueRange &x, const Corners &corners, bool negative,
                                       F &&diff) {
  Distance distance = Gudhi::multi_filtration::MF_T_inf<Distance>;

  for (typename Summand::Index g = 0; g < corners.num_generators(); ++g) {
    Distance tempDist = negative ? Gudhi::multi_filtration::MF_T_m_inf<Distance> : 0;
    for (typename Summand::Index p = 0; p < corners.num_parameters(); ++p) {
      auto d = std::forward<F>(diff)(corners(g, p), x[p]);
      if (tempDist < d) tempDist = d;
    }
    if (distance > tempDist) distance = tempDist;
  }

  return distance;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Computes the distance from the given point `x` to the given Summand `sum`.
 * TODO: proper definition of the distance.
 *
 * @tparam T First template argument of @ref Summand.
 * @tparam D Second template argument of @ref Summand.
 * @tparam RandomAccessValueRange Range of `T` with a size() and operator[] method.
 * @param sum Summand.
 * @param x Point.
 * @param negative If true, the distance is allowed to be signed.
 *
 * @note If `T` is unsigned, the values in the summand and in the point have to be small enough to fit in the
 * corresponding signed type, otherwise the behaviour is undefined. E.g., if `T` is `unsigned int`, the values
 * should not exceed `INT_MAX`.
 */
template <typename T, typename D, class RandomAccessValueRange>
inline T compute_summand_distance_to(const Summand<T, D> &sum, const RandomAccessValueRange &x, bool negative) {
  GUDHI_CHECK(x.size() >= static_cast<std::size_t>(sum.get_number_of_parameters()),
              std::invalid_argument("The given point does not have enough coordinates compared to the given Summand."));

  T lowerDist, upperDist;
  if constexpr (std::is_unsigned_v<T>) {
    // std::make_signed_t does not compile for T signed and std::conditional evaluates both possibilities,
    // so the case has to be separated...
    using stype = std::make_signed_t<T>;
    // This is a bit unsafe, as the unsigned value can be truncated
    // but I will just assume that the user will not use coordinates that big
    // otherwise I would need to do a quite long case study here, which seems overkill
    // for a case which will probably never happen
    lowerDist = _compute_distance_to_front<stype, Summand<T, D> >(
        x, sum.get_upset(), negative, [](stype cornerVal, stype xVal) { return cornerVal - xVal; });
    upperDist = _compute_distance_to_front<stype, Summand<T, D> >(
        x, sum.get_downset(), negative, [](stype cornerVal, stype xVal) { return xVal - cornerVal; });
  } else {
    lowerDist = _compute_distance_to_front<T, Summand<T, D> >(x, sum.get_upset(), negative,
                                                              [](T cornerVal, T xVal) { return cornerVal - xVal; });
    upperDist = _compute_distance_to_front<T, Summand<T, D> >(x, sum.get_downset(), negative,
                                                              [](T cornerVal, T xVal) { return xVal - cornerVal; });
  }

  return std::max(lowerDist, upperDist);
}

/**
 * @ingroup multi_persistence
 *
 * @private
 */
template <typename T, class RandomAccessValueRange1, class RandomAccessValueRange2>
inline T _get_summand_diagonal(const RandomAccessValueRange1 &birth, const RandomAccessValueRange2 &death,
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
 * @ingroup multi_persistence
 *
 * @brief For a birth and death corner in the summand, let the diagonal between those two be
 * \f$ min\{death[p] - birth[p]\} \f$ for all parameters \f$ p \f$. This method returns the maximal diagonal
 * of all birth-death pairs in the intersection between the summand and the box.
 *
 * @tparam T First template argument of @ref Summand.
 * @tparam D Second template argument of @ref Summand.
 * @param sum Summand.
 * @param box Box to intersect with. The box is ignored if trivial.
 */
template <typename T, typename D>
T compute_summand_interleaving(const Summand<T, D> &sum, const Box<T> &box) {
  T interleaving = 0;
  for (const auto &birth : sum.get_upset()) {
    for (const auto &death : sum.get_downset()) {
      // TODO: if the types of Births and Deaths in Summand changes (to become a template for example)
      // the input to _get_summand_diagonal has to get adapted to it, as it makes use of
      // Dynamic_multi_parameter_filtration::Generator working like a vector
      interleaving = std::max(interleaving, _get_summand_diagonal<T>(birth, death, box));
    }
  }
  return interleaving;
}

/**
 * @ingroup multi_persistence
 *
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

/**
 * @ingroup multi_persistence
 *
 * @brief Computes the local weight of the summand within the box centered at @p x and diagonal two times @p delta.
 * See description of @p delta.
 *
 * @tparam T First template argument of @ref Summand.
 * @tparam D Second template argument of @ref Summand.
 * @tparam RandomAccessValueRange Range of arithmetic type with a size() and operator[] method.
 * @param sum Summand.
 * @param x Center of the box.
 * @param delta Half diameter of the box. If positive, the weight computed is the interleaving distance to 0 of the
 * summand restricted to the box. If negative, the weight is the volume of the largest rectangle spanned by a birth
 * and a death corner of the summand intersected with the box.
 */
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
        // TODO: if the types of Births and Deaths in Summand changes (to become a template for example)
        // _summand_rectangle_volume has to get adapted to it, as it makes use of
        // Dynamic_multi_parameter_filtration::Generator working like a vector
        localWeight = std::max(localWeight, _summand_rectangle_volume(birth, death, threshold));
      }
    }
    return localWeight / std::pow(2 * std::abs(delta), x.size());
  }

  // local weight is interleaving to 0 of module restricted to the square
  localWeight = compute_summand_interleaving(sum, threshold);
  return localWeight / (2 * std::abs(delta));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Computes the local landscape value of the summand at the given point.
 *
 * @tparam T First template argument of @ref Summand.
 * @tparam D Second template argument of @ref Summand.
 * @tparam RandomAccessValueRange Range of arithmetic type with a size() and operator[] method.
 * @param sum Summand.
 * @param x Local point. Assumed to have as many coordinates than there are parameters in the summand.
 */
template <typename T, typename D, class RandomAccessValueRange>
inline T compute_summand_landscape_value(const Summand<T, D> &sum, const RandomAccessValueRange &x) {
  T landscapeValue = 0;
  for (const auto &birth : sum.get_upset()) {
    for (const auto &death : sum.get_downset()) {
      T value = std::min(_get_summand_diagonal<T>(birth, x), _get_summand_diagonal<T>(x, death));
      landscapeValue = std::max(landscapeValue, value);
    }
  }
  return landscapeValue;
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SUMMAND_HELPERS_H_
