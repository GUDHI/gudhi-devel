/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which
 * is released under MIT. See file LICENSE or go to
 * https://gudhi.inria.fr/licensing/ for full license details. Author(s): David
 * Loiseaux
 *
 *    tCopyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef LINE_FILTRATION_TRANSLATION_H_INCLUDED
#define LINE_FILTRATION_TRANSLATION_H_INCLUDED

#include "Box.h"
#include "Finitely_critical_filtrations.h"
#include <cstddef>

namespace Gudhi::multiparameter::multi_filtrations {

template <typename T> class Line {
public:
  using point_type = One_critical_filtration<T>;
  using kcritical_point_type = Multi_critical_filtration<T>;
  Line();
  Line(const point_type &x);
  Line(point_type &&x);
  Line(const point_type &x, const point_type &v);
  inline point_type push_forward(point_type x) const;
  template <typename U=T>
  inline U push_forward2(const point_type &x) const;
  template <typename U=T>
  inline U push_forward2(const kcritical_point_type &x) const;
  inline point_type push_back(point_type x) const;
  template <typename U=T>
  inline U push_back2(const point_type &x) const;
  template <typename U=T>
  inline U push_back2(const kcritical_point_type &x) const;
  inline int get_dim() const;
  std::pair<point_type, point_type> get_bounds(const Box<T> &box) const;
  std::pair<T,T> get_bounds2(const Box<T> &box) const;

  // translation
  inline friend Line &operator+=(Line &to_translate, const point_type &x) {
    to_translate.basepoint_ -= x;
    return to_translate;
  }

  inline point_type &basepoint() { return basepoint_; }
  inline point_type &direction() { return direction_; }
  inline const point_type &basepoint() const { return basepoint_; }
  inline const point_type &direction() const { return direction_; }

private:
  point_type basepoint_; // any point on the line
  point_type direction_; // direction of the line
};
template <typename T> Line<T>::Line() {}

template <typename T> Line<T>::Line(const point_type &x) : basepoint_(x) {}
template <typename T>
Line<T>::Line(point_type &&x) : basepoint_(std::move(x)) {}
template <typename T>
Line<T>::Line(const point_type &x, const point_type &v)
    : basepoint_(x), direction_(v) {}
template <typename T>
inline typename Line<T>::point_type
Line<T>::push_forward(point_type x) const { // TODO remove copy
  if (x.is_inf() || x.is_nan() || x.is_minus_inf())
    return x;
  T t = this->push_forward2<T>(x);
  if (direction_.size() > 0) {
    for (std::size_t i = 0; i < x.size(); i++)
      x[i] = basepoint_[i] + t * direction_[i];
  } else {
    for (std::size_t i = 0; i < x.size(); i++)
      x[i] = basepoint_[i] + t;
  }
  return x;
}
template <typename T>
template <typename U>
inline U Line<T>::push_forward2(const point_type &x) const {
  constexpr const U inf = std::numeric_limits<U>::infinity(); // This disable it for e.g.
                                                    // ints, but that's good
  if (x.is_inf() || x.is_nan())
    return inf;
  if (x.is_minus_inf())
    return -inf;
  // x -= basepoint_;
  U t = -inf;
  if (direction_.size()) {
    for (std::size_t i = 0; i < x.size(); i++) {
      if (direction_[i] == 0) [[unlikely]] {
        if (x[i] < basepoint_[i])
          continue;
        else {
          return inf;
        }
      } else [[likely]] {
        // the cast float -> float should not be a overhead (if compiler is not
        // stupid)
        t = std::max(t, (static_cast<U>(x[i]) - static_cast<U>(basepoint_[i])) / static_cast<U>((direction_[i])));
      }
    }
  } else {
    for (std::size_t i = 0; i < x.size(); i++)
      t = std::max(t, static_cast<U>(x[i]) - static_cast<U>(basepoint_[i]));
  }

  // for (std::size_t i = 0; i < x.size(); i++) {
  //   T dir = this->direction_.size() > i ? direction_[i] : 1;
  //   T scaled_coord;
  //   if (dir == 0) [[unlikely]] {
  //     scaled_coord = x[i] > basepoint_[i] ? inf : -inf;
  //   } else [[likely]] {
  //     scaled_coord = (x[i] - basepoint_[i]) / dir;
  //   }
  //   t = std::max(t, scaled_coord);
  // }
  return t;
}
template <typename T>
template <typename U>
inline U Line<T>::push_forward2(const kcritical_point_type &x) const {
  constexpr const U inf = std::numeric_limits<U>::infinity();
  if (x.is_inf() || x.is_nan())
    return inf;
  if (x.is_minus_inf())
    return -inf;
  U t = inf;
  for (const auto &y : x) {
    t = std::min(t, this->push_forward2<U>(y));
  }
  return t;
}

template <typename T>
inline typename Line<T>::point_type Line<T>::push_back(point_type x) const {
  if (x.is_inf() || x.is_nan() || x.is_minus_inf())
    return x;

  T t = this->push_back2(x);
  if (direction_.size() > 0) {
    for (std::size_t i = 0; i < x.size(); i++)
      x[i] = basepoint_[i] + t * direction_[i];
  } else
    for (std::size_t i = 0; i < x.size(); i++)
      x[i] = basepoint_[i] + t;

  // for (std::size_t i = 0; i < x.size(); i++)
  //   x[i] =
  //       basepoint_[i] + t * (this->direction_.size() > i ? direction_[i] :
  //       1);
  return x;
}

template <typename T>
template <typename U> inline U Line<T>::push_back2(const point_type &x) const {
  constexpr const  U inf = std::numeric_limits<U>::infinity();
  if (x.is_inf())
    return inf;
  if (x.is_minus_inf() || x.is_nan())
    return -inf;
  U t = inf;
  // x -= basepoint_;

  if (direction_.size()) {
    for (std::size_t i = 0; i < x.size(); i++) {
      if (direction_[i] == 0) [[unlikely]] {
        if (x[i] > basepoint_[i])
          continue;
        else {
          return -inf;
        }
      } else [[likely]] {
        t = std::min(t, (static_cast<U>(x[i])  - static_cast<U>(basepoint_[i])) / static_cast<U>(direction_[i]));
      }
    }
  } else {
    for (std::size_t i = 0; i < x.size(); i++)
      t = std::min(t, static_cast<U>(x[i] - basepoint_[i]));
  }
  // for (std::size_t i = 0; i < x.size(); i++) {
  //   T dir = this->direction_.size() > i ? direction_[i] : 1;
  //   T scaled_coord;
  //   if (dir == 0) [[unlikely]] {
  //     scaled_coord = x[i] > basepoint_[i] ? inf : -inf;
  //   } else [[likely]] {
  //     scaled_coord = (x[i] - basepoint_[i]) / dir;
  //   }
  //   t = std::min(t, scaled_coord);
  // }
  return t;
}

template <typename T>

template <typename U>
inline U Line<T>::push_back2(const kcritical_point_type &x) const {
  constexpr const U inf = std::numeric_limits<U>::infinity();
  if (x.is_inf())
    return inf;
  if (x.is_minus_inf() || x.is_nan())
    return -inf;
  U t = -inf;
  for (const auto &y : x) {
    t = std::max(t, this->push_back2<U>(y));
  }
  return t;
}

template <typename T> inline int Line<T>::get_dim() const {
  return basepoint_.size();
}

template <typename T>
inline std::pair<T,T>
Line<T>::get_bounds2(const Box<T> &box) const {
  return {this->push_forward2(box.get_bottom_corner()),
          this->push_back2(box.get_upper_corner())};
}


template <typename T>
inline std::pair<typename Line<T>::point_type, typename Line<T>::point_type>
Line<T>::get_bounds(const Box<T> &box) const {
  return {this->push_forward(box.get_bottom_corner()),
          this->push_back(box.get_upper_corner())};
}
} // namespace Gudhi::multiparameter::multi_filtrations

#endif // LINE_FILTRATION_TRANSLATION_H_INCLUDED
