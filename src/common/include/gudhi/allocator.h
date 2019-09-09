/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ALLOCATOR_H_
#define ALLOCATOR_H_

#include <memory>
#include <utility>

namespace Gudhi {

/** \private
 * An allocator that can be used to build an uninitialized vector.
 */
template <class T, class Base = std::allocator<T>>
struct no_init_allocator : Base {
  typedef std::allocator_traits<Base> Base_traits;
  template <class U> struct rebind {
    typedef no_init_allocator<U, typename Base_traits::template rebind_alloc<U>> other;
  };

  // Inherit constructors.
  using Base::Base;

  // Do nothing: that's the whole point!
  template<class P>
  void construct(P*) noexcept {}

  template<class P, class...U> void construct(P*p, U&&...u) {
    Base_traits::construct(*(Base*)this, p, std::forward<U>(u)...);
  }
};

}  // namespace Gudhi

#endif  // ALLOCATOR_H_
