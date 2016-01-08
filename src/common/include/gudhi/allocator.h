/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015  INRIA Saclay - Ile de France
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
