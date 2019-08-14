/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLE_OBJECT_POOL_H_
#define SIMPLE_OBJECT_POOL_H_

#include <boost/pool/pool.hpp>
#include <utility>

namespace Gudhi {

/** \private
 * This is a simpler version of boost::object_pool, that requires
 * that users explicitly destroy all objects. This lets the
 * performance scale much better, see
 * https://svn.boost.org/trac/boost/ticket/3789 .
 */
template <class T>
class Simple_object_pool : protected boost::pool<boost::default_user_allocator_malloc_free> {
 protected:
  typedef boost::pool<boost::default_user_allocator_malloc_free> Base;
  typedef T* pointer;

  Base& base() {
    return *this;
  }

  Base const& base()const {
    return *this;
  }

 public:
  typedef T element_type;
  typedef boost::default_user_allocator_malloc_free user_allocator;
  typedef typename Base::size_type size_type;
  typedef typename Base::difference_type difference_type;

  template<class...U>
  Simple_object_pool(U&&...u) : Base(sizeof (T), std::forward<U>(u)...) { }

  template<class...U>
  pointer construct(U&&...u) {
    void* p = base().malloc BOOST_PREVENT_MACRO_SUBSTITUTION();
    assert(p);
    try {
      new(p) T(std::forward<U>(u)...);
    }    catch (...) {
      base().free BOOST_PREVENT_MACRO_SUBSTITUTION(p);
      throw;
    }
    return static_cast<pointer> (p);
  }

  void destroy(pointer p) {
    p->~T();
    base().free BOOST_PREVENT_MACRO_SUBSTITUTION(p);
  }
};

}  // namespace Gudhi

#endif  // SIMPLE_OBJECT_POOL_H_
