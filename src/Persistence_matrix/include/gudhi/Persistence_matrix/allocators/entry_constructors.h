/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file entry_constructors.h
 * @author Hannah Schreiber
 * @brief Contains different versions of @ref Gudhi::persistence_matrix::Entry factories.
 */

#ifndef PM_COLUMN_ENTRY_CONSTRUCTORS_H
#define PM_COLUMN_ENTRY_CONSTRUCTORS_H

#include <utility>  //std::swap

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief @ref Entry factory. Constructs and destroys entry pointers with new and delete.
 * 
 * @tparam Entry @ref Entry with the right templates.
 */
template <class Entry>
struct New_entry_constructor 
{
  /**
   * @brief Default constructor.
   */
  New_entry_constructor() {}

  /**
   * @brief Constructs an entry with the given entry arguments.
   * 
   * @param u Arguments forwarded to the @ref Entry constructor.
   * @return @ref Entry pointer.
   */
  template <class... U>
  Entry* construct(U&&... u) const {
    return new Entry(std::forward<U>(u)...);
  }

  /**
   * @brief Destroys the given entry.
   * 
   * @param entry @ref Entry pointer.
   */
  void destroy(Entry* entry) const { delete entry; }

  /**
   * @brief Swap operator.
   */
  friend void swap(New_entry_constructor& col1, New_entry_constructor& col2) {}
};

/**
 * @private
 * @ingroup persistence_matrix
 *
 * @brief @ref Entry factory. Uses @ref Gudhi::Simple_object_pool, which is based on boost::object_pool,
 * to construct and destroy entry pointer.
 * 
 * @tparam Entry @ref Entry with the right templates.
 */
template <class Entry>
struct Pool_entry_constructor 
{
 public:
  /**
   * @brief Default constructor.
   * 
   */
  Pool_entry_constructor() : entryPool_() {}
  //TODO: what does happen when the pool is copied?
  /**
   * @brief Copy constructor.
   * 
   * @param col Factory to copy.
   */
  Pool_entry_constructor(const Pool_entry_constructor& col) : entryPool_(col.entryPool_) {}
  /**
   * @brief Move constructor.
   * 
   * @param col Factory to move.
   */
  Pool_entry_constructor(Pool_entry_constructor&& col) : entryPool_(std::move(col.entryPool_)) {}

  /**
   * @brief Constructs an entry with the given entry arguments.
   * 
   * @param u Arguments forwarded to the @ref Entry constructor.
   * @return @ref Entry pointer.
   */
  template <class... U>
  Entry* construct(U&&... u) {
    return entryPool_.construct(std::forward<U>(u)...);
  }

  /**
   * @brief Destroys the given entry.
   * 
   * @param entry @ref Entry pointer.
   */
  void destroy(Entry* entry) { entryPool_.destroy(entry); }

  //TODO: Again, what does it mean to copy the pool?
  /**
   * @brief Assign operator.
   */
  Pool_entry_constructor& operator=(const Pool_entry_constructor& other) {
    entryPool_ = other.entryPool_;
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Pool_entry_constructor& col1, Pool_entry_constructor& col2) {
    std::swap(col1.entryPool_, col2.entryPool_);
  }

 private:
  Simple_object_pool<Entry> entryPool_;   /**< Entry pool. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_COLUMN_ENTRY_CONSTRUCTORS_H
