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
 * @file cell_constructors.h
 * @author Hannah Schreiber
 * @brief Contains the @ref New_cell_constructor and @ref Pool_cell_constructor structures.
 */

#ifndef PM_COLUMN_CELL_CONSTRUCTORS_H
#define PM_COLUMN_CELL_CONSTRUCTORS_H

#include <utility>  //std::swap

#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief @ref Cell factory. Constructs and destroyes cell pointers with new and delete.
 * 
 * @tparam Cell @ref Cell with the right templates.
 */
template <class Cell>
struct New_cell_constructor 
{
  /**
   * @brief Default constructor.
   */
  New_cell_constructor() {}

  /**
   * @brief Constructs a cell with the given cell arguments.
   * 
   * @param u Arguments forwarded to the @ref Cell constructor.
   * @return @ref Cell pointer.
   */
  template <class... U>
  Cell* construct(U&&... u) const {
    return new Cell(std::forward<U>(u)...);
  }

  /**
   * @brief Destroyes the given cell.
   * 
   * @param cell @ref Cell pointer.
   */
  void destroy(Cell* cell) const { delete cell; }

  /**
   * @brief Swap operator.
   */
  friend void swap(New_cell_constructor& col1, New_cell_constructor& col2) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief @ref Cell factory. Uses @ref Gudhi::Simple_object_pool, which is based on boost::object_pool,
 * to construct and destroy cell pointer.
 * 
 * @tparam Cell @ref Cell with the right templates.
 */
template <class Cell>
struct Pool_cell_constructor 
{
 public:
  /**
   * @brief Default constructor.
   * 
   */
  Pool_cell_constructor() : cellPool_() {}
  //TODO: what does happen when the pool is copied?
  /**
   * @brief Copy constructor.
   * 
   * @param col Factory to copy.
   */
  Pool_cell_constructor(const Pool_cell_constructor& col) : cellPool_(col.cellPool_) {}
  /**
   * @brief Move constructor.
   * 
   * @param col Factory to move.
   */
  Pool_cell_constructor(Pool_cell_constructor&& col) : cellPool_(std::move(col.cellPool_)) {}

  /**
   * @brief Constructs a cell with the given cell arguments.
   * 
   * @param u Arguments forwarded to the @ref Cell constructor.
   * @return @ref Cell pointer.
   */
  template <class... U>
  Cell* construct(U&&... u) {
    return cellPool_.construct(std::forward<U>(u)...);
  }

  /**
   * @brief Destroyes the given cell.
   * 
   * @param cell @ref Cell pointer.
   */
  void destroy(Cell* cell) { cellPool_.destroy(cell); }

  //TODO: Again, what does it mean to copy the pool?
  /**
   * @brief Assign operator.
   */
  Pool_cell_constructor& operator=(const Pool_cell_constructor& other) {
    cellPool_ = other.cellPool_;
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Pool_cell_constructor& col1, Pool_cell_constructor& col2) {
    std::swap(col1.cellPool_, col2.cellPool_);
  }

 private:
  Simple_object_pool<Cell> cellPool_;   /**< Cell pool. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_COLUMN_CELL_CONSTRUCTORS_H
