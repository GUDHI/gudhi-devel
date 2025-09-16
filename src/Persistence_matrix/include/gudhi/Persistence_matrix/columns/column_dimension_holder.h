/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-24 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file column_dimension_holder.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Column_dimension_holder class and
 * @ref Gudhi::persistence_matrix::Dummy_dimension_holder structure.
 */

#ifndef PM_COLUMN_DIM_HOLDER_H
#define PM_COLUMN_DIM_HOLDER_H

#include <utility>  //std::swap

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Column_dimension_holder, when the columns are not storing a dimension.
 */
struct Dummy_dimension_holder {
  Dummy_dimension_holder() = default;

  template <typename Dimension>
  Dummy_dimension_holder([[maybe_unused]] Dimension dim)
  {}

  friend void swap([[maybe_unused]] Dummy_dimension_holder& col1,
                   [[maybe_unused]] Dummy_dimension_holder& col2) noexcept
  {}
};

/**
 * @class Column_dimension_holder column_dimension_holder.h gudhi/Persistence_matrix/columns/column_dimension_holder.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the dimension access of a column.
 *
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
struct Column_dimension_holder {
  using Dimension = typename Master_matrix::Dimension; /**< Dimension value type. */

  /**
   * @brief Default constructor. Sets the dimension to 0 for @ref boundarymatrix "boundary matrices" and to
   * @ref Matrix::get_null_value "null index" for @ref chainmatrix "chain matrices".
   */
  Column_dimension_holder()
      : dim_(Master_matrix::Option_list::is_of_boundary_type ? 0 : Master_matrix::template get_null_value<Dimension>())
  {}

  /**
   * @brief Constructor setting the dimension to the given value.
   *
   * @param dim Dimension of the column.
   */
  Column_dimension_holder(Dimension dim) : dim_(dim) {}

  /**
   * @brief Copy constructor.
   *
   * @param col Column to copy.
   */
  Column_dimension_holder(const Column_dimension_holder& col) = default;

  /**
   * @brief Move constructor.
   *
   * @param col Column to move.
   */
  Column_dimension_holder(Column_dimension_holder&& col) noexcept
      : dim_(std::exchange(col.dim_, Master_matrix::template get_null_value<Dimension>()))
  {}

  ~Column_dimension_holder() = default;

  /**
   * @brief Returns the dimension of the column.
   *
   * @return The dimension of the column.
   */
  Dimension get_dimension() const { return dim_; }

  /**
   * @brief Assign operator.
   */
  Column_dimension_holder& operator=(const Column_dimension_holder& other) = default;

  /**
   * @brief Move assign operator.
   */
  Column_dimension_holder& operator=(Column_dimension_holder&& other) noexcept
  {
    if (this == &other) return *this;

    dim_ = std::exchange(other.dim_, Master_matrix::template get_null_value<Dimension>());
    return *this;
  }

  /**
   * @brief Swap operator.
   */
  friend void swap(Column_dimension_holder& col1, Column_dimension_holder& col2) noexcept
  {
    std::swap(col1.dim_, col2.dim_);
  }

 protected:
  void _swap_dimension(Column_dimension_holder& other) { std::swap(dim_, other.dim_); }

 private:
  Dimension dim_; /**< Dimension of the column. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_COLUMN_DIM_HOLDER_H
