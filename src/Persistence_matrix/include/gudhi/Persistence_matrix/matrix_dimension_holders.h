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
 * @file matrix_dimension_holders.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Matrix_max_dimension_holder
 * @ref Gudhi::persistence_matrix::Matrix_all_dimension_holder classes
 * and the @ref Gudhi::persistence_matrix::Dummy_matrix_dimension_holder structure.
 */

#ifndef PM_MATRIX_DIM_HOLDER_H
#define PM_MATRIX_DIM_HOLDER_H

#include <utility>  //std::swap, std::move & std::exchange
#include <vector>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Matrix_max_dimension_holder or @ref Matrix_all_dimension_holder, when the maximal
 * dimension of a matrix is not stored.
 */
struct Dummy_matrix_dimension_holder
{
  template <typename Dimension>
  Dummy_matrix_dimension_holder([[maybe_unused]] Dimension maximalDimension) {}

  friend void swap([[maybe_unused]] Dummy_matrix_dimension_holder& d1,
                   [[maybe_unused]] Dummy_matrix_dimension_holder& d2) {}
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the maximal dimension of a cell represented in the inheriting matrix, when the option of
 * cell removal is not enabled.
 * 
 * @tparam Dimension Dimension value type. Has to be an integer type.
 * If unsigned, the maximal value of the type should not be attained during a run.
 */
template <typename Dimension>
class Matrix_max_dimension_holder 
{
 public:
  /**
   * @brief Default constructor. If a dimension is specified, stores it as the maximal value.
   * 
   * @param maximalDimension Value of the maximal dimension. Has to be either positive or -1. Default value: -1.
   */
  Matrix_max_dimension_holder(Dimension maximalDimension = -1) : maxDim_(maximalDimension){};
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Matrix to copy.
   */
  Matrix_max_dimension_holder(const Matrix_max_dimension_holder& toCopy) : maxDim_(toCopy.maxDim_){};
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Matrix_max_dimension_holder(Matrix_max_dimension_holder&& other) noexcept
      : maxDim_(std::exchange(other.maxDim_, -1)){};

  /**
   * @brief Returns the maximal dimension of a cell represented in the matrix.
   * 
   * @return The maximal dimension.
   */
  Dimension get_max_dimension() const { return maxDim_; };

  /**
   * @brief Assign operator.
   */
  Matrix_max_dimension_holder& operator=(const Matrix_max_dimension_holder& other) {
    std::swap(maxDim_, other.maxDim_);
    return *this;
  };
  /**
   * @brief Swap operator.
   */
  friend void swap(Matrix_max_dimension_holder& matrix1, Matrix_max_dimension_holder& matrix2) {
    std::swap(matrix1.maxDim_, matrix2.maxDim_);
  }

 protected:
  Dimension maxDim_; /**< Current maximal dimension. */

  void update_up(Dimension dimension) {
    if (maxDim_ == -1 || maxDim_ < dimension) maxDim_ = dimension;
  };
};

/**
 * @ingroup persistence_matrix
 *
 * @brief Class managing the maximal dimension of a cell represented in the inheriting matrix, when the option of
 * cell removal is enabled.
 * 
 * @tparam Dimension Dimension value type. Has to be an integer type.
 * If unsigned, the maximal value of the type should not be attained during a run.
 */
template <typename Dimension>
class Matrix_all_dimension_holder 
{
 public:
  /**
   * @brief Default constructor. If a dimension is specified, stores it as the maximal value.
   * 
   * @param maximalDimension Value of the maximal dimension. Has to be either positive or -1. Default value: -1.
   */
  Matrix_all_dimension_holder(Dimension maximalDimension = -1)
      : dimensions_(maximalDimension < 0 ? 0 : maximalDimension + 1, 0), maxDim_(maximalDimension) {
    if (maxDim_ != -1) dimensions_[maxDim_] = 1;
  };
  /**
   * @brief Copy constructor.
   * 
   * @param toCopy Matrix to copy.
   */
  Matrix_all_dimension_holder(const Matrix_all_dimension_holder& toCopy)
      : dimensions_(toCopy.dimensions_), maxDim_(toCopy.maxDim_){};
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Matrix_all_dimension_holder(Matrix_all_dimension_holder&& other) noexcept
      : dimensions_(std::move(other.dimensions_)), maxDim_(std::exchange(other.maxDim_, -1)){};

  /**
   * @brief Returns the maximal dimension of a cell represented in the matrix.
   * 
   * @return The maximal dimension.
   */
  Dimension get_max_dimension() const { return maxDim_; };

  /**
   * @brief Assign operator.
   */
  Matrix_all_dimension_holder& operator=(Matrix_all_dimension_holder other) {
    std::swap(maxDim_, other.maxDim_);
    dimensions_.swap(other.dimensions_);
    return *this;
  };
  /**
   * @brief Swap operator.
   */
  friend void swap(Matrix_all_dimension_holder& matrix1, Matrix_all_dimension_holder& matrix2) {
    std::swap(matrix1.maxDim_, matrix2.maxDim_);
    matrix1.dimensions_.swap(matrix2.dimensions_);
  }

 protected:
  std::vector<unsigned int> dimensions_;  /**< Number of cells by dimension. */
  Dimension maxDim_;                      /**< Current maximal dimension. */

  void update_up(unsigned int dimension) {
    if (dimensions_.size() <= dimension) dimensions_.resize(dimension + 1, 0);
    ++(dimensions_[dimension]);
    maxDim_ = dimensions_.size() - 1;
  };

  void update_down(unsigned int dimension) {
    --(dimensions_[dimension]);  // assumes dimension already exists and is not 0
    while (!dimensions_.empty() && dimensions_.back() == 0) dimensions_.pop_back();
    maxDim_ = dimensions_.size() - 1;
  };
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_MATRIX_DIM_HOLDER_H
