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
 * @file chain_column_extra_properties.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Chain_column_extra_properties class and @ref Dummy_chain_properties structure.
 */

#ifndef PM_CHAIN_COLUMN_PROP_H
#define PM_CHAIN_COLUMN_PROP_H

#include <utility>  //std::swap

namespace Gudhi {
namespace persistence_matrix {

/**
 * @brief Empty structure.
 * Inheritated instead of @ref Chain_column_extra_properties, when the columns are not meant for chain matrices.
 */
struct Dummy_chain_properties 
{
  // Dummy_chain_properties() {}
  // Dummy_chain_properties([[maybe_unused]] int pivot) {}
  Dummy_chain_properties([[maybe_unused]] int pivot = 0, [[maybe_unused]] int pair = 0) {}
  // Dummy_chain_properties([[maybe_unused]] const Dummy_chain_properties& col) {}
  // Dummy_chain_properties([[maybe_unused]] Dummy_chain_properties&& col) {}

  // Dummy_chain_properties& operator=([[maybe_unused]] const Dummy_chain_properties& other) { return *this; }

  friend void swap([[maybe_unused]] Dummy_chain_properties& col1, [[maybe_unused]] Dummy_chain_properties& col2) {}
};

/**
 * @brief Class managing the pivot and partitioning of columns in @ref Chain_matrix.
 *
 * The columns of a chain matrix are partitioned in three sets: \f$ F \f$, \f$ G \f$ and \f$ H \f$ with a 
 * bijection between \f$ G \f$ and \f$ H \f$. If a column is in \f$ F \f$, the value of
 * @ref Chain_column_extra_properties::pairedColumn_ is -1, while the value corresponds to the MatIdx index of 
 * the image of the bijection if the column is in either \f$ G \f$ or \f$ H \f$. See [TODO: zigzag paper] for
 * more details.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_column_extra_properties 
{
 public:
  using index = typename Master_matrix::index;        /**< MatIdx index type. */
  using id_index = typename Master_matrix::id_index;  /**< IDIdx index type. */

  /**
   * @brief Default constructor. Sets the pivot and pair to -1, which means "not existing".
   */
  Chain_column_extra_properties() : pivot_(-1), pairedColumn_(-1) {}
  /**
   * @brief Constructor setting the pivot at the given value and the pair to -1 (i.e. not paired).
   * 
   * @param pivot Row index of the pivot. Corresponds to the IDIdx index of the face represented by the column.
   */
  Chain_column_extra_properties(id_index pivot) : pivot_(pivot), pairedColumn_(-1) {}
  /**
   * @brief Constructor setting the pivot and the pair at the given values.
   * 
   * @param pivot Row index of the pivot. Corresponds to the IDIdx index of the face represented by the column.
   * @param pair MatIdx index of the pair of the column.
   */
  Chain_column_extra_properties(id_index pivot, index pair) : pivot_(pivot), pairedColumn_(pair) {}
  /**
   * @brief Copy constructor.
   * 
   * @param col Column to copy.
   */
  Chain_column_extra_properties(const Chain_column_extra_properties& col)
      : pivot_(col.pivot_), pairedColumn_(col.pairedColumn_) {}
  /**
   * @brief Move constructor.
   * 
   * @param col Column to move.
   */
  Chain_column_extra_properties(Chain_column_extra_properties&& col)
      : pivot_(std::exchange(col.pivot_, -1)), pairedColumn_(std::exchange(col.pairedColumn_, -1)) {}

  /**
   * @brief Returns -1 if the column is not paired, the matIdx of the pair otherwise.
   * 
   * @return -1 if the column is not paired, the matIdx of the pair otherwise.
   */
  index get_paired_chain_index() const { return pairedColumn_; }
  /**
   * @brief Indicates if the column is paired or not.
   * 
   * @return true If the column is paired.
   * @return false Otherwise.
   */
  bool is_paired() const { return pairedColumn_ != static_cast<index>(-1); }
  /**
   * @brief Sets the value of the pair.
   * 
   * @param other_col MatIdx of the pair column.
   */
  void assign_paired_chain(index other_col) { pairedColumn_ = other_col; }
  /**
   * @brief Unpairs a column.
   */
  void unassign_paired_chain() { pairedColumn_ = -1; };

  /**
   * @brief Assign operator.
   */
  Chain_column_extra_properties& operator=(const Chain_column_extra_properties& other) {
    pivot_ = other.pivot_;
    pairedColumn_ = other.pairedColumn_;
    return *this;
  }
  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_column_extra_properties& col1, Chain_column_extra_properties& col2) {
    std::swap(col1.pivot_, col2.pivot_);
    std::swap(col1.pairedColumn_, col2.pairedColumn_);
  }

 protected:
  id_index get_pivot() const { return pivot_; }
  void swap_pivots(Chain_column_extra_properties& other) { std::swap(pivot_, other.pivot_); }

 private:
  id_index pivot_;      /**< IDIdx index associated to the chain */
  index pairedColumn_;  /**< Represents the (F, G x H) partition of the columns.
                             -1 if in F, MatIdx of image of bijection otherwise. 
                             The pivot of a column in G will always be smaller than the pivot of its image in H. */
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_COLUMN_PROP_H
