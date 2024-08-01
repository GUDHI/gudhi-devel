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
 * @file chain_pairing.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::persistence_matrix::Chain_pairing class and
 * @ref Gudhi::persistence_matrix::Dummy_chain_pairing structure.
 */

#ifndef PM_CHAIN_PAIRING_H
#define PM_CHAIN_PAIRING_H

#include <utility>  //std::move

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inherited instead of @ref Chain_pairing, when the computation of the barcode was not enabled or if the pairing
 * is already managed by the vine update classes.
 */
struct Dummy_chain_pairing
{
  friend void swap([[maybe_unused]] Dummy_chain_pairing& d1, [[maybe_unused]] Dummy_chain_pairing& d2) {}
};

/**
 * @class Chain_pairing chain_pairing.h gudhi/Persistence_matrix/chain_pairing.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref Chain_matrix if the option was enabled.
 * 
 * @tparam Master_matrix An instantiation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class Chain_pairing 
{
 public:
  using Barcode = typename Master_matrix::Barcode;      /**< Barcode type. */
  using Dimension = typename Master_matrix::Dimension;  /**< Dimension value type. */

  /**
   * @brief Default constructor.
   */
  Chain_pairing();
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  Chain_pairing(const Chain_pairing& matrixToCopy);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  Chain_pairing(Chain_pairing&& other) noexcept;

  /**
   * @brief Returns the current barcode which is maintained at any insertion, removal or vine swap.
   * 
   * @return Const reference to the barcode.
   */
  const Barcode& get_current_barcode() const;

  /**
   * @brief Assign operator.
   */
  Chain_pairing& operator=(Chain_pairing other);
  /**
   * @brief Swap operator.
   */
  friend void swap(Chain_pairing& pairing1, Chain_pairing& pairing2) {
    pairing1.barcode_.swap(pairing2.barcode_);
    pairing1.indexToBar_.swap(pairing2.indexToBar_);
    std::swap(pairing1.nextPosition_, pairing2.nextPosition_);
  }

 protected:
  using Dictionary = typename Master_matrix::Bar_dictionary;
  using Pos_index = typename Master_matrix::Pos_index;

  Barcode barcode_;         /**< Bar container. */
  Dictionary indexToBar_;   /**< Map from @ref MatIdx index to bar index. */
  Pos_index nextPosition_;  /**< Next relative position in the filtration. */
};

template <class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing() : nextPosition_(0) 
{}

template <class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing(const Chain_pairing& matrixToCopy)
    : barcode_(matrixToCopy.barcode_),
      indexToBar_(matrixToCopy.indexToBar_),
      nextPosition_(matrixToCopy.nextPosition_) 
{}

template <class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing(Chain_pairing<Master_matrix>&& other) noexcept
    : barcode_(std::move(other.barcode_)),
      indexToBar_(std::move(other.indexToBar_)),
      nextPosition_(std::exchange(other.nextPosition_, 0)) 
{}

template <class Master_matrix>
inline const typename Chain_pairing<Master_matrix>::Barcode& Chain_pairing<Master_matrix>::get_current_barcode()
    const 
{
  return barcode_;
}

template <class Master_matrix>
inline Chain_pairing<Master_matrix>& Chain_pairing<Master_matrix>::operator=(Chain_pairing<Master_matrix> other) 
{
  barcode_.swap(other.barcode_);
  indexToBar_.swap(other.indexToBar_);
  std::swap(nextPosition_, other.nextPosition_);
  return *this;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_CHAIN_PAIRING_H
