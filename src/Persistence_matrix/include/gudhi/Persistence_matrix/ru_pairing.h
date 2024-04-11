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
 * @file ru_pairing.h
 * @author Hannah Schreiber
 * @brief Contains the @ref RU_pairing class and @ref Dummy_ru_pairing structure.
 */

#ifndef PM_RU_PAIRING_H
#define PM_RU_PAIRING_H

#include <utility>  //std::move

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Empty structure.
 * Inheritated instead of @ref RU_pairing, when the computation of the barcode was not enabled or if the pairing
 * is already managed by the vine update classes.
 */
struct Dummy_ru_pairing {
  friend void swap([[maybe_unused]] Dummy_ru_pairing& d1, [[maybe_unused]] Dummy_ru_pairing& d2) {}
};

/**
 * @class RU_pairing ru_pairing.h gudhi/Persistence_matrix/ru_pairing.h
 * @ingroup persistence_matrix
 *
 * @brief Class managing the barcode for @ref RU_matrix if the option was enabled.
 * 
 * @tparam Master_matrix An instanciation of @ref Matrix from which all types and options are deduced.
 */
template <class Master_matrix>
class RU_pairing 
{
 public:
  using barcode_type = typename Master_matrix::barcode_type;  /**< Barcode type. */

  /**
   * @brief Default constructor.
   */
  RU_pairing();
  /**
   * @brief Copy constructor.
   * 
   * @param matrixToCopy Matrix to copy.
   */
  RU_pairing(const RU_pairing& matrixToCopy);
  /**
   * @brief Move constructor.
   * 
   * @param other Matrix to move.
   */
  RU_pairing(RU_pairing&& other) noexcept;

  /**
   * @brief Returns the current barcode which is maintained at any insertion, removal or vine swap.
   * 
   * @return Const reference to the barcode.
   */
  const barcode_type& get_current_barcode() const;

  /**
   * @brief Assign operator.
   */
  RU_pairing& operator=(RU_pairing other);
  /**
   * @brief Swap operator.
   */
  friend void swap(RU_pairing& pairing1, RU_pairing& pairing2) {
    pairing1.barcode_.swap(pairing2.barcode_);
    pairing1.indexToBar_.swap(pairing2.indexToBar_);
  }

 protected:
  using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

  barcode_type barcode_;        /**< Bar container. */
  dictionnary_type indexToBar_; /**< Map from @ref MatIdx index to bar index. */
};

template <class Master_matrix>
inline RU_pairing<Master_matrix>::RU_pairing() 
{}

template <class Master_matrix>
inline RU_pairing<Master_matrix>::RU_pairing(const RU_pairing& matrixToCopy)
    : barcode_(matrixToCopy.barcode_), indexToBar_(matrixToCopy.indexToBar_) 
{}

template <class Master_matrix>
inline RU_pairing<Master_matrix>::RU_pairing(RU_pairing<Master_matrix>&& other) noexcept
    : barcode_(std::move(other.barcode_)), indexToBar_(std::move(other.indexToBar_)) 
{}

template <class Master_matrix>
inline const typename RU_pairing<Master_matrix>::barcode_type& RU_pairing<Master_matrix>::get_current_barcode() const 
{
  return barcode_;
}

template <class Master_matrix>
inline RU_pairing<Master_matrix>& RU_pairing<Master_matrix>::operator=(RU_pairing<Master_matrix> other) 
{
  barcode_.swap(other.barcode_);
  indexToBar_.swap(other.indexToBar_);
  return *this;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_RU_PAIRING_H
