/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CHAIN_PAIRING_H
#define CHAIN_PAIRING_H

#include <utility>

#include "../utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Chain_pairing
{
public:
	using barcode_type = typename Master_matrix::barcode_type;

	const barcode_type& get_current_barcode();

	Chain_pairing& operator=(Chain_pairing other);
	template<class Friend_matrix>
	friend void swap(Chain_pairing<Friend_matrix>& pairing1,
					 Chain_pairing<Friend_matrix>& pairing2);

	Chain_pairing();
	Chain_pairing(const Chain_pairing &matrixToCopy);
	Chain_pairing(Chain_pairing&& other) noexcept;

protected:
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	barcode_type barcode_;
	dictionnary_type indexToBar_;
	static constexpr bool isActive_ = true;

	dimension_type _get_dimension(index simplexIndex) const;		//to move
};

template<class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing()
{}

template<class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing(const Chain_pairing &matrixToCopy)
	: barcode_(matrixToCopy.barcode_),
	  indexToBar_(matrixToCopy.indexToBar_)
{}

template<class Master_matrix>
inline Chain_pairing<Master_matrix>::Chain_pairing(Chain_pairing<Master_matrix> &&other) noexcept
	: barcode_(std::move(other.barcode_)),
	  indexToBar_(std::move(other.indexToBar_))
{}

template<class Master_matrix>
inline const typename Chain_pairing<Master_matrix>::barcode_type &
Chain_pairing<Master_matrix>::get_current_barcode()
{
	return barcode_;
}

template<class Master_matrix>
inline Chain_pairing<Master_matrix> &Chain_pairing<Master_matrix>::operator=(Chain_pairing<Master_matrix> other)
{
	std::swap(barcode_, other.barcode_);
	std::swap(indexToBar_, other.indexToBar_);
	return *this;
}

template<class Master_matrix>
inline dimension_type Chain_pairing<Master_matrix>::_get_dimension(index pivot) const
{
	if (indexToBar_.at(pivot)->birth == static_cast<int>(pivot))
		return indexToBar_.at(pivot)->dim;
	else
		return indexToBar_.at(pivot)->dim + 1;
}

template<class Friend_matrix>
inline void swap(Chain_pairing<Friend_matrix>& pairing1,
				 Chain_pairing<Friend_matrix>& pairing2)
{
	pairing1.barcode_.swap(pairing2.barcode_);
	pairing1.indexToBar_.swap(pairing2.indexToBar_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // CHAIN_PAIRING_H
