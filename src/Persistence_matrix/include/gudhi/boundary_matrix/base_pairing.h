/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BASE_PAIRING_H
#define BASE_PAIRING_H

#include <utility>

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_pairing
{
public:
	using barcode_type = typename Master_matrix::barcode_type;

	const barcode_type& get_current_barcode();

	Base_pairing& operator=(Base_pairing other);
	template<class Friend_matrix>
	friend void swap(Base_pairing<Friend_matrix>& pairing1,
					 Base_pairing<Friend_matrix>& pairing2);

protected:
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	Base_pairing();
	Base_pairing(Base_pairing &matrixToCopy);
	Base_pairing(Base_pairing&& other) noexcept;

	barcode_type barcode_;
	dictionnary_type indexToBar_;
};

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing()
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(
		Base_pairing &matrixToCopy)
	: barcode_(matrixToCopy.barcode_),
	  indexToBar_(matrixToCopy.indexToBar_)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(Base_pairing<Master_matrix> &&other) noexcept
	: barcode_(std::move(other.barcode_)),
	  indexToBar_(std::move(other.indexToBar_))
{}

template<class Master_matrix>
inline const typename Base_pairing<Master_matrix>::barcode_type &
Base_pairing<Master_matrix>::get_current_barcode()
{
	return barcode_;
}

template<class Master_matrix>
inline Base_pairing<Master_matrix> &Base_pairing<Master_matrix>::operator=(Base_pairing<Master_matrix> other)
{
	std::swap(barcode_, other.barcode_);
	std::swap(indexToBar_, other.indexToBar_);
	return *this;
}

template<class Friend_matrix>
inline void swap(Base_pairing<Friend_matrix>& pairing1,
				 Base_pairing<Friend_matrix>& pairing2)
{
	pairing1.barcode_.swap(pairing2.barcode_);
	pairing1.indexToBar_.swap(pairing2.indexToBar_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_PAIRING_H
