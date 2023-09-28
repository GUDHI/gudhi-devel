/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_CHAIN_PAIRING_H
#define PM_CHAIN_PAIRING_H

#include <utility>	//std::move

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_chain_pairing{
	friend void swap([[maybe_unused]] Dummy_chain_pairing& d1, [[maybe_unused]] Dummy_chain_pairing& d2){}

	Dummy_chain_pairing(){}
};

template<class Master_matrix>
class Chain_pairing
{
public:
	using barcode_type = typename Master_matrix::barcode_type;
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;

	Chain_pairing();
	Chain_pairing(const Chain_pairing &matrixToCopy);
	Chain_pairing(Chain_pairing&& other) noexcept;

	const barcode_type& get_current_barcode() const;

	Chain_pairing& operator=(Chain_pairing other);
	friend void swap(Chain_pairing& pairing1,
					 Chain_pairing& pairing2){
		pairing1.barcode_.swap(pairing2.barcode_);
		pairing1.indexToBar_.swap(pairing2.indexToBar_);
	}

protected:
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	barcode_type barcode_;
	dictionnary_type indexToBar_;

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
Chain_pairing<Master_matrix>::get_current_barcode() const
{
	return barcode_;
}

template<class Master_matrix>
inline Chain_pairing<Master_matrix> &Chain_pairing<Master_matrix>::operator=(Chain_pairing<Master_matrix> other)
{
	barcode_.swap(other.barcode_);
	indexToBar_.swap(other.indexToBar_);
	return *this;
}

template<class Master_matrix>
inline typename Chain_pairing<Master_matrix>::dimension_type Chain_pairing<Master_matrix>::_get_dimension(index pivot) const
{
	if (indexToBar_.at(pivot)->birth == static_cast<int>(pivot))
		return indexToBar_.at(pivot)->dim;
	else
		return indexToBar_.at(pivot)->dim + 1;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_CHAIN_PAIRING_H
