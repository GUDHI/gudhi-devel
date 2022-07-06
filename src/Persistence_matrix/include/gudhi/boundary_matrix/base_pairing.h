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
#include <unordered_map>

#include "../utilities.h"

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
	using matrix_type = typename Master_matrix::column_container_type;
	using column_type = typename Master_matrix::Column_type;
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	Base_pairing(matrix_type& matrix, dimension_type& maxDim);
	Base_pairing(Base_pairing& matrixToCopy);
	Base_pairing(Base_pairing&& other) noexcept;

	matrix_type& matrix_;
	dimension_type& maxDim_;
	barcode_type barcode_;
	dictionnary_type indexToBar_;
	bool isReduced_;
	static constexpr bool isActive_ = true;

	void _reduce();
};

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(matrix_type &matrix, dimension_type &maxDim)
	: matrix_(matrix), maxDim_(maxDim), isReduced_(false)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(
		Base_pairing &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  barcode_(matrixToCopy.barcode_),
	  maxDim_(matrixToCopy.maxDim_),
	  isReduced_(matrixToCopy.isReduced_)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(Base_pairing<Master_matrix> &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  barcode_(std::move(other.barcode_)),
	  maxDim_(std::move(other.maxDim_)),
	  isReduced_(std::move(other.isReduced_))
{}

template<class Master_matrix>
inline void Base_pairing<Master_matrix>::_reduce()
{
	std::unordered_map<index, index> pivotsToColumn;

	for (int d = maxDim_; d > 0; d--){
		for (unsigned int i = 0; i < matrix_.size(); i++){
			if (!(matrix_.at(i).is_empty()) && matrix_.at(i).get_dimension() == d)
			{
				column_type &curr = matrix_.at(i);
				int pivot = curr.get_pivot();

				while (pivot != -1 && pivotsToColumn.find(pivot) != pivotsToColumn.end()){
					if constexpr (Master_matrix::Field_type::get_characteristic() == 2){
						curr += matrix_.at(pivotsToColumn.at(pivot));
					} else {
						column_type &toadd = matrix_.at(pivotsToColumn.at(pivot));
						typename Master_matrix::Field_type coef = curr.get_pivot_value();
						coef = coef.get_inverse();
						coef *= (Master_matrix::Field_type::get_characteristic() - toadd.get_pivot_value());
						curr *= coef;
						curr += toadd;
					}

					pivot = curr.get_pivot();
				}

				if (pivot != -1){
					pivotsToColumn.emplace(pivot, i);
					matrix_.at(pivot).clear();
					barcode_.push_back(Bar(d - 1, pivot, i));
				} else {
					matrix_.at(i).clear();
					barcode_.push_back(Bar(d, i, -1));
				}
			}
		}
	}

	isReduced_ = true;
}

template<class Master_matrix>
inline const typename Base_pairing<Master_matrix>::barcode_type &
Base_pairing<Master_matrix>::get_current_barcode()
{
	if (!isReduced_) _reduce();
	return barcode_;
}

template<class Master_matrix>
inline Base_pairing<Master_matrix> &Base_pairing<Master_matrix>::operator=(Base_pairing<Master_matrix> other)
{
	std::swap(barcode_, other.barcode_);
	std::swap(matrix_, other.matrix_);
	std::swap(maxDim_, other.maxDim_);
	return *this;
}

template<class Friend_matrix>
inline void swap(Base_pairing<Friend_matrix>& pairing1,
				 Base_pairing<Friend_matrix>& pairing2)
{
	std::swap(pairing1.matrix_, pairing2.matrix_);
	pairing1.barcode_.swap(pairing2.barcode_);
	std::swap(pairing1.maxDim_, pairing2.maxDim_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_PAIRING_H
