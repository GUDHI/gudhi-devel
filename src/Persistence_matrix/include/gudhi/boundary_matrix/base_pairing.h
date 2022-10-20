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

#include "../utilities/utilities.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Master_matrix>
class Base_pairing
{
public:
	using barcode_type = typename Master_matrix::barcode_type;
	using matrix_type = typename Master_matrix::column_container_type;

	Base_pairing(matrix_type& matrix, dimension_type& maxDim);
	Base_pairing(const Base_pairing& matrixToCopy);
	Base_pairing(Base_pairing&& other) noexcept;

	const barcode_type& get_current_barcode();

	Base_pairing& operator=(Base_pairing other);
	friend void swap(Base_pairing& pairing1, Base_pairing& pairing2){
		std::swap(pairing1.matrix_, pairing2.matrix_);
		pairing1.barcode_.swap(pairing2.barcode_);
		std::swap(pairing1.maxDim_, pairing2.maxDim_);
		pairing1.indexToBar_.swap(pairing2.indexToBar_);
		std::swap(pairing1.isReduced_, pairing2.isReduced_);
	}

protected:
	using column_type = typename Master_matrix::Column_type;
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;

	matrix_type* matrix_;
	dimension_type* maxDim_;
	barcode_type barcode_;
	dictionnary_type indexToBar_;
	bool isReduced_;
	static constexpr bool isActive_ = true;

	void _reduce();
};

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(matrix_type &matrix, dimension_type &maxDim)
	: matrix_(&matrix), maxDim_(&maxDim), isReduced_(false)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(const Base_pairing &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  maxDim_(matrixToCopy.maxDim_),
	  barcode_(matrixToCopy.barcode_),
	  indexToBar_(matrixToCopy.indexToBar_),
	  isReduced_(matrixToCopy.isReduced_)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(Base_pairing<Master_matrix> &&other) noexcept
	: matrix_(std::exchange(other.matrix_, nullptr)),
	  maxDim_(std::exchange(other.maxDim_, nullptr)),
	  barcode_(std::move(other.barcode_)),
	  indexToBar_(std::move(other.indexToBar_)),
	  isReduced_(std::move(other.isReduced_))
{}

template<class Master_matrix>
inline void Base_pairing<Master_matrix>::_reduce()
{
	std::unordered_map<index, index> pivotsToColumn;

	for (int d = *maxDim_; d > 0; d--){
		for (unsigned int i = 0; i < matrix_->size(); i++){
			if (!(matrix_->at(i).is_empty()) && matrix_->at(i).get_dimension() == d)
			{
				column_type &curr = matrix_->at(i);
				int pivot = curr.get_pivot();

				while (pivot != -1 && pivotsToColumn.find(pivot) != pivotsToColumn.end()){
					if constexpr (Master_matrix::Field_type::get_characteristic() == 2){
						curr += matrix_->at(pivotsToColumn.at(pivot));
					} else {
						column_type &toadd = matrix_->at(pivotsToColumn.at(pivot));
						typename Master_matrix::Field_type coef = curr.get_pivot_value();
						coef = coef.get_inverse();
						coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));
						curr *= coef;
						curr += toadd;
					}

					pivot = curr.get_pivot();
				}

				if (pivot != -1){
					pivotsToColumn.emplace(pivot, i);
					matrix_->at(pivot).clear();
					barcode_.push_back(Bar(d - 1, pivot, i));
				} else {
					matrix_->at(i).clear();
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
	std::swap(matrix_, other.matrix_);
	barcode_.swap(other.barcode_);
	std::swap(maxDim_, other.maxDim_);
	indexToBar_.swap(other.indexToBar_);
	std::swap(isReduced_, other.isReduced_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // BASE_PAIRING_H
