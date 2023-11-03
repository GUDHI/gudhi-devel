/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BASE_PAIRING_H
#define PM_BASE_PAIRING_H

#include <utility>	//std::swap & std::move
#include <unordered_map>

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_base_pairing{
	Dummy_base_pairing& operator=([[maybe_unused]] Dummy_base_pairing other){return *this;}
	friend void swap([[maybe_unused]] Dummy_base_pairing& d1, [[maybe_unused]] Dummy_base_pairing& d2){}

	Dummy_base_pairing(){}
	Dummy_base_pairing([[maybe_unused]] const Dummy_base_pairing& matrixToCopy){}
	Dummy_base_pairing([[maybe_unused]] Dummy_base_pairing&& other) noexcept{}
};

template<class Master_matrix>
class Base_pairing
{
public:
	using barcode_type = typename Master_matrix::barcode_type;
	using matrix_type = typename Master_matrix::column_container_type;
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;
	using Bar = typename Master_matrix::Bar;

	Base_pairing();
	Base_pairing(const Base_pairing& matrixToCopy);
	Base_pairing(Base_pairing&& other) noexcept;

	const barcode_type& get_current_barcode();

	Base_pairing& operator=(Base_pairing other);
	friend void swap(Base_pairing& pairing1, Base_pairing& pairing2){
		pairing1.barcode_.swap(pairing2.barcode_);
		pairing1.indexToBar_.swap(pairing2.indexToBar_);
		std::swap(pairing1.isReduced_, pairing2.isReduced_);
	}

protected:
	using column_type = typename Master_matrix::Column_type;
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using base_matrix = typename Master_matrix::Boundary_matrix_type;

	barcode_type barcode_;
	dictionnary_type indexToBar_;
	bool isReduced_;

	void _reduce();
	void _remove_maximal(index columnIndex);

	constexpr base_matrix* _matrix() { return static_cast<base_matrix*>(this); }
	constexpr const base_matrix* _matrix() const { return static_cast<const base_matrix*>(this); }
};

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing()
	: isReduced_(false)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(const Base_pairing &matrixToCopy)
	: barcode_(matrixToCopy.barcode_),
	  indexToBar_(matrixToCopy.indexToBar_),
	  isReduced_(matrixToCopy.isReduced_)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(Base_pairing<Master_matrix> &&other) noexcept
	: barcode_(std::move(other.barcode_)),
	  indexToBar_(std::move(other.indexToBar_)),
	  isReduced_(std::move(other.isReduced_))
{}

template<class Master_matrix>
inline const typename Base_pairing<Master_matrix>::barcode_type &
Base_pairing<Master_matrix>::get_current_barcode()
{
	if (!isReduced_) _reduce();
	return barcode_;
}

template<class Master_matrix>
inline void Base_pairing<Master_matrix>::_reduce()
{
	std::unordered_map<index, index> pivotsToColumn;

	for (int d = _matrix()->get_max_dimension(); d > 0; d--){
		for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++){
			auto& curr = _matrix()->get_column(i);
			if (!(curr.is_empty()) && curr.get_dimension() == d)
			{
				int pivot = curr.get_pivot();

				while (pivot != -1 && pivotsToColumn.find(pivot) != pivotsToColumn.end()){
					if constexpr (Master_matrix::Option_list::is_z2){
						curr += _matrix()->get_column(pivotsToColumn.at(pivot));
					} else {
						column_type &toadd = _matrix()->get_column(pivotsToColumn.at(pivot));
						typename Master_matrix::Field_type coef = curr.get_pivot_value();
						coef = coef.get_inverse();
						coef *= (Master_matrix::Field_type::get_characteristic() - static_cast<unsigned int>(toadd.get_pivot_value()));
						curr.multiply_and_add(coef, toadd);
					}

					pivot = curr.get_pivot();
				}

				if (pivot != -1){
					pivotsToColumn.emplace(pivot, i);
					_matrix()->get_column(pivot).clear();
					barcode_.push_back(Bar(d - 1, pivot, i));
					if constexpr (Master_matrix::Option_list::has_removable_columns){
						indexToBar_.emplace(pivot,std::prev(barcode_.end()));
						indexToBar_.emplace(i,std::prev(barcode_.end()));
					}
				} else {
					curr.clear();
					barcode_.push_back(Bar(d, i, -1));
					if constexpr (Master_matrix::Option_list::has_removable_columns){
						indexToBar_.emplace(i,std::prev(barcode_.end()));
					}
				}
			}
		}
	}
	for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++){
		if (_matrix()->get_column(i).get_dimension() == 0 && pivotsToColumn.find(i) == pivotsToColumn.end()){
			barcode_.push_back(Bar(0, i, -1));
			if constexpr (Master_matrix::Option_list::has_removable_columns){
				indexToBar_.emplace(i,std::prev(barcode_.end()));
			}
		}
	}

	isReduced_ = true;
}

template<class Master_matrix>
inline void Base_pairing<Master_matrix>::_remove_maximal(index columnIndex){
	if (isReduced_){
		auto bar = indexToBar_.at(columnIndex);

		if (bar->death == -1) barcode_.erase(bar);
		else bar->death = -1;

		indexToBar_.erase(columnIndex);
	}
}

template<class Master_matrix>
inline Base_pairing<Master_matrix> &Base_pairing<Master_matrix>::operator=(Base_pairing<Master_matrix> other)
{
	barcode_.swap(other.barcode_);
	indexToBar_.swap(other.indexToBar_);
	std::swap(isReduced_, other.isReduced_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_PAIRING_H
