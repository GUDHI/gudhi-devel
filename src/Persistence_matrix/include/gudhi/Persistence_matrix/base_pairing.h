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
#include <algorithm>
#include <vector>

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
	using Bar = typename Master_matrix::Bar;
	using barcode_type = typename Master_matrix::barcode_type;
	using matrix_type = typename Master_matrix::column_container_type;
	using index = typename Master_matrix::index;
	using dimension_type = typename Master_matrix::dimension_type;

	Base_pairing();
	Base_pairing(const Base_pairing& matrixToCopy);
	Base_pairing(Base_pairing&& other) noexcept;

	const barcode_type& get_current_barcode();

	Base_pairing& operator=(Base_pairing other);
	friend void swap(Base_pairing& pairing1, Base_pairing& pairing2){
		pairing1.barcode_.swap(pairing2.barcode_);
		pairing1.deathToBar_.swap(pairing2.deathToBar_);
		std::swap(pairing1.isReduced_, pairing2.isReduced_);
	}

protected:
	using pos_index = typename Master_matrix::pos_index;
	using dictionnary_type = typename Master_matrix::bar_dictionnary_type;
	using base_matrix = typename Master_matrix::Boundary_matrix_type;

	barcode_type barcode_;
	dictionnary_type deathToBar_;	//records deaths only
	bool isReduced_;

	void _reduce();
	void _remove_last(pos_index columnIndex);

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
	  deathToBar_(matrixToCopy.deathToBar_),
	  isReduced_(matrixToCopy.isReduced_)
{}

template<class Master_matrix>
inline Base_pairing<Master_matrix>::Base_pairing(Base_pairing<Master_matrix> &&other) noexcept
	: barcode_(std::move(other.barcode_)),
	  deathToBar_(std::move(other.deathToBar_)),
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
	using id_index = typename Master_matrix::index;
	std::unordered_map<id_index, index> pivotsToColumn(_matrix()->get_number_of_columns());

	auto dim = _matrix()->get_max_dimension();
	std::vector<std::vector<index> > columnsByDim(dim + 1);
	for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++){
		columnsByDim[dim - _matrix()->get_column_dimension(i)].push_back(i);
	}

	for (auto cols : columnsByDim){
		for (auto i : cols){
			auto& curr = _matrix()->get_column(i);
			if (curr.is_empty()){
				if (pivotsToColumn.find(i) == pivotsToColumn.end()){
					barcode_.emplace_back(dim, i, -1);
				}
			} else {
				id_index pivot = curr.get_pivot();

				while (pivot != static_cast<id_index>(-1) && pivotsToColumn.find(pivot) != pivotsToColumn.end()){
					if constexpr (Master_matrix::Option_list::is_z2){
						curr += _matrix()->get_column(pivotsToColumn.at(pivot));
					} else {
						auto &toadd = _matrix()->get_column(pivotsToColumn.at(pivot));
						typename Master_matrix::element_type coef = curr.get_pivot_value();
						coef = _matrix()->operators_->get_inverse(coef);
						coef = _matrix()->operators_->multiply(coef, _matrix()->operators_->get_characteristic() - toadd.get_pivot_value());
						curr.multiply_and_add(coef, toadd);
					}

					pivot = curr.get_pivot();
				}

				if (pivot != static_cast<id_index>(-1)){
					pivotsToColumn.emplace(pivot, i);
					_matrix()->get_column(pivot).clear();
					barcode_.emplace_back(dim - 1, pivot, i);
				} else {
					curr.clear();
					barcode_.emplace_back(dim, i, -1);
				}
			}
		}
		--dim;
	}

	if constexpr (Master_matrix::Option_list::has_removable_columns){
		//sort barcode by birth such that a removal is trivial
		std::sort(barcode_.begin(), barcode_.end(), [](const Bar& b1, const Bar& b2){ return b1.birth < b2.birth; });
		//map can only be constructed once barcode is sorted
		for (index i = 0; i < barcode_.size(); ++i){
			auto d = barcode_[i].death;
			if (d != static_cast<pos_index>(-1)){
				deathToBar_.emplace(d, i);
			}
		}
	}

	isReduced_ = true;
}

template<class Master_matrix>
inline void Base_pairing<Master_matrix>::_remove_last(pos_index columnIndex){
	static_assert(Master_matrix::Option_list::has_removable_columns, "remove_last not available.");

	if (isReduced_){
		auto it = deathToBar_.find(columnIndex);

		if (it == deathToBar_.end()) {	//birth
			barcode_.pop_back();	//sorted by birth and columnIndex has to be the heighest one
		} else {						//death
			barcode_[it->second].death = -1;
			deathToBar_.erase(it);
		};
	}
}

template<class Master_matrix>
inline Base_pairing<Master_matrix> &Base_pairing<Master_matrix>::operator=(Base_pairing<Master_matrix> other)
{
	barcode_.swap(other.barcode_);
	deathToBar_.swap(other.deathToBar_);
	std::swap(isReduced_, other.isReduced_);
	return *this;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_PAIRING_H
