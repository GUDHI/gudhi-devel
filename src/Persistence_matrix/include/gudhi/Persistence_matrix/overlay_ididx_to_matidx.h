/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_ID_TO_POS_TRANSLATION_H
#define PM_ID_TO_POS_TRANSLATION_H

#include <vector>
#include <cassert>
#include <utility>		//std::swap, std::move & std::exchange
#include <algorithm>	//std::transform

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Id_to_index_overlay
{
public:
	using index = typename Master_matrix_type::index;
	using id_index = typename Master_matrix_type::id_index;
	using dimension_type = typename Master_matrix_type::dimension_type;
	using Field_operators = typename Master_matrix_type::Field_operators;
	using Field_element_type = typename Master_matrix_type::element_type;
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;
	using Row_type = typename Master_matrix_type::Row_type;
	using bar_type = typename Master_matrix_type::Bar;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;
	using Cell_constructor = typename Master_matrix_type::Cell_constructor;

	Id_to_index_overlay(Field_operators* operators, Cell_constructor* cellConstructor);
	template<class Boundary_type = boundary_type>
	Id_to_index_overlay(const std::vector<Boundary_type>& boundaries, Field_operators* operators, Cell_constructor* cellConstructor);
	Id_to_index_overlay(unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor);
	Id_to_index_overlay(const Id_to_index_overlay& matrixToCopy, Field_operators* operators = nullptr, Cell_constructor* cellConstructor = nullptr);
	Id_to_index_overlay(Id_to_index_overlay&& other) noexcept;
	~Id_to_index_overlay();

	//boundary: does not update barcode as it needs reduction
	//ru
	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary, dimension_type dim = -1);
	//chain: new simplex = new ID even if the same simplex was already inserted and then removed, ie., an ID cannot come back.
	template<class Boundary_type = boundary_type>
	void insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim = -1);
	//boundary
	//ru: inR = true forced
	//chain
	Column_type& get_column(id_index faceID);
	//boundary
	//ru: inR = true forced
	//chain
	// const Column_type& get_column(id_index faceID) const;
	//get_row(rowIndex) --> simplex ID (=/= faceID)
	//boundary
	//ru: inR = true forced
	//chain
	Row_type& get_row(id_index rowIndex);
	//boundary
	//ru: inR = true forced
	//chain
	// const Row_type& get_row(id_index rowIndex) const;
	//boundary: indirect
	//ru: indirect
	//chain: indirect
	void erase_row(id_index rowIndex);
	//boundary: update barcode if already computed, does not verify if it really was maximal
	//ru
	//chain
	void remove_maximal_face(id_index faceID);
	void remove_maximal_face(id_index faceID, const std::vector<id_index>& columnsToSwap);
	void remove_last();

	//boundary: indirect
	//ru
	//chain: indirect
	dimension_type get_max_dimension() const;
	//boundary
	//ru
	//chain
	index get_number_of_columns() const;
	//boundary
	//ru
	//chain
	dimension_type get_column_dimension(id_index faceID) const;

	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	void add_to(id_index sourceFaceID, id_index targetFaceID);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	void multiply_target_and_add_to(id_index sourceFaceID, const Field_element_type& coefficient, id_index targetFaceID);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: avoid calling with specialized options or make it such that it makes sense for persistence
	void multiply_source_and_add_to(const Field_element_type& coefficient, id_index sourceFaceID, id_index targetFaceID);

	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	void zero_cell(id_index faceID, id_index rowIndex);
	//boundary: avoid calling with pairing option or make it such that it makes sense for persistence
	//ru: inR = true forced, avoid calling with specialized options or make it such that it makes sense for persistence
	void zero_column(id_index faceID);
	//boundary
	//ru: inR = true forced
	//chain
	bool is_zero_cell(id_index faceID, id_index rowIndex) const;
	//boundary
	//ru: inR = true forced
	//chain: just for sanity checks as a valid chain matrix never has an empty column.
	bool is_zero_column(id_index faceID);

	//ru: assumes that pivot exists
	//chain
	id_index get_column_with_pivot(id_index faceIndex) const;	//assumes that pivot exists
	//boundary
	//ru
	//chain
	id_index get_pivot(id_index faceID);

	void reset(Field_operators* operators, Cell_constructor* cellConstructor){
		matrix_.reset(operators, cellConstructor);
		nextIndex_ = 0;
	}

	// void set_operators(Field_operators* operators){ 
	// 	matrix_.set_operators(operators);
	// }

	Id_to_index_overlay& operator=(const Id_to_index_overlay& other);
	friend void swap(Id_to_index_overlay& matrix1,
					 Id_to_index_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		if (Master_matrix_type::Option_list::is_of_boundary_type) 
			std::swap(matrix1.idToIndex_, matrix2.idToIndex_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
	}

	void print();  //for debug

	//access to optionnal methods

	//boundary
	//ru
	//chain
	const barcode_type& get_current_barcode();
	//boundary
	void swap_columns(id_index faceID1, id_index faceID2);
	//boundary
	void swap_rows(index rowIndex1, index rowIndex2);
	//chain
	//ru
	void update_representative_cycles();
	//chain
	//ru
	const std::vector<cycle_type>& get_representative_cycles();
	//chain
	//ru
	const cycle_type& get_representative_cycle(const bar_type& bar);
	//chain: returns id_index which was not modified, ie new i+1
	id_index vine_swap_with_z_eq_1_case(id_index faceID1, id_index faceID2);	//by column id with potentielly unordered columns
	//chain: returns id_index which was not modified, ie new i+1
	id_index vine_swap(id_index faceID1, id_index faceID2);					//by column id with potentielly unordered columns;

private:
	using dictionnary_type = typename Master_matrix_type::template dictionnary_type<index>;

	Matrix_type matrix_;
	dictionnary_type* idToIndex_;
	index nextIndex_;
	
	void _initialize_map(unsigned int size);
	index _id_to_index(id_index id) const;
};

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::Id_to_index_overlay(Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(operators, cellConstructor), idToIndex_(nullptr), nextIndex_(0)
{
	_initialize_map(0);
}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::Id_to_index_overlay(
		const std::vector<Boundary_type> &boundaries, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(boundaries, operators, cellConstructor),
	  idToIndex_(nullptr),
	  nextIndex_(boundaries.size())
{
	_initialize_map(boundaries.size());
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		for (unsigned int i = 0; i < boundaries.size(); i++){
			idToIndex_->operator[](i) = i;
		}
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::Id_to_index_overlay(
		unsigned int numberOfColumns, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(numberOfColumns, operators, cellConstructor), idToIndex_(nullptr), nextIndex_(0)
{
	_initialize_map(numberOfColumns);
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::Id_to_index_overlay(
		const Id_to_index_overlay &matrixToCopy, Field_operators* operators, Cell_constructor* cellConstructor)
	: matrix_(matrixToCopy.matrix_, operators, cellConstructor),
	  idToIndex_(nullptr),
	  nextIndex_(matrixToCopy.nextIndex_)
{
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		idToIndex_ = new dictionnary_type(*matrixToCopy.idToIndex_);
	} else {
		idToIndex_ = &matrix_.pivotToColumnIndex_;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::Id_to_index_overlay(
		Id_to_index_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  idToIndex_(std::exchange(other.idToIndex_, nullptr)),
	  nextIndex_(std::exchange(other.nextIndex_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>::~Id_to_index_overlay()
{
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		if (idToIndex_ != nullptr) delete idToIndex_;
	}
}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type& boundary, dimension_type dim)
{
	matrix_.insert_boundary(boundary, dim);
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			idToIndex_->emplace(nextIndex_, nextIndex_);
		} else {
			if (idToIndex_->size() == nextIndex_) {
				idToIndex_->push_back(nextIndex_);
			} else {
				idToIndex_->operator[](nextIndex_) = nextIndex_;
			}
		}
		++nextIndex_;
	}
}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::insert_boundary(id_index faceIndex, const Boundary_type& boundary, dimension_type dim)
{
	if constexpr (Master_matrix_type::Option_list::has_map_column_container){
		assert(idToIndex_->find(faceIndex) == idToIndex_->end() && "Index for simplex already chosen!");
	} else {
		assert((idToIndex_->size() <= faceIndex || idToIndex_[faceIndex] == -1)  && "Index for simplex already chosen!");
	}
	matrix_.insert_boundary(faceIndex, boundary, dim);
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			idToIndex_->emplace(faceIndex, nextIndex_);
		} else {
			if (idToIndex_->size() <= faceIndex) {
				idToIndex_->resize(faceIndex + 1, -1);
			}
			idToIndex_->operator[](faceIndex) = nextIndex_;
		}
		++nextIndex_;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::Column_type &
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_column(id_index faceID)
{
	return matrix_.get_column(_id_to_index(faceID));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::Row_type &
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_row(id_index rowIndex)
{
	return matrix_.get_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::erase_row(id_index rowIndex)
{
	return matrix_.erase_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::remove_maximal_face(id_index faceID)
{
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		std::vector<id_index> indexToID(nextIndex_);
		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			for (auto& p : *idToIndex_){
				indexToID[p.second] = p.first;
			}
		} else {
			for (id_index i = 0; i < idToIndex_->size(); ++i){
				if (idToIndex_->operator[](i) != -1) indexToID[idToIndex_->operator[](i)] = i;
			}
		}
		--nextIndex_;
		for (index curr = _id_to_index(faceID); curr < nextIndex_; ++curr) {
			matrix_.vine_swap(curr);
			std::swap(idToIndex_->at(indexToID[curr]), idToIndex_->at(indexToID[curr + 1]));
		}
		matrix_.remove_last();
		assert(_id_to_index(faceID) == nextIndex_ && "Indexation problem.");

		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			idToIndex_->erase(faceID);
		} else {
			idToIndex_->operator[](faceID) = -1;
		}
	} else {
		matrix_.remove_maximal_face(faceID);
	}
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::remove_maximal_face(id_index faceID, const std::vector<id_index>& columnsToSwap)
{
	static_assert(!Master_matrix_type::Option_list::is_of_boundary_type, "'remove_maximal_face(id_index,const std::vector<index>&)' is not available for the chosen options.");
	std::vector<index> translatedIndices;
	std::transform(columnsToSwap.cbegin(), columnsToSwap.cend(), std::back_inserter(translatedIndices), 
					[&](id_index id) { return _id_to_index(id); });
	matrix_.remove_maximal_face(faceID, translatedIndices);
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::remove_last()
{
	matrix_.remove_last();
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		--nextIndex_;
		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			auto it = idToIndex_->begin();
			while (it != idToIndex_->end() && it->second != nextIndex_) ++it;
			assert(it != idToIndex_->end());
			idToIndex_->erase(it);
		} else {
			assert(idToIndex_->size() != 0);
			index id = idToIndex_->size() - 1;
			while (idToIndex_->operator[](id) == -1) --id;	//should always stop before reaching -1
			assert(idToIndex_->operator[](id) == nextIndex_);
			idToIndex_->operator[](id) = -1;
		}
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::index Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::dimension_type 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(id_index faceID) const
{
	return matrix_.get_column_dimension(_id_to_index(faceID));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::add_to(id_index sourceFaceID, id_index targetFaceID)
{
	return matrix_.add_to(_id_to_index(sourceFaceID), _id_to_index(targetFaceID));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::multiply_target_and_add_to(id_index sourceFaceID, const Field_element_type& coefficient, id_index targetFaceID)
{
	return matrix_.multiply_target_and_add_to(_id_to_index(sourceFaceID), coefficient, _id_to_index(targetFaceID));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::multiply_source_and_add_to(const Field_element_type& coefficient, id_index sourceFaceID, id_index targetFaceID)
{
	return matrix_.multiply_source_and_add_to(coefficient, _id_to_index(sourceFaceID), _id_to_index(targetFaceID));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::zero_cell(id_index faceID, id_index rowIndex)
{
	return matrix_.zero_cell(_id_to_index(faceID), rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::zero_column(id_index faceID)
{
	return matrix_.zero_column(_id_to_index(faceID));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_index_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(id_index faceID, id_index rowIndex) const
{
	return matrix_.is_zero_cell(_id_to_index(faceID), rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Id_to_index_overlay<Matrix_type,Master_matrix_type>::is_zero_column(id_index faceID)
{
	return matrix_.is_zero_column(_id_to_index(faceID));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::id_index 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(id_index simplexIndex) const
{
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		int pos = matrix_.get_column_with_pivot(simplexIndex);
		unsigned int i = 0;
		while (_id_to_index(i) != pos) ++i;
		return i;
	} else {
		return simplexIndex;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::id_index Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_pivot(id_index faceID)
{
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		return matrix_.get_pivot(_id_to_index(faceID));
	} else {
		return faceID;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Id_to_index_overlay<Matrix_type,Master_matrix_type>&
Id_to_index_overlay<Matrix_type,Master_matrix_type>::operator=(const Id_to_index_overlay& other)
{
	matrix_ = other.matrix_;
	if (Master_matrix_type::Option_list::is_of_boundary_type) 
		idToIndex_ = other.idToIndex_;
	else
		idToIndex_ = &matrix_.pivotToColumnIndex_;
	nextIndex_ = other.nextIndex_;

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::barcode_type& 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_current_barcode()
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::cycle_type>& 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::cycle_type& 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const bar_type& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::swap_columns(id_index faceID1, id_index faceID2)
{
	matrix_.swap_columns(_id_to_index(faceID1), _id_to_index(faceID2));
	std::swap(idToIndex_->at(faceID1), idToIndex_->at(faceID2));
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::swap_rows(index rowIndex1, index rowIndex2)
{
	matrix_.swap_rows(rowIndex1, rowIndex2);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::id_index 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(id_index faceID1, id_index faceID2)
{
	index first = _id_to_index(faceID1);
	index second = _id_to_index(faceID2);
	if (first > second) std::swap(first, second);

	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		assert(second - first == 1 && "The columns to swap are not contiguous.");

		bool change = matrix_.vine_swap_with_z_eq_1_case(first);

		std::swap(idToIndex_->at(faceID1), idToIndex_->at(faceID2));

		if (change){
			return faceID1;
		}
		return faceID2;
	} else {
		return matrix_.vine_swap_with_z_eq_1_case(first, second);
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::id_index 
Id_to_index_overlay<Matrix_type,Master_matrix_type>::vine_swap(id_index faceID1, id_index faceID2)
{
	index first = _id_to_index(faceID1);
	index second = _id_to_index(faceID2);
	if (first > second) std::swap(first, second);

	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		assert(second - first == 1 && "The columns to swap are not contiguous.");
		
		bool change = matrix_.vine_swap(first);

		std::swap(idToIndex_->at(faceID1), idToIndex_->at(faceID2));

		if (change){
			return faceID1;
		}
		return faceID2;
	} else {
		return matrix_.vine_swap(first, second);
	}
}

template<class Matrix_type, class Master_matrix_type>
inline void Id_to_index_overlay<Matrix_type,Master_matrix_type>::_initialize_map([[maybe_unused]] unsigned int size){
	if constexpr (Master_matrix_type::Option_list::is_of_boundary_type){
		if constexpr (Master_matrix_type::Option_list::has_map_column_container){
			idToIndex_ = new dictionnary_type(size);
		} else {
			idToIndex_ = new dictionnary_type(size, -1);
		}
	} else {
		idToIndex_ = &matrix_.pivotToColumnIndex_;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline typename Id_to_index_overlay<Matrix_type,Master_matrix_type>::index Id_to_index_overlay<Matrix_type,Master_matrix_type>::_id_to_index(id_index id) const{
	if constexpr (Master_matrix_type::Option_list::has_map_column_container){
		return idToIndex_->at(id);
	} else {
		return idToIndex_->operator[](id);
	}
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_ID_TO_POS_TRANSLATION_H
