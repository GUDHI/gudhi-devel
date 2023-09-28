/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022-23 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_BASE_SWAP_H
#define PM_BASE_SWAP_H

#include <utility>	//std::swap, std::move & std::exchange
#include <vector>

namespace Gudhi {
namespace persistence_matrix {

struct Dummy_base_swap{
	Dummy_base_swap& operator=([[maybe_unused]] Dummy_base_swap other){return *this;}
	friend void swap([[maybe_unused]] Dummy_base_swap& d1, [[maybe_unused]] Dummy_base_swap& d2){}

	Dummy_base_swap(){}
	Dummy_base_swap([[maybe_unused]] unsigned int numberOfColumns){}
	Dummy_base_swap([[maybe_unused]] const Dummy_base_swap& matrixToCopy){}
	Dummy_base_swap([[maybe_unused]] Dummy_base_swap&& other) noexcept{}
};

template<class Master_matrix, class Base_matrix>
class Base_swap
{
public:
	using matrix_type = typename Master_matrix::column_container_type;
	using index = typename Master_matrix::index;

	Base_swap();
	Base_swap(unsigned int numberOfColumns);
	Base_swap(const Base_swap& matrixToCopy);
	Base_swap(Base_swap&& other) noexcept;

	void swap_columns(index columnIndex1, index columnIndex2);
	void swap_rows(index rowIndex1, index rowIndex2);
	void swap_at_indices(index index1, index index2);

	Base_swap& operator=(Base_swap other);
	friend void swap(Base_swap& base1, Base_swap& base2){
		base1.indexToRow_.swap(base2.indexToRow_);
		base1.rowToIndex_.swap(base2.rowToIndex_);
		std::swap(base1.rowSwapped_, base2.rowSwapped_);
	}

protected:
	using index_dictionnary_type = typename Master_matrix::template dictionnary_type<unsigned int>;
	using row_dictionnary_type = typename Master_matrix::template dictionnary_type<index>;

	index_dictionnary_type indexToRow_;
	row_dictionnary_type rowToIndex_;
	bool rowSwapped_;

	void _orderRows();

	constexpr Base_matrix* _matrix() { return static_cast<Base_matrix*>(this); }
	constexpr const Base_matrix* _matrix() const { return static_cast<const Base_matrix*>(this); }
};

template<class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix,Base_matrix>::Base_swap()
	: rowSwapped_(false)
{}

template<class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix,Base_matrix>::Base_swap(unsigned int numberOfColumns)
	: indexToRow_(numberOfColumns),
	  rowToIndex_(numberOfColumns),
	  rowSwapped_(false)
{
	for (unsigned int i = 0; i < numberOfColumns; i++){
		indexToRow_[i] = i;
		rowToIndex_[i] = i;
	}
}

template<class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix,Base_matrix>::Base_swap(const Base_swap<Master_matrix,Base_matrix>& matrixToCopy)
	: indexToRow_(matrixToCopy.indexToRow_),
	  rowToIndex_(matrixToCopy.rowToIndex_),
	  rowSwapped_(matrixToCopy.rowSwapped_)
{}

template<class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix,Base_matrix>::Base_swap(Base_swap<Master_matrix,Base_matrix> &&other) noexcept
	: indexToRow_(std::move(other.indexToRow_)),
	  rowToIndex_(std::move(other.rowToIndex_)),
	  rowSwapped_(std::exchange(other.rowSwapped_, 0))
{}

template<class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix,Base_matrix>::swap_columns(index columnIndex1, index columnIndex2)
{
	swap(_matrix()->matrix_.at(columnIndex1), _matrix()->matrix_.at(columnIndex2));
}

template<class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix,Base_matrix>::swap_rows(index rowIndex1, index rowIndex2)
{
	rowSwapped_ = true;
	std::swap(rowToIndex_[indexToRow_[rowIndex1]], rowToIndex_[indexToRow_[rowIndex2]]);
	std::swap(indexToRow_[rowIndex1], indexToRow_[rowIndex2]);
}

template<class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix,Base_matrix>::swap_at_indices(index index1, index index2)
{
	swap_columns(index1, index2);
	swap_rows(index1, index2);
}

template<class Master_matrix, class Base_matrix>
inline Base_swap<Master_matrix,Base_matrix> &Base_swap<Master_matrix,Base_matrix>::operator=(Base_swap other)
{
	indexToRow_.swap(other.indexToRow_);
	rowToIndex_.swap(other.rowToIndex_);
	std::swap(rowSwapped_, other.rowSwapped_);
	return *this;
}

template<class Master_matrix, class Base_matrix>
inline void Base_swap<Master_matrix,Base_matrix>::_orderRows()
{
	for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++){
		_matrix()->matrix_.at(i).reorder(rowToIndex_);
	}
	for (unsigned int i = 0; i < _matrix()->get_number_of_columns(); i++){
		indexToRow_[i] = i;
		rowToIndex_[i] = i;
	}
	rowSwapped_ = false;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // PM_BASE_SWAP_H
