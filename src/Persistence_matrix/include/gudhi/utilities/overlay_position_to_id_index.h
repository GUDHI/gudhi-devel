#ifndef POS_TO_ID_TRANSLATION_H
#define POS_TO_ID_TRANSLATION_H

#include <vector>

#include "utilities.h"  //type definitions

namespace Gudhi {
namespace persistence_matrix {

template<class Matrix_type, class Master_matrix_type>
class Position_to_id_indexation_overlay
{
public:
	using boundary_type = typename Master_matrix_type::boundary_type;
	using Column_type = typename Master_matrix_type::Column_type;
	using Row_type = typename Master_matrix_type::Row_type;
	using barcode_type = typename Master_matrix_type::barcode_type;
	using cycle_type = typename Master_matrix_type::cycle_type;

	Position_to_id_indexation_overlay();
	template<class Boundary_type = boundary_type>
	Position_to_id_indexation_overlay(const std::vector<Boundary_type>& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Position_to_id_indexation_overlay(int numberOfColumns);
	Position_to_id_indexation_overlay(const Position_to_id_indexation_overlay &matrixToCopy);
	Position_to_id_indexation_overlay(Position_to_id_indexation_overlay&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(const Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	Row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex) const;
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex) const;
	index get_pivot(index columnIndex);

	Position_to_id_indexation_overlay& operator=(Position_to_id_indexation_overlay other);
	friend void swap(Position_to_id_indexation_overlay& matrix1,
					 Position_to_id_indexation_overlay& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
		matrix1.columnPositionToID_.swap(matrix2.columnPositionToID_);
		std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
		std::swap(matrix1.nextID_, matrix2.nextID_);
	}

	void print();  //for debug

	//access to optionnal methods
	const barcode_type& get_current_barcode();
	void update_representative_cycles();
	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);
	bool vine_swap_with_z_eq_1_case(index position);							//by column position with ordered columns
	bool vine_swap(index position);												//by column position with ordered columns

private:
	Matrix_type matrix_;
	std::vector<index> columnPositionToID_;
	index nextIndex_;
	index nextID_;
};

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay()
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const std::vector<Boundary_type> &boundaries)
	: matrix_(boundaries),
	  columnPositionToID_(boundaries.size()),
	  nextIndex_(boundaries.size()),
	  nextID_(boundaries.size())
{
	for (unsigned int i = 0; i < boundaries.size(); i++){
		columnPositionToID_[i] = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		int numberOfColumns)
	: matrix_(numberOfColumns), columnPositionToID_(numberOfColumns), nextIndex_(0), nextID_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const Position_to_id_indexation_overlay &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  columnPositionToID_(matrixToCopy.columnPositionToID_),
	  nextIndex_(matrixToCopy.nextIndex_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		Position_to_id_indexation_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  columnPositionToID_(std::move(other.columnPositionToID_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(const Boundary_type &boundary)
{
	matrix_.insert_boundary(boundary);
	if (columnPositionToID_.size() <= nextIndex_) {
		columnPositionToID_.resize(nextIndex_ * 2 + 1);
	}
	columnPositionToID_[nextIndex_++] = nextID_++;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex)
{
	return matrix_.get_column(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex) const
{
	return matrix_.get_column(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex)
{
	return matrix_.get_row(matrix_.get_pivot(columnPositionToID_.at(rowIndex)));
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex) const
{
	return matrix_.get_row(matrix_.get_pivot(columnPositionToID_.at(rowIndex)));
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::erase_last()
{
	matrix_.erase_last();
	--nextIndex_;
	if (columnPositionToID_.size() > nextIndex_ * 2)
		columnPositionToID_.resize(nextIndex_ + 1);
}

template<class Matrix_type, class Master_matrix_type>
inline dimension_type Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Matrix_type, class Master_matrix_type>
inline unsigned int Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Matrix_type, class Master_matrix_type>
inline dimension_type Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_dimension(index columnIndex) const
{
	return matrix_.get_column_dimension(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(columnPositionToID_.at(sourceColumnIndex), columnPositionToID_.at(targetColumnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.zero_cell(columnPositionToID_.at(columnIndex), columnPositionToID_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::zero_column(index columnIndex)
{
	return matrix_.zero_column(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return matrix_.is_zero_cell(columnPositionToID_.at(columnIndex), columnPositionToID_.at(rowIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline index Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(index simplexIndex) const
{
	index id = matrix_.get_column_with_pivot(simplexIndex);
	unsigned int i = 0;
	while (columnPositionToID_[i] != id) ++i;
	return i;
}

template<class Matrix_type, class Master_matrix_type>
inline index Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(columnPositionToID_.at(columnIndex));
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::operator=(Position_to_id_indexation_overlay other)
{
	swap(matrix_, other.matrix_);
	columnPositionToID_.swap(other.columnPositionToID_);
	std::swap(nextIndex_, other.nextIndex_);
	std::swap(nextID_, other.nextID_);

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::barcode_type&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_current_barcode()
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type>&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const Bar& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(index position)
{
	index next = matrix_.vine_swap_with_z_eq_1_case(columnPositionToID_.at(position), columnPositionToID_.at(position + 1));
	if (next == columnPositionToID_.at(position)){
		std::swap(columnPositionToID_.at(position),
				  columnPositionToID_.at(position + 1));
		return true;
	}

	return false;
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type, Master_matrix_type>::vine_swap(index position)
{
	index next = matrix_.vine_swap(columnPositionToID_.at(position), columnPositionToID_.at(position + 1));
	if (next == columnPositionToID_.at(position)){
		std::swap(columnPositionToID_.at(position),
				  columnPositionToID_.at(position + 1));
		return true;
	}

	return false;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // POS_TO_ID_TRANSLATION_H
