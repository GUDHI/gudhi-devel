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
	Position_to_id_indexation_overlay(std::vector<Boundary_type>& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Position_to_id_indexation_overlay(int numberOfColumns);
	Position_to_id_indexation_overlay(const Position_to_id_indexation_overlay &matrixToCopy);
	Position_to_id_indexation_overlay(Position_to_id_indexation_overlay&& other) noexcept;

	template<class Boundary_type = boundary_type>
	void insert_boundary(Boundary_type& boundary);
	Column_type& get_column(index columnIndex);
	Row_type& get_row(index rowIndex);
	void erase_last();

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	void add_to(index sourceColumnIndex, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex);
	bool is_zero_column(index columnIndex);

	index get_column_with_pivot(index simplexIndex);
	index get_pivot(index columnIndex);

	Position_to_id_indexation_overlay& operator=(Position_to_id_indexation_overlay other);
	template<class Friend_matrix_type, class Friend_master_matrix_type>
	friend void swap(Position_to_id_indexation_overlay<Friend_matrix_type,Friend_master_matrix_type>& matrix1,
					 Position_to_id_indexation_overlay<Friend_matrix_type,Friend_master_matrix_type>& matrix2);

	void print();  //for debug

	//access to optionnal methods
	const barcode_type& get_current_barcode();
	void update_representative_cycles();
	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);
	void vine_swap_with_z_eq_1_case(index index);								//by column position with ordered columns
	void vine_swap(index position);												//by column position with ordered columns

private:
	Matrix_type matrix_;
	std::vector<index> positionToColumnIndex_;
	index nextIndex_;
};

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay()
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		std::vector<Boundary_type> &boundaries)
	: matrix_(boundaries), positionToColumnIndex_(boundaries.size()), nextIndex_(boundaries.size())
{
	for (unsigned int i = 0; i < boundaries.size(); i++){
		positionToColumnIndex_.at(i) = i;
	}
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		int numberOfColumns)
	: matrix_(numberOfColumns), positionToColumnIndex_(numberOfColumns), nextIndex_(0)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		const Position_to_id_indexation_overlay &matrixToCopy)
	: matrix_(matrixToCopy.matrix_),
	  positionToColumnIndex_(matrixToCopy.positionToColumnIndex_),
	  nextIndex_(matrixToCopy.nextIndex_)
{}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Position_to_id_indexation_overlay(
		Position_to_id_indexation_overlay &&other) noexcept
	: matrix_(std::move(other.matrix_)),
	  positionToColumnIndex_(std::move(other.positionToColumnIndex_)),
	  nextIndex_(std::exchange(other.nextIndex_, 0))
{}

template<class Matrix_type, class Master_matrix_type>
template<class Boundary_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::insert_boundary(Boundary_type &boundary)
{
	matrix_.insert_boundary(boundary);
	if (positionToColumnIndex_.size() <= nextIndex_) {
		positionToColumnIndex_.resize(nextIndex_ + 1);
	}
	positionToColumnIndex_[nextIndex_] = nextIndex_;
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Column_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column(index columnIndex)
{
	return matrix_.get_column(columnIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::Row_type &
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_row(index rowIndex)
{
	return matrix_.get_row(rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::erase_last()
{
	matrix_.erase_last();
	--nextIndex_;
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
	return matrix_.get_column_dimension(matrix_.get_pivot(positionToColumnIndex_.at(columnIndex)));
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.zero_cell(columnIndex, rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::zero_column(index columnIndex)
{
	return matrix_.zero_column(columnIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline bool Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline index Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_column_with_pivot(index simplexIndex)
{
	return matrix_.get_column_with_pivot(simplexIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline index Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(columnIndex);
}

template<class Matrix_type, class Master_matrix_type>
inline Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>&
Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::operator=(Position_to_id_indexation_overlay other)
{
	std::swap(matrix_, other.matrix_);
	std::swap(positionToColumnIndex_, other.positionToColumnIndex_);
	std::swap(nextIndex_, other.nextIndex_);

	return *this;
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::print()
{
	return matrix_.print();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::barcode_type& Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_current_barcode()
{
	return matrix_.get_current_barcode();
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::update_representative_cycles()
{
	matrix_.update_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const std::vector<typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type>& Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycles()
{
	return matrix_.get_representative_cycles();
}

template<class Matrix_type, class Master_matrix_type>
inline const typename Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::cycle_type& Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::get_representative_cycle(const Bar& bar)
{
	return matrix_.get_representative_cycle(bar);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap_with_z_eq_1_case(index index)
{
	matrix_.vine_swap_with_z_eq_1_case(index);
}

template<class Matrix_type, class Master_matrix_type>
inline void Position_to_id_indexation_overlay<Matrix_type,Master_matrix_type>::vine_swap(index position)
{
	index next = matrix_.vine_swap(positionToColumnIndex_.at(position), positionToColumnIndex_.at(position + 1));
	if (next == positionToColumnIndex_.at(position)){
		std::swap(positionToColumnIndex_.at(position),
				  positionToColumnIndex_.at(position + 1));
	}
}

template<class Friend_matrix_type, class Friend_master_matrix_type>
void swap(Position_to_id_indexation_overlay<Friend_matrix_type,Friend_master_matrix_type>& matrix1,
		  Position_to_id_indexation_overlay<Friend_matrix_type,Friend_master_matrix_type>& matrix2)
{
	std::swap(matrix1.matrix_, matrix2.matrix_);
	std::swap(matrix1.positionToColumnIndex_, matrix2.positionToColumnIndex_);
	std::swap(matrix1.nextIndex_, matrix2.nextIndex_);
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // POS_TO_ID_TRANSLATION_H
