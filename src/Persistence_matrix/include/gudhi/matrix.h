/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MASTER_MATRIX_H
#define MASTER_MATRIX_H

#include <type_traits>
#include <vector>
#include <unordered_map>

#include "utilities/utilities.h"
#include "utilities/overlay_id_to_position_index.h"
#include "utilities/overlay_position_to_id_index.h"
#include "options.h"

#include "boundary_matrix/base_swap.h"
#include "boundary_matrix/base_pairing.h"
#include "boundary_matrix/ru_vine_swap.h"
#include "boundary_matrix/ru_rep_cycles.h"
#include "chain_matrix/chain_pairing.h"
#include "chain_matrix/chain_vine_swap.h"
#include "chain_matrix/chain_rep_cycles.h"

#include "boundary_matrix/base_matrix_0000.h"
#include "boundary_matrix/base_matrix_0001.h"
#include "boundary_matrix/base_matrix_0010.h"
#include "boundary_matrix/base_matrix_0011.h"
#include "boundary_matrix/ru_matrix_0000.h"
#include "boundary_matrix/ru_matrix_0001.h"
#include "boundary_matrix/ru_matrix_0010.h"
#include "boundary_matrix/ru_matrix_0011.h"
#include "chain_matrix/chain_matrix_0000.h"
#include "chain_matrix/chain_matrix_0001.h"
#include "chain_matrix/chain_matrix_0010.h"
#include "chain_matrix/chain_matrix_0011.h"
#include "chain_matrix/chain_matrix_1010.h"
#include "chain_matrix/chain_matrix_1011.h"

#include "column_types/list_column.h"
#include "column_types/vector_column.h"
#include "column_types/set_column.h"
#include "column_types/unordered_set_column.h"
#include "column_types/reduced_cell_list_column_with_row.h"
#include "column_types/reduced_cell_set_column_with_row.h"
#include "column_types/z2_heap_column.h"
#include "column_types/z2_list_column.h"
#include "column_types/z2_vector_column.h"
#include "column_types/z2_set_column.h"
#include "column_types/z2_unordered_set_column.h"
#include "column_types/z2_reduced_cell_list_column_with_row.h"
#include "column_types/z2_reduced_cell_set_column_with_row.h"
#include "column_types/cell.h"
#include "column_types/column_pairing.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Options = Default_options<> >
class Matrix
{
public:
	using Option_list = Options;
	using Field_type = typename Options::field_coeff_type;

	struct Dummy_column_pairing{
		Dummy_column_pairing& operator=([[maybe_unused]] Dummy_column_pairing other){return *this;}
		friend void swap([[maybe_unused]] Dummy_column_pairing& d1, [[maybe_unused]] Dummy_column_pairing& d2){}

		Dummy_column_pairing(){}
		Dummy_column_pairing([[maybe_unused]] const Dummy_column_pairing& toCopy){}
		Dummy_column_pairing([[maybe_unused]] Dummy_column_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Column_pairing_option = typename std::conditional<
											Options::is_of_boundary_type,
											Dummy_column_pairing,
											Column_pairing
										>::type;

	using Heap_column_type = Z2_heap_column<Column_pairing_option>;

	using List_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								typename std::conditional<
									Options::has_row_access,
									Z2_reduced_cell_list_column_with_row<Matrix<Options> >,
									Z2_list_column<Column_pairing_option>
								>::type,
								typename std::conditional<
									Options::has_row_access,
									Reduced_cell_list_column_with_row<Matrix<Options> >,
									List_column<Field_type,Column_pairing_option>
								>::type
							>::type;

	using Vector_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								Z2_vector_column<Column_pairing_option>,
								Vector_column<Field_type,Column_pairing_option>
							>::type;

	using Set_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								typename std::conditional<
									Options::has_row_access,
									Z2_reduced_cell_set_column_with_row<Matrix<Options> >,
									Z2_set_column<Column_pairing_option>
								>::type,
								typename std::conditional<
									Options::has_row_access,
									Reduced_cell_set_column_with_row<Matrix<Options> >,
									Set_column<Field_type,Column_pairing_option>
								>::type
							>::type;

	using Unordered_set_column_type = typename std::conditional<
										Field_type::get_characteristic() == 2,
										Z2_unordered_set_column<Column_pairing_option>,
										Unordered_set_column<Field_type,Column_pairing_option>
									>::type;

	using Column_type = typename std::conditional<
								Options::column_type == Column_types::HEAP,
								Heap_column_type,
								typename std::conditional<
									Options::column_type == Column_types::LIST,
									List_column_type,
									typename std::conditional<
										Options::column_type == Column_types::SET,
										Set_column_type,
										typename std::conditional<
											Options::column_type == Column_types::UNORDERED_SET,
											Unordered_set_column_type,
											Vector_column_type
										>::type
									>::type
								>::type
							>::type;

	using barcode_type = typename std::conditional<
										Options::has_removable_columns,
										std::list<Bar>,
										std::vector<Bar>
									>::type;
	using bar_dictionnary_type = typename std::conditional<
											Options::has_removable_columns,
											std::unordered_map<index,typename barcode_type::iterator>,
											std::vector<unsigned int>
										>::type;
	template<typename value_type>
	using dictionnary_type = typename std::conditional<
											Options::has_removable_columns,
											std::unordered_map<unsigned int,value_type>,
											std::vector<value_type>
										>::type;

	using column_container_type = typename std::conditional<
											Options::has_removable_columns,
											std::unordered_map<index,Column_type>,
											std::vector<Column_type>
										>::type;

	using boundary_type = typename std::conditional<
									Field_type::get_characteristic() == 2,
									std::vector<index>,
									std::vector<std::pair<index,Field_type> >
								>::type;
	using boundary_matrix = std::vector<boundary_type>;

	using Row_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								typename std::conditional<
									Options::column_type == Column_types::LIST,
									typename Z2_reduced_cell_list_column_with_row<Matrix<Options> >::Row_type,
									typename Z2_reduced_cell_set_column_with_row<Matrix<Options> >::Row_type
								>::type,
								typename std::conditional<
									Options::column_type == Column_types::LIST,
									typename Reduced_cell_list_column_with_row<Matrix<Options> >::Row_type,
									typename Reduced_cell_set_column_with_row<Matrix<Options> >::Row_type
								>::type
							>::type;

	using Base_matrix_type = typename std::conditional<
											Options::has_removable_columns,
											typename std::conditional<
												Options::has_row_access,
												Base_matrix_with_row_acces_with_removals<Matrix<Options> >,
												Base_matrix_with_removals<Matrix<Options> >
											>::type,
											typename std::conditional<
												Options::has_row_access,
												Base_matrix_with_row_access<Matrix<Options> >,
												Base_matrix<Matrix<Options> >
											>::type
										>::type;

	using RU_matrix_type = typename std::conditional<
											Options::has_removable_columns,
											typename std::conditional<
												Options::has_row_access,
												RU_matrix_with_row_access_with_removals<Matrix<Options> >,
												RU_matrix_with_removals<Matrix<Options> >
											>::type,
											typename std::conditional<
												Options::has_row_access,
												RU_matrix_with_row_access<Matrix<Options> >,
												RU_matrix<Matrix<Options> >
											>::type
										>::type;

	using Chain_matrix_type = typename std::conditional<
											Options::has_removable_columns,
											typename std::conditional<
												Options::has_row_access,
												typename std::conditional<
													Options::has_column_compression,
													Chain_matrix_with_column_compression_with_row_access_with_removals<Matrix<Options> >,
													Chain_matrix_with_row_access_with_removals<Matrix<Options> >
												>::type,
												Chain_matrix_with_removals<Matrix<Options> >
											>::type,
											typename std::conditional<
												Options::has_row_access,
												typename std::conditional<
													Options::has_column_compression,
													Chain_matrix_with_column_compression_with_row_access<Matrix<Options> >,
													Chain_matrix_with_row_access<Matrix<Options> >
												>::type,
												Chain_matrix<Matrix<Options> >
											>::type
										>::type;

	struct Dummy_base_swap{
		Dummy_base_swap& operator=([[maybe_unused]] Dummy_base_swap other){return *this;}
		friend void swap([[maybe_unused]] Dummy_base_swap& d1, [[maybe_unused]] Dummy_base_swap& d2){}

		Dummy_base_swap([[maybe_unused]] column_container_type &matrix){}
		Dummy_base_swap([[maybe_unused]] column_container_type &matrix, [[maybe_unused]] unsigned int numberOfColumns){}
		Dummy_base_swap([[maybe_unused]] const Dummy_base_swap& matrixToCopy){}
		Dummy_base_swap([[maybe_unused]] Dummy_base_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Base_swap_option = typename std::conditional<
											Options::has_vine_update,
											Base_swap<Matrix<Options> >,
											Dummy_base_swap
										>::type;

	struct Dummy_base_pairing{
		Dummy_base_pairing& operator=([[maybe_unused]] Dummy_base_pairing other){return *this;}
		friend void swap([[maybe_unused]] Dummy_base_pairing& d1, [[maybe_unused]] Dummy_base_pairing& d2){}

		Dummy_base_pairing([[maybe_unused]] column_container_type& matrix, [[maybe_unused]] dimension_type& maxDim){}
		Dummy_base_pairing([[maybe_unused]] const Dummy_base_pairing& matrixToCopy){}
		Dummy_base_pairing([[maybe_unused]] Dummy_base_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Base_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update && !Options::can_retrieve_representative_cycles,
											Base_pairing<Matrix<Options> >,
											Dummy_base_pairing
										>::type;

	struct Dummy_ru_pairing{
		Dummy_ru_pairing& operator=([[maybe_unused]] Dummy_ru_pairing other){return *this;}
		friend void swap([[maybe_unused]] Dummy_ru_pairing& d1, [[maybe_unused]] Dummy_ru_pairing& d2){}

		Dummy_ru_pairing(){}
		Dummy_ru_pairing([[maybe_unused]] const Dummy_base_pairing& matrixToCopy){}
		Dummy_ru_pairing([[maybe_unused]] Dummy_base_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											RU_pairing<Matrix<Options> >,
											Dummy_ru_pairing
										>::type;

	struct Dummy_ru_vine_swap{
		Dummy_ru_vine_swap& operator=([[maybe_unused]] Dummy_ru_vine_swap other){return *this;}
		friend void swap([[maybe_unused]] Dummy_ru_vine_swap& d1, [[maybe_unused]] Dummy_ru_vine_swap& d2){}

		Dummy_ru_vine_swap([[maybe_unused]] Base_matrix_type &matrixR, [[maybe_unused]] Base_matrix_type &matrixU, [[maybe_unused]] dictionnary_type<int> &pivotToColumn){}
		Dummy_ru_vine_swap([[maybe_unused]] const Dummy_ru_vine_swap& matrixToCopy){}
		Dummy_ru_vine_swap([[maybe_unused]] Dummy_ru_vine_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											RU_vine_swap<Matrix<Options>>,
											Dummy_ru_vine_swap
										>::type;

	struct Dummy_chain_pairing{
		Dummy_chain_pairing& operator=([[maybe_unused]] Dummy_chain_pairing other){return *this;}
		friend void swap([[maybe_unused]] Dummy_chain_pairing& d1, [[maybe_unused]] Dummy_chain_pairing& d2){}

		Dummy_chain_pairing(){}
		Dummy_chain_pairing([[maybe_unused]] const Dummy_chain_pairing& matrixToCopy){}
		Dummy_chain_pairing([[maybe_unused]] Dummy_chain_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											Chain_pairing<Matrix<Options> >,
											Dummy_chain_pairing
										>::type;

	struct Dummy_chain_vine_swap{
		Dummy_chain_vine_swap& operator=([[maybe_unused]] Dummy_chain_vine_swap other){return *this;}
		friend void swap([[maybe_unused]] Dummy_chain_vine_swap& d1, [[maybe_unused]] Dummy_chain_vine_swap& d2){}

		Dummy_chain_vine_swap([[maybe_unused]] column_container_type& matrix){}
		Dummy_chain_vine_swap([[maybe_unused]] const Dummy_chain_vine_swap& matrixToCopy){}
		Dummy_chain_vine_swap([[maybe_unused]] Dummy_chain_vine_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											Chain_vine_swap<Matrix<Options>>,
											Dummy_chain_vine_swap
										>::type;

	struct Dummy_ru_representative_cycles{
		Dummy_ru_representative_cycles& operator=([[maybe_unused]] Dummy_ru_representative_cycles other){return *this;}
		friend void swap([[maybe_unused]] Dummy_ru_representative_cycles& d1, [[maybe_unused]] Dummy_ru_representative_cycles& d2){}

		Dummy_ru_representative_cycles([[maybe_unused]] Base_matrix_type &matrixR, [[maybe_unused]] Base_matrix_type &matrixU){}
		Dummy_ru_representative_cycles([[maybe_unused]] const Dummy_ru_representative_cycles& matrixToCopy){}
		Dummy_ru_representative_cycles([[maybe_unused]] Dummy_ru_representative_cycles&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_representative_cycles_option = typename std::conditional<
												Options::can_retrieve_representative_cycles,
												RU_representative_cycles<Matrix<Options>>,
												Dummy_ru_representative_cycles
											>::type;

	struct Dummy_chain_representative_cycles{
		Dummy_chain_representative_cycles& operator=([[maybe_unused]] Dummy_chain_representative_cycles other){return *this;}
		friend void swap([[maybe_unused]] Dummy_chain_representative_cycles& d1, [[maybe_unused]] Dummy_chain_representative_cycles& d2){}

		Dummy_chain_representative_cycles([[maybe_unused]] column_container_type& matrix, [[maybe_unused]] dictionnary_type<index>& pivotToPosition){}
		Dummy_chain_representative_cycles([[maybe_unused]] const Dummy_chain_representative_cycles& matrixToCopy){}
		Dummy_chain_representative_cycles([[maybe_unused]] Dummy_chain_representative_cycles&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_representative_cycles_option = typename std::conditional<
													Options::can_retrieve_representative_cycles,
													Chain_representative_cycles<Matrix<Options>>,
													Dummy_chain_representative_cycles
												>::type;

	using cycle_type = std::vector<index>;

	Matrix();
	Matrix(const boundary_matrix& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Matrix(int numberOfColumns);
	Matrix(const Matrix &matrixToCopy);
	Matrix(Matrix&& other) noexcept;

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

	Matrix& operator=(Matrix other);
	friend void swap(Matrix& matrix1, Matrix& matrix2){
		swap(matrix1.matrix_, matrix2.matrix_);
	}

	void print();  //for debug

	//access to optionnal methods
	const barcode_type& get_current_barcode();
	void update_representative_cycles();
	const std::vector<cycle_type>& get_representative_cycles();
	const cycle_type& get_representative_cycle(const Bar& bar);
	bool vine_swap_with_z_eq_1_case(index index);								//by column position with ordered columns
	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//by column id with potentielly unordered columns
	bool vine_swap(index index);												//by column position with ordered columns
	index vine_swap(index columnIndex1, index columnIndex2);					//by column id with potentielly unordered columns

private:
	using matrix_type = typename std::conditional<
							Options::is_of_boundary_type,
							typename std::conditional<
								Options::has_vine_update || Options::can_retrieve_representative_cycles,
								typename std::conditional<
									Options::is_indexed_by_position,
									RU_matrix_type,
									Id_to_position_indexation_overlay<RU_matrix_type,Matrix<Options>>
								>::type,
								typename std::conditional<
									Options::is_indexed_by_position,
									Base_matrix_type,
									Id_to_position_indexation_overlay<Base_matrix_type,Matrix<Options>>
								>::type
							>::type,
							typename std::conditional<
								Options::is_indexed_by_position,
								Position_to_id_indexation_overlay<Chain_matrix_type,Matrix<Options>>,
								Chain_matrix_type
							>::type
						>::type;

	matrix_type matrix_;

	static constexpr void _assert_options();
};

template<class Options>
inline Matrix<Options>::Matrix()
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(const boundary_matrix &boundaries) : matrix_(boundaries)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(int numberOfColumns) : matrix_(numberOfColumns)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(const Matrix &matrixToCopy) : matrix_(matrixToCopy.matrix_)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(Matrix &&other) noexcept : matrix_(std::move(other.matrix_))
{
	_assert_options();
}

template<class Options>
template<class Boundary_type>
inline void Matrix<Options>::insert_boundary(const Boundary_type &boundary)
{
	matrix_.insert_boundary(boundary);
}

template<class Options>
inline typename Matrix<Options>::Column_type &Matrix<Options>::get_column(index columnIndex)
{
	return matrix_.get_column(columnIndex);
}

template<class Options>
inline const typename Matrix<Options>::Column_type &Matrix<Options>::get_column(index columnIndex) const
{
	return matrix_.get_column(columnIndex);
}

template<class Options>
inline typename Matrix<Options>::Row_type &Matrix<Options>::get_row(index rowIndex)
{
	static_assert(Options::has_row_access, "'get_row' is not implemented for the chosen options.");

	return matrix_.get_row(rowIndex);
}

template<class Options>
inline const typename Matrix<Options>::Row_type &Matrix<Options>::get_row(index rowIndex) const
{
	static_assert(Options::has_row_access, "'get_row' is not implemented for the chosen options.");

	return matrix_.get_row(rowIndex);
}

template<class Options>
inline void Matrix<Options>::erase_last()
{
	static_assert(Options::has_removable_columns, "'erase_last' is not implemented for the chosen options.");

	matrix_.erase_last();
}

template<class Options>
inline dimension_type Matrix<Options>::get_max_dimension() const
{
	return matrix_.get_max_dimension();
}

template<class Options>
inline unsigned int Matrix<Options>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Options>
inline dimension_type Matrix<Options>::get_column_dimension(index columnIndex) const
{
	return matrix_.get_column_dimension(columnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(index sourceColumnIndex, index targetColumnIndex)
{
	return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.zero_cell(columnIndex, rowIndex);
}

template<class Options>
inline void Matrix<Options>::zero_column(index columnIndex)
{
	return matrix_.zero_column(columnIndex);
}

template<class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, index rowIndex) const
{
	return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template<class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnIndex);
}

template<class Options>
inline index Matrix<Options>::get_column_with_pivot(index simplexIndex) const
{
	return matrix_.get_column_with_pivot(simplexIndex);
}

template<class Options>
inline index Matrix<Options>::get_pivot(index columnIndex)
{
	return matrix_.get_pivot(columnIndex);
}

template<class Options>
inline Matrix<Options>& Matrix<Options>::operator=(Matrix other)
{
	swap(matrix_, other.matrix_);

	return *this;
}

template<class Options>
inline void Matrix<Options>::print()
{
	return matrix_.print();
}

template<class Options>
inline const typename Matrix<Options>::barcode_type& Matrix<Options>::get_current_barcode()
{
	static_assert(Options::has_column_pairings || Options::has_vine_update, "This method was not enabled.");
	return matrix_.get_current_barcode();
}

template<class Options>
inline void Matrix<Options>::update_representative_cycles()
{
	static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
	matrix_.update_representative_cycles();
}

template<class Options>
inline const std::vector<typename Matrix<Options>::cycle_type>& Matrix<Options>::get_representative_cycles()
{
	static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
	return matrix_.get_representative_cycles();
}

template<class Options>
inline const typename Matrix<Options>::cycle_type& Matrix<Options>::get_representative_cycle(const Bar& bar)
{
	static_assert(Options::can_retrieve_representative_cycles, "This method was not enabled.");
	return matrix_.get_representative_cycle(bar);
}

template<class Options>
inline bool Matrix<Options>::vine_swap_with_z_eq_1_case(index index)
{
	static_assert(Options::has_vine_update && Options::is_indexed_by_position, "This method was not enabled.");
	return matrix_.vine_swap_with_z_eq_1_case(index);
}

template<class Options>
inline index Matrix<Options>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
{
	static_assert(Options::has_vine_update && !Options::is_indexed_by_position, "This method was not enabled.");
	return matrix_.vine_swap_with_z_eq_1_case(columnIndex1, columnIndex2);
}

template<class Options>
inline bool Matrix<Options>::vine_swap(index index)
{
	static_assert(Options::has_vine_update && Options::is_indexed_by_position, "This method was not enabled.");
	return matrix_.vine_swap(index);
}

template<class Options>
inline index Matrix<Options>::vine_swap(index columnIndex1, index columnIndex2)
{
	static_assert(Options::has_vine_update && !Options::is_indexed_by_position, "This method was not enabled.");
	return matrix_.vine_swap(columnIndex1, columnIndex2);
}

template<class Options>
inline constexpr void Matrix<Options>::_assert_options()
{

	static_assert(!Options::has_column_compression || !Options::is_of_boundary_type, "Column compression only exists for chain type matrices.");
	static_assert(!Options::has_vine_update || Options::has_column_pairings, "Vine update requires computation of the barcode (column pairing).");
	static_assert(!Options::can_retrieve_representative_cycles || Options::has_column_pairings, "Representative cycles requires computation of the barcode (column pairing).");
	static_assert(!Options::has_vine_update || Field_type::get_characteristic() == 2, "Vine update currently works only for Z_2 coefficients.");
	static_assert(!Options::has_row_access || Options::column_type == Column_types::LIST || Options::column_type == Column_types::SET, "Row access is currently implemented only for set and list type of columns.");
	static_assert(Options::column_type != Column_types::HEAP || Field_type::get_characteristic() == 2, "Heap column currently works only for Z_2 coefficients.");
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // MASTER_MATRIX_H
