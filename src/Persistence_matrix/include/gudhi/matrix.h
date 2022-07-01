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

#include "utilities.h"
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
#include "cell.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Options = Default_options<> >
class Matrix
{
public:
	using Field_type = typename Options::field_coeff_type;
//	using Options::is_double_linked;			//later

//	using Options::is_of_boundary_type;
//	using Options::is_of_chain_type;

//	static const bool is_indexed_by_column_index = false;

//	using Options::is_separated_by_dimension;	//later
//	using Options::is_parallelizable;			//later
//	using Options::has_column_compression;		//later
//	using Options::has_row_access;
//	using Options::has_removable_columns;

//	using Options::has_vine_update;
//	using Options::can_retrieve_representative_cycles;
//	using Options::has_column_pairings;

//	using Cell_type = typename std::conditional<
//							Field_type::get_characteristic() == 2,
//							typename std::conditional<
//								Options::has_row_access,
//								Z2_row_cell,
//								Z2_base_cell
//							>::type,
//							typename std::conditional<
//								Options::has_row_access,
//								Row_cell<Field_type>,
//								Base_cell<Field_type>
//							>::type
//						>::type;

	using Heap_column_type = Z2_heap_column;

	using List_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								typename std::conditional<
									Options::has_row_access,
									Z2_reduced_cell_list_column_with_row,
									Z2_list_column
								>::type,
								typename std::conditional<
									Options::has_row_access,
									Reduced_cell_list_column_with_row<Field_type>,
									List_column<Field_type>
								>::type
							>::type;

	using Vector_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								Z2_vector_column,
								Vector_column<Field_type>
							>::type;

	using Set_column_type = typename std::conditional<
								Field_type::get_characteristic() == 2,
								typename std::conditional<
									Options::has_row_access,
									Z2_reduced_cell_set_column_with_row,
									Z2_set_column
								>::type,
								typename std::conditional<
									Options::has_row_access,
									Reduced_cell_set_column_with_row<Field_type>,
									Set_column<Field_type>
								>::type
							>::type;

	using Unordered_set_column_type = typename std::conditional<
										Field_type::get_characteristic() == 2,
										Z2_unordered_set_column,
										Unordered_set_column<Field_type>
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

	using Base_matrix = typename std::conditional<
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

	using RU_matrix = typename std::conditional<
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

	using Chain_matrix = typename std::conditional<
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
	protected:
		Dummy_base_swap(column_container_type &matrix){}
		Dummy_base_swap(column_container_type &matrix, std::vector<boundary_type>& orderedBoundaries){}
		Dummy_base_swap(column_container_type &matrix, unsigned int numberOfColumns){}
		Dummy_base_swap(Dummy_base_swap& matrixToCopy){}
		Dummy_base_swap(Dummy_base_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Base_swap_option = typename std::conditional<
											Options::has_vine_update,
											Base_swap<Matrix<Options> >,
											Dummy_base_swap
										>::type;

	struct Dummy_base_pairing{
	protected:
		Dummy_base_pairing(){}
		Dummy_base_pairing(Dummy_base_pairing& matrixToCopy){}
		Dummy_base_pairing(Dummy_base_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Base_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update && !Options::can_retrieve_representative_cycles,
											Base_pairing<Matrix<Options> >,
											Dummy_base_pairing
										>::type;

	struct Dummy_ru_pairing{
	protected:
		Dummy_ru_pairing(){}
		Dummy_ru_pairing(Dummy_base_pairing& matrixToCopy){}
		Dummy_ru_pairing(Dummy_base_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											RU_pairing<Matrix<Options> >,
											Dummy_ru_pairing
										>::type;

	struct Dummy_ru_vine_swap{
	protected:
		Dummy_ru_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU){}
		Dummy_ru_vine_swap(Dummy_ru_vine_swap& matrixToCopy){}
		Dummy_ru_vine_swap(Dummy_ru_vine_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											RU_vine_swap<Matrix<Options>>,
											Dummy_ru_vine_swap
										>::type;

	struct Dummy_chain_pairing{
	protected:
		Dummy_chain_pairing(){}
		Dummy_chain_pairing(Dummy_chain_pairing& matrixToCopy){}
		Dummy_chain_pairing(Dummy_chain_pairing&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											Chain_pairing<Matrix<Options> >,
											Dummy_chain_pairing
										>::type;

	struct Dummy_chain_vine_swap{
	protected:
		Dummy_chain_vine_swap(column_container_type& matrix, dictionnary_type<index>& pivotToPosition){}
		Dummy_chain_vine_swap(Dummy_chain_vine_swap& matrixToCopy){}
		Dummy_chain_vine_swap(Dummy_chain_vine_swap&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											Chain_vine_swap<Matrix<Options>>,
											Dummy_chain_vine_swap
										>::type;

	struct Dummy_ru_representative_cycles{
	protected:
		Dummy_ru_representative_cycles(Base_matrix &matrixR, Base_matrix &matrixU){}
		Dummy_ru_representative_cycles(Dummy_ru_representative_cycles& matrixToCopy){}
		Dummy_ru_representative_cycles(Dummy_ru_representative_cycles&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using RU_representative_cycles_option = typename std::conditional<
												Options::can_retrieve_representative_cycles,
												RU_representative_cycles<Matrix<Options>>,
												Dummy_ru_representative_cycles
											>::type;

	struct Dummy_chain_representative_cycles{
	protected:
		Dummy_chain_representative_cycles(column_container_type& matrix, dictionnary_type<index>& pivotToPosition){}
		Dummy_chain_representative_cycles(Dummy_chain_representative_cycles& matrixToCopy){}
		Dummy_chain_representative_cycles(Dummy_chain_representative_cycles&& other) noexcept{}

		static constexpr bool isActive_ = false;
	};

	using Chain_representative_cycles_option = typename std::conditional<
													Options::can_retrieve_representative_cycles,
													Chain_representative_cycles<Matrix<Options>>,
													Dummy_chain_representative_cycles
												>::type;

	Matrix();
	Matrix(boundary_matrix& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Matrix(int numberOfColumns);
	Matrix(Matrix &matrixToCopy);
	Matrix(Matrix&& other) noexcept;

	//	template<class Chain_type>
	//	void insert_boundary(index simplexIndex, Chain_type& boundary);
	//	template<class Chain_type>
	//	void insert_boundary(index simplexIndex, Chain_type& boundary, std::vector<index>& chains_in_F);
	//	dimension_type get_dimension(index simplexIndex) const;
	//	unsigned int get_number_of_simplices() const;

	//	index vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2);	//returns index which was not modified, ie new i+1
	//	index vine_swap(index columnIndex1, index columnIndex2);					//returns index which was not modified, ie new i+1
	//	const barcode_type& get_current_barcode();
	//	void erase_last(index simplexIndex);

	//	void add_to(index sourceColumnIndex, index targetColumnIndex);

	//	index get_column_with_pivot(index simplexIndex);
	//	void get_row(index rowIndex, std::vector<index>& row);
	//	index get_pivot(index columnIndex);

	//	bool is_paired(index columnIndex);
	//	index get_paired_column(index columnIndex);

	//	void print_matrices();  //for debug

	//	Matrix& operator=(Matrix other);
	//	template<class Friend_options>
	//	friend void swap(Matrix& matrix1, Matrix& matrix2);

private:
	using matrix_type = typename std::conditional<
							Options::is_of_boundary_type,
							typename std::conditional<
								Options::has_vine_update || Options::can_retrieve_representative_cycles,
								RU_matrix,
								Base_matrix
							>::type,
							Chain_matrix
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
inline Matrix<Options>::Matrix(boundary_matrix &boundaries) : matrix_(boundaries)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(int numberOfColumns) : matrix_(numberOfColumns)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(Matrix &matrixToCopy) : matrix_(matrixToCopy.matrix_)
{
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(Matrix &&other) noexcept : matrix_(std::move(other.matrix_))
{
	_assert_options();
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

//	using Field_type = typename Options::field_coeff_type;
//	using Column_type = typename Options::column_type;
//	using Options::is_double_linked;			//later

//	using Options::is_of_boundary_type;

//	static const bool is_indexed_by_column_index = false;

//	using Options::is_separated_by_dimension;	//later
//	using Options::is_parallelizable;			//later
//	using Options::has_column_compression;		//later
//	using Options::has_row;
//	using Options::has_removable_columns;

//	using Options::has_vine_update;
//	using Options::can_retrieve_representative_cycles;
	//	using Options::has_column_pairings;
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // MASTER_MATRIX_H
