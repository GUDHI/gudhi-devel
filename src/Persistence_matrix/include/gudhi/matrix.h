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

#include <cstddef>
#include <type_traits>
#include <vector>
#include <unordered_map>
#include <map>
#include <functional>
#include <assert.h>
#include <initializer_list>

#include <boost/intrusive/list.hpp>

#include "utilities/utilities.h"
#include "utilities/overlay_id_to_position_index.h"
#include "utilities/overlay_position_to_id_index.h"
#include "options.h"
// #include "column_types/cell.h"
#include "Persistence_matrix/columns/cell_types.h"
// #include "column_types/row_access.h"
#include "Persistence_matrix/columns/row_access.h"

#include "Persistence_matrix/matrix_dimension_holders.h"
#include "Persistence_matrix/base_matrix_row_access.h"
// #include "boundary_matrix/base_swap.h"
// #include "boundary_matrix/base_pairing.h"
#include "Persistence_matrix/base_swap.h"
#include "Persistence_matrix/base_pairing.h"
// #include "boundary_matrix/ru_vine_swap.h"
// #include "boundary_matrix/custom_ru_vine_swap.h"
// #include "boundary_matrix/ru_rep_cycles.h"
#include "Persistence_matrix/ru_pairing.h"
#include "Persistence_matrix/ru_vine_swap.h"
#include "Persistence_matrix/ru_rep_cycles.h"
// #include "chain_matrix/chain_pairing.h"
// #include "chain_matrix/chain_vine_swap.h"
// #include "chain_matrix/custom_chain_vine_swap.h"
// #include "chain_matrix/chain_rep_cycles.h"
#include "Persistence_matrix/chain_pairing.h"
#include "Persistence_matrix/chain_vine_swap.h"
#include "Persistence_matrix/chain_rep_cycles.h"

// #include "base_matrix/base_matrix_0000.h"
// #include "base_matrix/base_matrix_0001.h"
// #include "base_matrix/base_matrix_0010.h"
// #include "base_matrix/base_matrix_0011.h"
#include "Persistence_matrix/base_matrix.h"
// #include "base_matrix/base_matrix_1000.h"
// #include "base_matrix/base_matrix_1010.h"
#include "Persistence_matrix/base_matrix_with_column_compression.h"
// #include "boundary_matrix/boundary_matrix_0000.h"
// #include "boundary_matrix/boundary_matrix_0001.h"
// #include "boundary_matrix/boundary_matrix_0010.h"
// #include "boundary_matrix/boundary_matrix_0011.h"
#include "Persistence_matrix/boundary_matrix.h"
// #include "boundary_matrix/ru_matrix_0000.h"
// #include "boundary_matrix/ru_matrix_0001.h"
// #include "boundary_matrix/ru_matrix_0010.h"
// #include "boundary_matrix/ru_matrix_0011.h"
#include "Persistence_matrix/ru_matrix.h"
// #include "chain_matrix/chain_matrix_0000.h"
// #include "chain_matrix/chain_matrix_0001.h"
// #include "chain_matrix/chain_matrix_0010.h"
// #include "chain_matrix/chain_matrix_0011.h"
#include "Persistence_matrix/chain_matrix.h"

#include "Persistence_matrix/columns/column_dimension_holder.h"
#include "Persistence_matrix/columns/chain_column_extra_properties.h"
#include "Persistence_matrix/columns/intrusive_list_column.h"
#include "Persistence_matrix/columns/intrusive_set_column.h"
#include "Persistence_matrix/columns/list_column.h"
#include "Persistence_matrix/columns/set_column.h"
// #include "gudhi/column_types/list_column.h"
// #include "gudhi/column_types/set_column.h"
#include "gudhi/column_types/unordered_set_column.h"
#include "gudhi/column_types/vector_column.h"
#include "gudhi/column_types/z2_heap_column.h"
// #include "gudhi/column_types/z2_list_column.h"
// #include "gudhi/column_types/z2_set_column.h"
#include "gudhi/column_types/z2_unordered_set_column.h"
#include "gudhi/column_types/z2_vector_column.h"
// #include "gudhi/column_types/intrusive_list_column.h"
// #include "gudhi/column_types/intrusive_set_column.h"
// #include "gudhi/column_types/z2_intrusive_list_column.h"
// #include "gudhi/column_types/z2_intrusive_set_column.h"
// #include "gudhi/column_types/boundary_columns/boundary_list_column.h"
// #include "gudhi/column_types/boundary_columns/boundary_set_column.h"
#include "gudhi/column_types/boundary_columns/boundary_unordered_set_column.h"
#include "gudhi/column_types/boundary_columns/boundary_vector_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_heap_column.h"
// #include "gudhi/column_types/boundary_columns/z2_boundary_list_column.h"
// #include "gudhi/column_types/boundary_columns/z2_boundary_set_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_unordered_set_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_vector_column.h"
// #include "gudhi/column_types/boundary_columns/boundary_intrusive_list_column.h"
// #include "gudhi/column_types/boundary_columns/boundary_intrusive_set_column.h"
// #include "gudhi/column_types/boundary_columns/z2_boundary_intrusive_list_column.h"
// #include "gudhi/column_types/boundary_columns/z2_boundary_intrusive_set_column.h"
// #include "gudhi/column_types/chain_columns/chain_list_column.h"
// #include "gudhi/column_types/chain_columns/chain_set_column.h"
#include "gudhi/column_types/chain_columns/chain_unordered_set_column.h"
#include "gudhi/column_types/chain_columns/chain_vector_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_heap_column.h"
// #include "gudhi/column_types/chain_columns/z2_chain_list_column.h"
// #include "gudhi/column_types/chain_columns/z2_chain_set_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_unordered_set_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_vector_column.h"
// #include "gudhi/column_types/chain_columns/chain_intrusive_list_column.h"
// #include "gudhi/column_types/chain_columns/chain_intrusive_set_column.h"
// #include "gudhi/column_types/chain_columns/z2_chain_intrusive_list_column.h"
// #include "gudhi/column_types/chain_columns/z2_chain_intrusive_set_column.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Options = Default_options<> >
class Matrix
{
public:
	using Option_list = Options;
	using Field_type = typename Options::field_coeff_type;
	using index = typename Options::index_type;
	using dimension_type = typename Options::dimension_type;

	//TODO: move outside
	struct Bar{
		Bar() : dim(-1), birth(-1), death(-1)
		{}

		Bar(dimension_type dim, int birth, int death)
			: dim(dim), birth(birth), death(death)
		{}

		dimension_type dim;
		int birth;
		int death;
	};

	// struct matrix_row_tag;
	// struct matrix_column_tag;

	using base_hook_matrix_row = boost::intrusive::list_base_hook<
				boost::intrusive::tag < matrix_row_tag >
			, boost::intrusive::link_mode < boost::intrusive::auto_unlink >
		>;
	using base_hook_matrix_list_column = boost::intrusive::list_base_hook <
				boost::intrusive::tag < matrix_column_tag >
			, boost::intrusive::link_mode < boost::intrusive::safe_link >
			>;
	using base_hook_matrix_set_column = boost::intrusive::set_base_hook <
				boost::intrusive::tag < matrix_column_tag >
			, boost::intrusive::link_mode < boost::intrusive::safe_link >
			>;
	struct Dummy_row_hook{};
	struct Dummy_column_hook{};

	using row_hook_type = typename std::conditional<
								Options::has_row_access && Options::has_intrusive_rows,
								base_hook_matrix_row,
								Dummy_row_hook
							>::type;
	using column_hook_type = typename std::conditional<
								Options::column_type == Column_types::INTRUSIVE_LIST,
								base_hook_matrix_list_column,
								typename std::conditional<
									Options::column_type == Column_types::INTRUSIVE_SET,
									base_hook_matrix_set_column,
									Dummy_column_hook
								>::type
							>::type;

	using Cell_column_index_option = typename std::conditional<
										Options::has_row_access,
										Cell_column_index<index>,
										Dummy_cell_column_index_mixin
									>::type;
	using Cell_field_element_option = typename std::conditional<
										Options::is_z2,
										Dummy_cell_field_element_mixin,
										Cell_field_element<Field_type>
									>::type;
	using Cell_type = Cell<Matrix<Options> >;

	using cell_rep_type = typename std::conditional<
								Options::is_z2,
								index,
								std::pair<index,Field_type>
							>::type;

	struct CellPairComparator {
		bool operator()(const std::pair<index,Field_type>& p1, const std::pair<index,Field_type>& p2) const
		{
			return p1.first < p2.first;
		};
	};

	template<class Cell_type>
	struct RowCellComp {
		bool operator()(const Cell_type& c1, const Cell_type& c2) const
		{
			return c1.get_column_index() < c2.get_column_index();
		}
	};

	template<class Cell_type>
	struct CellPointerComp {
		bool operator()(const Cell_type* c1, const Cell_type* c2) const
		{
			return *c1 < *c2;
		}
	};

	using Row_type = typename std::conditional<
								Options::has_intrusive_rows,
								boost::intrusive::list<Cell_type, boost::intrusive::constant_time_size<false>, boost::intrusive::base_hook<base_hook_matrix_row> >,
								std::set<Cell_type,RowCellComp<Cell_type> >
							>::type;

	using row_container_type = typename std::conditional<
											Options::has_removable_rows,
											std::map<index,Row_type>,
//											std::unordered_map<index,Row_type>,
											std::vector<Row_type>
										>::type;

	using Row_access_option = typename std::conditional<
											Options::has_row_access,
											Row_access<Matrix<Options> >,
											Dummy_row_access
										>::type;

	using Matrix_row_access_option = typename std::conditional<
											Options::has_row_access,
											Base_matrix_row_access<Row_type, row_container_type, Options::has_removable_rows, index>,
											Dummy_matrix_row_access
										>::type;

	template<typename value_type>
	using dictionnary_type = typename std::conditional<
											Options::has_removable_columns,
											std::unordered_map<unsigned int,value_type>,
											std::vector<value_type>
										>::type;

	static const bool isNonBasic = Options::has_column_pairings || Options::has_vine_update || Options::can_retrieve_representative_cycles;

	using Column_dimension_option = typename std::conditional<
											isNonBasic,
											Column_dimension_holder<Matrix<Options> >,
											Dummy_dimension_holder
										>::type;

	using Chain_column_option = typename std::conditional<
											isNonBasic && !Options::is_of_boundary_type,
											Chain_column_extra_properties<Matrix<Options> >,
											Dummy_chain_properties
										>::type;

	using Heap_column_type = typename std::conditional<
											isNonBasic,
											typename std::conditional<
												Options::is_of_boundary_type,
												Z2_heap_boundary_column,
												Z2_heap_chain_column<dictionnary_type<index>>
											>::type,
											Z2_heap_column
										>::type;

	using List_column_type = List_column<Matrix<Options> >;

	using Vector_column_type = typename std::conditional<
											isNonBasic,
											typename std::conditional<
												Options::is_of_boundary_type,
												typename std::conditional<
													Options::is_z2,
													Z2_vector_boundary_column<Cell_type, Row_access_option>,
													Vector_boundary_column<Field_type, Cell_type, Row_access_option>
												>::type,
												typename std::conditional<
													Options::is_z2,
													Z2_vector_chain_column<dictionnary_type<index>, Cell_type, Row_access_option>,
													Vector_chain_column<dictionnary_type<index>, Field_type, Cell_type, Row_access_option>
												>::type
											>::type,
											typename std::conditional<
												Options::is_z2,
												Z2_vector_column<Cell_type, Row_access_option>,
												Vector_column<Field_type, Cell_type, Row_access_option>
											>::type
										>::type;

	using Set_column_type = Set_column<Matrix<Options> >;

	using Unordered_set_column_type = typename std::conditional<
											isNonBasic,
											typename std::conditional<
												Options::is_of_boundary_type,
												typename std::conditional<
													Options::is_z2,
													Z2_unordered_set_boundary_column<Cell_type, Row_access_option>,
													Unordered_set_boundary_column<Field_type, Cell_type, Row_access_option>
												>::type,
												typename std::conditional<
													Options::is_z2,
													Z2_unordered_set_chain_column<dictionnary_type<index>, Cell_type, Row_access_option>,
													Unordered_set_chain_column<dictionnary_type<index>, Field_type, Cell_type, Row_access_option>
												>::type
											>::type,
											typename std::conditional<
												Options::is_z2,
												Z2_unordered_set_column<Cell_type, Row_access_option>,
												Unordered_set_column<Field_type, Cell_type, Row_access_option>
											>::type
										>::type;

	using Intrusive_list_column_type = Intrusive_list_column<Matrix<Options> >;

	using Intrusive_set_column_type = Intrusive_set_column<Matrix<Options> >;

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
											typename std::conditional<
												Options::column_type == Column_types::VECTOR,
												Vector_column_type,
												typename std::conditional<
													Options::column_type == Column_types::INTRUSIVE_LIST,
													Intrusive_list_column_type,
													Intrusive_set_column_type
												>::type
											>::type
										>::type
									>::type
								>::type
							>::type;

	using column_container_type = typename std::conditional<
											Options::has_removable_columns,
											std::unordered_map<index,Column_type>,
											std::vector<Column_type>
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

	using boundary_type = typename std::conditional<
									Options::is_z2,
									std::initializer_list<index>,
									std::initializer_list<std::pair<index,Field_type> >
								>::type;

	static const bool dimensionIsNeeded = Options::has_column_pairings && !Options::has_vine_update && !Options::can_retrieve_representative_cycles;

	using Matrix_dimension_option = typename std::conditional<
											Options::has_matrix_maximal_dimension_access || dimensionIsNeeded,
											typename std::conditional<
												Options::has_removable_columns,
												Matrix_all_dimension_holder<dimension_type>,
												Matrix_max_dimension_holder<dimension_type>
											>::type,
											Dummy_matrix_dimension_holder
										>::type;

	using Base_matrix_type = typename std::conditional<
												Options::has_column_compression,
												Base_matrix_with_column_compression<Matrix<Options> >,
												Base_matrix<Matrix<Options> >
											>::type;
	using Boundary_matrix_type = Boundary_matrix<Matrix<Options> >;
	using RU_matrix_type = RU_matrix<Matrix<Options> >;
	using Chain_matrix_type = Chain_matrix<Matrix<Options> >;

	template<class Base>
	using Base_swap_option = typename std::conditional<
											Options::has_vine_update || Options::has_column_and_row_swaps,
											Base_swap<Matrix<Options>,Base>,
											Dummy_base_swap
										>::type;
	using Base_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update && !Options::can_retrieve_representative_cycles,
											Base_pairing<Matrix<Options> >,
											Dummy_base_pairing
										>::type;

	using RU_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											RU_pairing<Matrix<Options> >,
											Dummy_ru_pairing
										>::type;
	using RU_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											RU_vine_swap<Matrix<Options> >,
											Dummy_ru_vine_swap
										>::type;
	using RU_representative_cycles_option = typename std::conditional<
												Options::can_retrieve_representative_cycles,
												RU_representative_cycles<Matrix<Options>>,
												Dummy_ru_representative_cycles
											>::type;

	using Chain_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											Chain_pairing<Matrix<Options> >,
											Dummy_chain_pairing
										>::type;
	using Chain_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											Chain_vine_swap<Matrix<Options> >,
											Dummy_chain_vine_swap
										>::type;
	using Chain_representative_cycles_option = typename std::conditional<
													Options::can_retrieve_representative_cycles,
													Chain_representative_cycles<Matrix<Options>>,
													Dummy_chain_representative_cycles
												>::type;

	using cycle_type = std::vector<index>;

	using returned_column_type = typename std::conditional<
									Options::has_column_compression,
									const Column_type,
									Column_type
								>::type;
	using returned_row_type = typename std::conditional<
									Options::has_column_compression && Options::has_row_access,
									const Row_type,
									Row_type
								>::type;
	using insertion_return_type = typename std::conditional<
									Options::is_of_boundary_type || !isNonBasic,
									void,
									std::vector<cell_rep_type>
								>::type;

	Matrix();
	template<class Boundary_type = boundary_type>
	Matrix(const std::vector<Boundary_type>& boundaries);	//simplex indices have to start at 0 and be consecutifs
	Matrix(int numberOfColumns);
	Matrix(std::function<bool(index,index)> birthComparator, 
		   std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	template<class Boundary_type = boundary_type>
	Matrix(const std::vector<Boundary_type>& boundaries,
		   std::function<bool(index,index)> birthComparator, 
		   std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Matrix(int numberOfColumns,
		   std::function<bool(index,index)> birthComparator, 
		   std::function<bool(index,index)> deathComparator = _no_G_death_comparator);
	Matrix(const Matrix &matrixToCopy);
	Matrix(Matrix&& other) noexcept;

	template<class Container_type>
	void insert_column(const Container_type& column);
	template<class Boundary_type = boundary_type>
	insertion_return_type insert_boundary(const Boundary_type& boundary);
	template<class Boundary_type = boundary_type>
	std::vector<cell_rep_type> insert_boundary(index simplexIndex, const Boundary_type& boundary);
	returned_column_type& get_column(index columnIndex);
	const Column_type& get_column(index columnIndex) const;
	//Warning: the get_column_index() function of the row cells returns not
	//the expected type of indices: for boundary matrices, it will returns
	//the simplex number and for chain matrices, it will return the effectiv
	//column index, independently of the indexing chosen in the options.
	returned_row_type& get_row(index rowIndex);
	const Row_type& get_row(index rowIndex) const;
	void remove_maximal_simplex(index columnIndex);		//for boundary or chain matrix
	void erase_column(index columnIndex);	//for basic matrix, where a direct link between column and row is not guaranteed.
	void erase_row(index rowIndex);			// /!\ assumes row is empty

	dimension_type get_max_dimension() const;
	unsigned int get_number_of_columns() const;

	dimension_type get_column_dimension(index columnIndex) const;

	//targetColumn += sourceColumn
	// void add_to(index sourceColumnIndex, index targetColumnIndex);
	// void add_to(int sourceColumnIndex, int targetColumnIndex);
	template<typename Index_type>
	std::enable_if_t<std::is_integral_v<Index_type> > add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex);

	//targetColumn += sourceColumn
	void add_to(Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, index targetColumnIndex);

	//targetColumn = targetColumn * coefficient + sourceColumn
	void add_to(Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex);
	void add_to(const Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Cell_range& sourceColumn, const Field_type& coefficient, index targetColumnIndex);

	//targetColumn += (coefficient * sourceColumn)
	void add_to(const Field_type& coefficient, Column_type& sourceColumn, index targetColumnIndex);
	void add_to(const Field_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex);
	template<class Cell_range>
	void add_to(const Field_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex);

	void zero_cell(index columnIndex, index rowIndex);
	void zero_column(index columnIndex);
	bool is_zero_cell(index columnIndex, index rowIndex);
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

	const Bar& get_bar_from_pivot(index simplexIndex){
		return matrix_.get_bar_from_pivot(simplexIndex);
	}

private:
	using matrix_type = typename std::conditional<
							isNonBasic,
							typename std::conditional<
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
										Boundary_matrix_type,
										Id_to_position_indexation_overlay<Boundary_matrix_type,Matrix<Options>>
									>::type
								>::type,
								typename std::conditional<
									Options::is_indexed_by_position,
									Position_to_id_indexation_overlay<Chain_matrix_type,Matrix<Options>>,
									Chain_matrix_type
								>::type
							>::type,
							Base_matrix_type
						>::type;

	matrix_type matrix_;

	static constexpr void _assert_options();
};

template<class Options>
inline Matrix<Options>::Matrix()
{
	static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings, 
					"When no barcode is recorded with vine swaps, comparaison functions for the columns have to be provided.");
	_assert_options();
}

template<class Options>
template<class Boundary_type>
inline Matrix<Options>::Matrix(const std::vector<Boundary_type> &boundaries) : matrix_(boundaries)
{
	static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings, 
					"When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns have to be provided.");
	assert(Field_type::get_characteristic() != 0 &&
			"Columns cannot be initialized if the coefficient field characteristic is not specified. "
			"Use a compile-time characteristic initialized field type or use another constructor and call coefficient initializer of the chosen field class.");
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(int numberOfColumns) : matrix_(numberOfColumns)
{
	static_assert(Options::is_of_boundary_type || !Options::has_vine_update || Options::has_column_pairings, 
					"When no barcode is recorded with vine swaps for chain matrices, comparaison functions for the columns have to be provided.");
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
{
	static_assert(!Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings, 
					"Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
	_assert_options();
}

template<class Options>
template<class Boundary_type>
inline Matrix<Options>::Matrix(
		const std::vector<Boundary_type> &boundaries,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator)
	: matrix_(boundaries, birthComparator, deathComparator)
{
	static_assert(!Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings, 
					"Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
	assert(Field_type::get_characteristic() != 0 &&
			"Columns cannot be initialized if the coefficient field characteristic is not specified. "
			"Use a compile-time characteristic initialized field type or use another constructor and call coefficient initializer of the chosen field class.");
	_assert_options();
}

template<class Options>
inline Matrix<Options>::Matrix(
		int numberOfColumns,
		std::function<bool(index,index)> birthComparator, 
		std::function<bool(index,index)> deathComparator) 
	: matrix_(numberOfColumns, birthComparator, deathComparator)
{
	static_assert(!Options::is_of_boundary_type && Options::has_vine_update && !Options::has_column_pairings, 
					"Constructor only available for chain matrices when vine swaps are enabled, but barcodes are not recorded.");
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
template<class Container_type>
inline void Matrix<Options>::insert_column(const Container_type& column)
{
	assert(Field_type::get_characteristic() != 0 &&
			"Columns cannot be initialized if the coefficient field characteristic is not specified. "
			"Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field class.");
	static_assert(!isNonBasic,
			"'insert_column' not available for the chosen options. The input has to be in the form of a simplex boundary.");
	matrix_.insert_column(column);
}

template<class Options>
template<class Boundary_type>
inline typename Matrix<Options>::insertion_return_type Matrix<Options>::insert_boundary(const Boundary_type &boundary)
{
	assert(Field_type::get_characteristic() != 0 &&
			"Columns cannot be initialized if the coefficient field characteristic is not specified. "
			"Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field class.");
	return matrix_.insert_boundary(boundary);
}

template<class Options>
template<class Boundary_type>
inline std::vector<typename Matrix<Options>::cell_rep_type> Matrix<Options>::insert_boundary(index simplexIndex, const Boundary_type &boundary)
{
	assert(Field_type::get_characteristic() != 0 &&
			"Columns cannot be initialized if the coefficient field characteristic is not specified. "
			"Use a compile-time characteristic initialized field type or call coefficient initializer of the chosen field class.");
	static_assert(!Options::is_of_boundary_type, "Only enabled for chain matrices.");
	return matrix_.insert_boundary(simplexIndex, boundary);
}

template<class Options>
inline typename Matrix<Options>::returned_column_type &Matrix<Options>::get_column(index columnIndex)
{
	return matrix_.get_column(columnIndex);
}

template<class Options>
inline const typename Matrix<Options>::Column_type &Matrix<Options>::get_column(index columnIndex) const
{
	static_assert(!Options::has_column_compression, "'get_column' cannot be a const method with column compression.");

	return matrix_.get_column(columnIndex);
}

template<class Options>
inline typename Matrix<Options>::returned_row_type &Matrix<Options>::get_row(index rowIndex)
{
	static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

	return matrix_.get_row(rowIndex);
}

template<class Options>
inline const typename Matrix<Options>::Row_type &Matrix<Options>::get_row(index rowIndex) const
{
	static_assert(Options::has_row_access, "'get_row' is not available for the chosen options.");

	return matrix_.get_row(rowIndex);
}

template<class Options>
inline void Matrix<Options>::remove_maximal_simplex(index columnIndex)
{
	static_assert(Options::has_removable_columns && isNonBasic, "'remove_maximal_simplex' is not available for the chosen options.");

	matrix_.remove_maximal_simplex(columnIndex);
}

template<class Options>
inline void Matrix<Options>::erase_column(index columnIndex)
{
	static_assert(Options::has_removable_columns && !isNonBasic, "'erase_column' is not available for the chosen options.");

	matrix_.erase_column(columnIndex);
}

template<class Options>
inline void Matrix<Options>::erase_row(index rowIndex)
{
	static_assert(Options::has_removable_rows && !isNonBasic, "'erase_row' is not available for the chosen options.");

	matrix_.erase_row(rowIndex);
}

template<class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_max_dimension() const
{
	static_assert(isNonBasic, "'get_max_dimension' is not available for the chosen options.");

	return matrix_.get_max_dimension();
}

template<class Options>
inline unsigned int Matrix<Options>::get_number_of_columns() const
{
	return matrix_.get_number_of_columns();
}

template<class Options>
inline typename Matrix<Options>::dimension_type Matrix<Options>::get_column_dimension(index columnIndex) const
{
	static_assert(isNonBasic, "'get_column_dimension' is not available for the chosen options.");

	return matrix_.get_column_dimension(columnIndex);
}

// template<class Options>
// inline void Matrix<Options>::add_to(index sourceColumnIndex, index targetColumnIndex)
// {
// 	return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
// }

// template<class Options>
// inline void Matrix<Options>::add_to(int sourceColumnIndex, int targetColumnIndex)
// {
// 	return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
// }

template<class Options>
template<typename Index_type>
inline std::enable_if_t<std::is_integral_v<Index_type> > Matrix<Options>::add_to(Index_type sourceColumnIndex, Index_type targetColumnIndex)
{
	return matrix_.add_to(sourceColumnIndex, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(Column_type& sourceColumn, index targetColumnIndex)
{
	return matrix_.add_to(sourceColumn, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(const Column_type& sourceColumn, index targetColumnIndex)
{
	static_assert(Options::column_type != Column_types::VECTOR && (!isNonBasic || Options::is_of_boundary_type), "Addition non constant for chain or vector source columns");
	
	return matrix_.add_to(sourceColumn, targetColumnIndex);
}

template<class Options>
template<class Cell_range>
inline void Matrix<Options>::add_to(const Cell_range& sourceColumn, index targetColumnIndex)
{
	static_assert(!isNonBasic, "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain algebraic consistency.");

	return matrix_.add_to(sourceColumn, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");

	return matrix_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(const Column_type& sourceColumn, const Field_type& coefficient, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");
	static_assert(Options::column_type != Column_types::VECTOR && (!isNonBasic || Options::is_of_boundary_type), "Addition non constant for chain or vector source columns");

	return matrix_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template<class Options>
template<class Cell_range>
inline void Matrix<Options>::add_to(const Cell_range& sourceColumn, const Field_type& coefficient, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");
	static_assert(!isNonBasic, "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain algebraic consistency.");

	return matrix_.add_to(sourceColumn, coefficient, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(const Field_type& coefficient, Column_type& sourceColumn, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");

	return matrix_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::add_to(const Field_type& coefficient, const Column_type& sourceColumn, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");
	static_assert(Options::column_type != Column_types::VECTOR && (!isNonBasic || Options::is_of_boundary_type), "Addition non constant for chain or vector source columns");

	return matrix_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template<class Options>
template<class Cell_range>
inline void Matrix<Options>::add_to(const Field_type& coefficient, const Cell_range& sourceColumn, index targetColumnIndex)
{
	static_assert(!Options::is_z2, "addition with a coefficient multiplication is not available for the chosen options.");
	static_assert(!isNonBasic, "For boundary or chain matrices, only additions with columns inside the matrix is allowed to maintain algebraic consistency.");

	return matrix_.add_to(coefficient, sourceColumn, targetColumnIndex);
}

template<class Options>
inline void Matrix<Options>::zero_cell(index columnIndex, index rowIndex)
{
	static_assert(Options::is_of_boundary_type
			&& Options::has_column_pairings
			&& !Options::has_vine_update
			&& !Options::can_retrieve_representative_cycles,
			"'zero_cell' is not available for the chosen options.");

	return matrix_.zero_cell(columnIndex, rowIndex);
}

template<class Options>
inline void Matrix<Options>::zero_column(index columnIndex)
{
	static_assert(Options::is_of_boundary_type
			&& Options::has_column_pairings
			&& !Options::has_vine_update
			&& !Options::can_retrieve_representative_cycles,
			"'zero_column' is not available for the chosen options.");

	return matrix_.zero_column(columnIndex);
}

template<class Options>
inline bool Matrix<Options>::is_zero_cell(index columnIndex, index rowIndex)
{
	return matrix_.is_zero_cell(columnIndex, rowIndex);
}

template<class Options>
inline bool Matrix<Options>::is_zero_column(index columnIndex)
{
	return matrix_.is_zero_column(columnIndex);
}

template<class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_column_with_pivot(index simplexIndex) const
{
	static_assert(isNonBasic && 
			(!Options::is_of_boundary_type
			|| Options::has_vine_update
			|| Options::can_retrieve_representative_cycles),
			"'get_pivot' is not available for the chosen options.");

	return matrix_.get_column_with_pivot(simplexIndex);
}

template<class Options>
inline typename Matrix<Options>::index Matrix<Options>::get_pivot(index columnIndex)
{
	static_assert(!Options::is_of_boundary_type || isNonBasic, "'get_pivot' is not available for the chosen options.");

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
inline typename Matrix<Options>::index Matrix<Options>::vine_swap_with_z_eq_1_case(index columnIndex1, index columnIndex2)
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
inline typename Matrix<Options>::index Matrix<Options>::vine_swap(index columnIndex1, index columnIndex2)
{
	static_assert(Options::has_vine_update && !Options::is_indexed_by_position, "This method was not enabled.");
	return matrix_.vine_swap(columnIndex1, columnIndex2);
}

template<class Options>
inline constexpr void Matrix<Options>::_assert_options()
{
	static_assert(!Options::has_row_access || (Options::column_type != Column_types::SET && Options::column_type != Column_types::UNORDERED_SET) || !Options::has_intrusive_rows, "Intrusive row access is not compatible with column types storing const elements.");
	static_assert(!Options::can_retrieve_representative_cycles || Options::has_column_pairings, "Representative cycles requires computation of the barcode (column pairing).");
	static_assert(!Options::has_vine_update || Options::is_z2, "Vine update currently works only for Z_2 coefficients.");
	static_assert(Options::column_type != Column_types::HEAP || !Options::has_row_access, "Row access is not possible for heap columns.");
	static_assert(Options::column_type != Column_types::HEAP || Options::is_z2, "Heap column works only for Z_2 coefficients.");
	static_assert(Options::column_type != Column_types::HEAP || Options::is_of_boundary_type, "Heap column works only for boundary matrices.");
	static_assert(Options::column_type != Column_types::HEAP || !Options::has_column_compression, "Column compression not compatible with heap columns.");

	//This should be warnings instead, as Options::has_column_compression is just ignored in those cases and don't produces errors.
	static_assert(!Options::has_column_compression || !Options::has_column_pairings, "Column compression not available to compute persistence homology (it would bring no advantages; use it for co-homology instead).");
	static_assert(!Options::has_column_compression || !Options::has_vine_update, "Column compression not available for vineyards.");
	static_assert(!Options::has_column_compression || !Options::can_retrieve_representative_cycles, "Column compression not available to retrieve representative cycles.");
	//Would column removal while column compression be useful? If yes, should erase() remove a single column or the class of columns identical to the input?
	//For a single column, I have an implementation for union-find (not the current one) which allows deleting a single element in constant time, but find becomes log n in worst case.
	//For a column class, we either just ignore the removed class (constant time), but this results in memory residues, or, we have an erase method which is at least linear in the size of the class.
	static_assert(!Options::has_column_compression || !Options::has_removable_columns, "When column compression is used, the removal of columns is not implemented yet.");
}

} //namespace persistence_matrix
} //namespace Gudhi

#endif // MASTER_MATRIX_H
