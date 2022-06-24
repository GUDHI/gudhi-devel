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
#include "chain_matrix/chain_pairing.h"
#include "chain_matrix/chain_vine_swap.h"

#include "boundary_matrix/base_matrix_0000.h"
#include "boundary_matrix/base_matrix_0001.h"
#include "boundary_matrix/base_matrix_0010.h"
#include "boundary_matrix/base_matrix_0011.h"
#include "boundary_matrix/ru_matrix_0000.h"
#include "boundary_matrix/ru_matrix_0001.h"
#include "boundary_matrix/ru_matrix_0010.h"
#include "boundary_matrix/ru_matrix_0011.h"

namespace Gudhi {
namespace persistence_matrix {

template<class Options = Default_options<> >
class Matrix
{
public:
	using Field_type = typename Options::field_coeff_type;
	using Column_type = typename Options::column_type;
//	using Options::is_double_linked;

//	using Options::is_of_boundary_type;
//	using Options::is_of_chain_type;

//	static const bool is_indexed_by_column_index = false;

//	using Options::is_separated_by_dimension;	//later
//	using Options::is_parallelizable;			//later
//	using Options::has_column_compression;		//later
//	using Options::has_row;
//	using Options::has_removable_columns;

//	using Options::has_vine_update;
//	using Options::can_retrieve_representative_cycles;
//	using Options::has_column_pairings;

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
											std::unordered_map<index,typename barcode_type::iterator>,
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

	struct Dummy_base_swap{
	protected:
		Dummy_base_swap(column_container_type &matrix){}
		Dummy_base_swap(column_container_type &matrix, std::vector<boundary_type>& orderedBoundaries){}
		Dummy_base_swap(column_container_type &matrix, unsigned int numberOfColumns){}
		Dummy_base_swap(Dummy_base_swap& matrixToCopy){}
		Dummy_base_swap(Dummy_base_swap&& other) noexcept{}
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
	};

	using Base_pairing_option = typename std::conditional<
											Options::has_column_pairings && !Options::has_vine_update,
											Base_pairing<Matrix<Options> >,
											Dummy_base_pairing
										>::type;

	struct Dummy_ru_vine_swap{
	protected:
		Dummy_ru_vine_swap(Base_matrix &matrixR, Base_matrix &matrixU){}
		Dummy_ru_vine_swap(Dummy_ru_vine_swap& matrixToCopy){}
		Dummy_ru_vine_swap(Dummy_ru_vine_swap&& other) noexcept{}
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
	};

	using Chain_vine_swap_option = typename std::conditional<
											Options::has_vine_update,
											Chain_vine_swap<Matrix<Options>>,
											Dummy_chain_vine_swap
										>::type;


//	using index = unsigned int;
//	using boundary_type = std::vector<index>;
//	using boundary_matrix = std::vector<boundary_type>;

//	Matrix();
//	Matrix(boundary_matrix& boundaries);	//simplex indices have to start at 0 and be consecutifs
//	Matrix(int numberOfColumns);
//	Matrix(Matrix &matrixToCopy);
//	Matrix(Matrix&& other) noexcept;

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

};

} //namespace persistence_matrix
} //namespace Gudhi

#endif // MASTER_MATRIX_H
