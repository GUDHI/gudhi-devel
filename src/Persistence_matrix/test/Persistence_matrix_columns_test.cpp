/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_matrix"
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 30
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <utility>

#include "gudhi/utilities/Zp_field.h"
#include "gudhi/utilities/utilities.h"
#include "gudhi/matrix.h"
#include "gudhi/options.h"
#include "gudhi/column_types/list_column.h"
#include "gudhi/column_types/set_column.h"
#include "gudhi/column_types/unordered_set_column.h"
#include "gudhi/column_types/vector_column.h"
#include "gudhi/column_types/reduced_cell_list_column_with_row.h"
#include "gudhi/column_types/reduced_cell_set_column_with_row.h"
#include "gudhi/column_types/z2_heap_column.h"
#include "gudhi/column_types/z2_list_column.h"
#include "gudhi/column_types/z2_set_column.h"
#include "gudhi/column_types/z2_unordered_set_column.h"
#include "gudhi/column_types/z2_vector_column.h"
#include "gudhi/column_types/z2_reduced_cell_list_column_with_row.h"
#include "gudhi/column_types/z2_reduced_cell_set_column_with_row.h"
#include "gudhi/column_types/column_pairing.h"

using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::List_column;
using Gudhi::persistence_matrix::Set_column;
using Gudhi::persistence_matrix::Unordered_set_column;
using Gudhi::persistence_matrix::Vector_column;
using Gudhi::persistence_matrix::Reduced_cell_list_column_with_row;
using Gudhi::persistence_matrix::Reduced_cell_set_column_with_row;
using Gudhi::persistence_matrix::Z2_heap_column;
using Gudhi::persistence_matrix::Z2_list_column;
using Gudhi::persistence_matrix::Z2_set_column;
using Gudhi::persistence_matrix::Z2_unordered_set_column;
using Gudhi::persistence_matrix::Z2_vector_column;
using Gudhi::persistence_matrix::Z2_reduced_cell_list_column_with_row;
using Gudhi::persistence_matrix::Z2_reduced_cell_set_column_with_row;
using Gudhi::persistence_matrix::Column_pairing;
using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Default_options;

struct set_test_options_with_pairing : Default_options<Zp_field_element<5>, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct set_test_options_without_pairing : Default_options<Zp_field_element<5>, Column_types::SET, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

struct list_test_options_with_pairing : Default_options<Zp_field_element<5>, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_row_access = true;
};

struct list_test_options_without_pairing : Default_options<Zp_field_element<5>, Column_types::LIST, false, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_row_access = true;
};

using Set_matrix_with_pairing = Matrix<set_test_options_with_pairing>;
using Set_matrix_without_pairing = Matrix<set_test_options_without_pairing>;
using List_matrix_with_pairing = Matrix<list_test_options_with_pairing>;
using List_matrix_without_pairing = Matrix<list_test_options_without_pairing>;
using Dummy_pairing = Set_matrix_without_pairing::Dummy_column_pairing;

void build_boundary_matrix(std::vector<std::vector<unsigned int> >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back(std::vector<unsigned int>{0,1});
	boundaries.push_back(std::vector<unsigned int>{1,2});
	boundaries.push_back(std::vector<unsigned int>{0,2});
	boundaries.push_back(std::vector<unsigned int>{3,4,5});
	boundaries.emplace_back();
	boundaries.push_back(std::vector<unsigned int>{1,7});
}

template<typename Field_type>
void build_boundary_matrix(std::vector<std::vector<std::pair<unsigned int,Field_type> > >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(1)},{1,Field_type(1)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{1,Field_type(1)},{2,Field_type(1)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{0,Field_type(4)},{2,Field_type(1)}});
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{3,Field_type(1)},{4,Field_type(1)},{5,Field_type(1)}});
	boundaries.emplace_back();
	boundaries.push_back(std::vector<std::pair<unsigned int,Field_type> >{{1,Field_type(1)},{7,Field_type(1)}});
}

typedef boost::mpl::list<List_column<Zp_field_element<5>, Column_pairing>,
							List_column<Zp_field_element<5>, Dummy_pairing>,
							Set_column<Zp_field_element<5>, Column_pairing>,
							Set_column<Zp_field_element<5>, Dummy_pairing>,
							Unordered_set_column<Zp_field_element<5>, Column_pairing>,
							Unordered_set_column<Zp_field_element<5>, Dummy_pairing>,
							Vector_column<Zp_field_element<5>, Column_pairing>,
							Vector_column<Zp_field_element<5>, Dummy_pairing>,
							Reduced_cell_list_column_with_row<List_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Reduced_cell_list_column_with_row<List_matrix_without_pairing, Zp_field_element<5>, Dummy_pairing>,
							Reduced_cell_set_column_with_row<Set_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Reduced_cell_set_column_with_row<Set_matrix_without_pairing, Zp_field_element<5>, Dummy_pairing>,
							Z2_heap_column<Column_pairing>,
							Z2_heap_column<Dummy_pairing>,
							Z2_list_column<Column_pairing>,
							Z2_list_column<Dummy_pairing>,
							Z2_set_column<Column_pairing>,
							Z2_set_column<Dummy_pairing>,
							Z2_unordered_set_column<Column_pairing>,
							Z2_unordered_set_column<Dummy_pairing>,
							Z2_vector_column<Column_pairing>,
							Z2_vector_column<Dummy_pairing>,
							Z2_reduced_cell_list_column_with_row<List_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_list_column_with_row<List_matrix_without_pairing, Dummy_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_matrix_without_pairing, Dummy_pairing>
						> list_of_columns;

//List_column& operator+=(List_column const &column);
//friend List_column operator+(List_column column1, List_column const& column2);
//List_column& operator*=(unsigned int v);
//friend List_column operator*(List_column column, unsigned int const& v);
//friend List_column operator*(unsigned int const& v, List_column const column);

//iterator begin() noexcept;
//const_iterator begin() const noexcept;
//iterator end() noexcept;
//const_iterator end() const noexcept;

//int get_pivot() const;
//Field_element_type get_pivot_value() const;
//dimension_type get_dimension() const;
//bool is_empty() const;
//bool is_non_zero(index rowIndex) const;

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common, Column, list_of_columns) {
}

typedef boost::mpl::list<List_column<Zp_field_element<5>, Column_pairing>,
							Set_column<Zp_field_element<5>, Column_pairing>,
							Unordered_set_column<Zp_field_element<5>, Column_pairing>,
							Vector_column<Zp_field_element<5>, Column_pairing>,
							Reduced_cell_list_column_with_row<List_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Reduced_cell_set_column_with_row<Set_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Z2_heap_column<Column_pairing>,
							Z2_list_column<Column_pairing>,
							Z2_set_column<Column_pairing>,
							Z2_unordered_set_column<Column_pairing>,
							Z2_vector_column<Column_pairing>,
							Z2_reduced_cell_list_column_with_row<List_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_matrix_with_pairing, Column_pairing>
						> list_of_columns_with_pairing;

//index get_paired_chain_index() const;
//bool is_paired() const;
//void assign_paired_chain(index other_col);
//void unassign_paired_chain();

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_pairing_option, Column, list_of_columns_with_pairing) {
}

typedef boost::mpl::list<List_column<Zp_field_element<5>, Column_pairing>,
							List_column<Zp_field_element<5>, Dummy_pairing>,
							Set_column<Zp_field_element<5>, Column_pairing>,
							Set_column<Zp_field_element<5>, Dummy_pairing>,
							Unordered_set_column<Zp_field_element<5>, Column_pairing>,
							Unordered_set_column<Zp_field_element<5>, Dummy_pairing>,
							Vector_column<Zp_field_element<5>, Column_pairing>,
							Vector_column<Zp_field_element<5>, Dummy_pairing>,
							Z2_heap_column<Column_pairing>,
							Z2_heap_column<Dummy_pairing>,
							Z2_list_column<Column_pairing>,
							Z2_list_column<Dummy_pairing>,
							Z2_set_column<Column_pairing>,
							Z2_set_column<Dummy_pairing>,
							Z2_unordered_set_column<Column_pairing>,
							Z2_unordered_set_column<Dummy_pairing>,
							Z2_vector_column<Column_pairing>,
							Z2_vector_column<Dummy_pairing>
						> list_of_columns_without_row;

//List_column();
//template<class Boundary_type>
//List_column(const Boundary_type& boundary);
//template<class Boundary_type>
//List_column(const Boundary_type& boundary, dimension_type dimension);
//List_column(const List_column& column);
//List_column(List_column&& column) noexcept;

//std::vector<Field_element_type> get_content(unsigned int columnLength) const;
//void clear();
//void clear(index rowIndex);
//void reorder(std::vector<index>& valueMap);

//List_column& operator=(List_column other);
//friend void swap(List_column& col1, List_column& col2);

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_without_row_methods, Column, list_of_columns_without_row) {
}

typedef boost::mpl::list<Reduced_cell_list_column_with_row<List_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Reduced_cell_list_column_with_row<List_matrix_without_pairing, Zp_field_element<5>, Dummy_pairing>,
							Reduced_cell_set_column_with_row<Set_matrix_with_pairing, Zp_field_element<5>, Column_pairing>,
							Reduced_cell_set_column_with_row<Set_matrix_without_pairing, Zp_field_element<5>, Dummy_pairing>,
							Z2_reduced_cell_list_column_with_row<List_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_list_column_with_row<List_matrix_without_pairing, Dummy_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_matrix_with_pairing, Column_pairing>,
							Z2_reduced_cell_set_column_with_row<Set_matrix_without_pairing, Dummy_pairing>
						> list_of_columns_with_row;

//Reduced_cell_list_column_with_row(matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
//template<class Chain_type>
//Reduced_cell_list_column_with_row(index chainIndex, const Chain_type& chain, dimension_type dimension, matrix_type& matrix, dictionnary_type& pivotToColumnIndex);
//Reduced_cell_list_column_with_row(const Reduced_cell_list_column_with_row& other);

//Column_type& get_column();
//const Column_type& get_column() const;
//Row_type& get_row();
//const Row_type& get_row() const;
//int get_lowest_simplex_index() const;

//void swap_rows(Reduced_cell_column_with_row& other);
//void swap_lowest_simplex_index(Reduced_cell_column_with_row& other);

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_methods, Column, list_of_columns_with_row) {
}


