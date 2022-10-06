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

#include "gudhi/matrix.h"
#include "gudhi/options.h"
#include "gudhi/utilities/Z2_field.h"
#include "gudhi/utilities/Zp_field.h"
#include "gudhi/utilities/utilities.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::Matrix;
using Gudhi::persistence_matrix::Representative_cycles_options;
using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Zigzag_options;
using Gudhi::persistence_matrix::Multi_persistence_options;
using Gudhi::persistence_matrix::Cohomology_persistence_options;
using Gudhi::persistence_matrix::Column_types;

template<Column_types column_type = Column_types::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct test_options1 : Default_options<Z2_field_element, column_type, separated_by_dimension, parallelizable>{
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<Column_types column_type = Column_types::SET, bool separated_by_dimension = false, bool parallelizable = false>
struct test_options2 : Default_options<Z2_field_element, column_type, separated_by_dimension, parallelizable>{
	static const bool has_row_access = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_of_boundary_type = false;
	static const bool has_removable_columns = true;
	static const bool is_indexed_by_position = true;
};

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

typedef boost::mpl::list<Matrix<Default_options<Zp_field_element<5> > >,
							Matrix<Default_options<Zp_field_element<2> > >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::HEAP> >
						> list_of_options_with_base_matrix;


BOOST_AUTO_TEST_CASE_TEMPLATE(Zero_options, Matrix, list_of_options_with_base_matrix) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	BOOST_CHECK(!mb.is_zero_cell(3,0));
	mb.zero_cell(3,0);
	BOOST_CHECK(mb.is_zero_cell(3,0));

	BOOST_CHECK(!mb.is_zero_column(6));
	mb.zero_column(6);
	BOOST_CHECK(mb.is_zero_column(6));
}

typedef boost::mpl::list<Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<Cohomology_persistence_options<Zp_field_element<5> > >,
							Matrix<Cohomology_persistence_options<Zp_field_element<2> > >,

							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_options_with_row_access_and_remove;


BOOST_AUTO_TEST_CASE_TEMPLATE(Row_access_and_removable_columns_options, Matrix, list_of_options_with_row_access_and_remove) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	ordered_boundaries.pop_back();
	ordered_boundaries.pop_back();

	Matrix mb(ordered_boundaries);

	BOOST_CHECK_EQUAL(mb.get_row(0).size(), 3);
	unsigned int i = 0;
	for (auto& c : mb.get_row(0)){
		BOOST_CHECK_EQUAL(c.get_row_index(), 0);
		BOOST_CHECK_EQUAL(c.get_column_index(), i++);
	}

	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 2);
	if constexpr (Matrix::Option_list::has_vine_update){
		const auto& barcode = mb.get_current_barcode();
		BOOST_CHECK_EQUAL(barcode.back().death, 6);
	}
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	mb.erase_last();
	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 1);
	if constexpr (Matrix::Option_list::has_vine_update){
		const auto& barcode2 = mb.get_current_barcode();
		BOOST_CHECK_EQUAL(barcode2.back().death, -1);
	}
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 6);
}

typedef boost::mpl::list<Matrix<Representative_cycles_options<Zp_field_element<5> > >,
							Matrix<Representative_cycles_options<Zp_field_element<2> > >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::VECTOR> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::HEAP> >,

							Matrix<Multi_persistence_options<> >,
							Matrix<Multi_persistence_options<Column_types::LIST> >,
							Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> >,
							Matrix<Multi_persistence_options<Column_types::VECTOR> >,
							Matrix<Multi_persistence_options<Column_types::HEAP> >,

							Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<test_options1<> >,
							Matrix<test_options1<Column_types::LIST> >,
							Matrix<test_options1<Column_types::UNORDERED_SET> >,
							Matrix<test_options1<Column_types::VECTOR> >,
							Matrix<test_options1<Column_types::HEAP> >,

							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_options_with_barcode_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Barcode_options, Matrix, list_of_options_with_barcode_access) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	const auto& barcode = mb.get_current_barcode();
	auto it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 3);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 5);
	BOOST_CHECK_EQUAL(it->death, 6);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 7);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK(it == barcode.end());

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 8);
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
}

typedef boost::mpl::list<Matrix<Representative_cycles_options<Zp_field_element<5> > >,
							Matrix<Representative_cycles_options<Zp_field_element<2> > >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::VECTOR> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::HEAP> >
						> list_of_options_with_rep_cycles;

BOOST_AUTO_TEST_CASE_TEMPLATE(Representative_cycle_options, Matrix, list_of_options_with_rep_cycles) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.update_representative_cycles();
	const auto& cycles = mb.get_representative_cycles();
	const auto& barcode = mb.get_current_barcode();
	auto it = barcode.begin();
	BOOST_CHECK_EQUAL(cycles.size(), 5);
	for (auto& cycle : cycles){
		BOOST_CHECK(cycle == mb.get_representative_cycle(*it));
		++it;
	}

	std::vector<unsigned int> tmp;
	tmp.push_back(0);
	BOOST_CHECK(cycles[0] == tmp);
	tmp[0] = 1;
	BOOST_CHECK(cycles[1] == tmp);
	tmp[0] = 2;
	BOOST_CHECK(cycles[2] == tmp);
	tmp[0] = 3;
	tmp.push_back(4);
	tmp.push_back(5);
	BOOST_CHECK(cycles[3] == tmp);
	tmp.clear();
	tmp.push_back(7);
	BOOST_CHECK(cycles[4] == tmp);
}

typedef boost::mpl::list<Matrix<Multi_persistence_options<> >,
							Matrix<Multi_persistence_options<Column_types::LIST> >,
							Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> >,
							Matrix<Multi_persistence_options<Column_types::VECTOR> >,
							Matrix<Multi_persistence_options<Column_types::HEAP> >,
							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_matrix_options_with_vine_and_position_index;

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_vine_option_with_position_indexing, Matrix, list_of_matrix_options_with_vine_and_position_index) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);
	bool change;

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(5));
		BOOST_CHECK(mb.is_zero_column(7));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	const auto& barcode = mb.get_current_barcode();
	auto it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 3);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 5);
	BOOST_CHECK_EQUAL(it->death, 6);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 7);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK(it == barcode.end());

	change = mb.vine_swap(6);
	BOOST_CHECK(change);
	change = mb.vine_swap(5);
	BOOST_CHECK(change);
	change = mb.vine_swap(4);
	BOOST_CHECK(change);
	change = mb.vine_swap(3);
	BOOST_CHECK(change);
	change = mb.vine_swap(7);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(3));
		BOOST_CHECK(mb.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 7);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	change = mb.vine_swap(0);
	BOOST_CHECK(!change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(3));
		BOOST_CHECK(mb.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 7);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	change = mb.vine_swap_with_z_eq_1_case(4);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(3));
		BOOST_CHECK(mb.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 7);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	change = mb.vine_swap(5);
	BOOST_CHECK(!change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(3));
		BOOST_CHECK(mb.is_zero_column(6));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 7);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	change = mb.vine_swap_with_z_eq_1_case(6);
	BOOST_CHECK(change);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(3));
		BOOST_CHECK(mb.is_zero_column(7));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 6);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 7);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 6);
	++it;
	BOOST_CHECK(it == barcode.end());
}

typedef boost::mpl::list<Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<test_options1<> >,
							Matrix<test_options1<Column_types::LIST> >,
							Matrix<test_options1<Column_types::UNORDERED_SET> >,
							Matrix<test_options1<Column_types::VECTOR> >,
							Matrix<test_options1<Column_types::HEAP> >
						> list_of_options_with_vine_and_id_index;

BOOST_AUTO_TEST_CASE_TEMPLATE(Vine_option_with_id_indexing, Matrix, list_of_options_with_vine_and_id_index) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);
	unsigned int next;

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 8);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(5));
		BOOST_CHECK(mb.is_zero_column(7));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	const auto& barcode = mb.get_current_barcode();
	auto it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 3);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 5);
	BOOST_CHECK_EQUAL(it->death, 6);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 7);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK(it == barcode.end());

	next = mb.vine_swap(6, 7);
	BOOST_CHECK_EQUAL(next, 6);
	next = mb.vine_swap(5, 7);
	BOOST_CHECK_EQUAL(next, 5);
	next = mb.vine_swap(4, 7);
	BOOST_CHECK_EQUAL(next, 4);
	next = mb.vine_swap(3, 7);
	BOOST_CHECK_EQUAL(next, 3);
	next = mb.vine_swap(6, 8);
	BOOST_CHECK_EQUAL(next, 6);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	next = mb.vine_swap(0, 1);
	BOOST_CHECK_EQUAL(next, 1);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	next = mb.vine_swap_with_z_eq_1_case(3, 4);
	BOOST_CHECK_EQUAL(next, 3);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	next = mb.vine_swap(3, 5);
	BOOST_CHECK_EQUAL(next, 5);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 6);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 7);
	++it;
	BOOST_CHECK(it == barcode.end());

	next = mb.vine_swap_with_z_eq_1_case(5, 8);
	BOOST_CHECK_EQUAL(next, 5);

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(5));
	} else {
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(0), 1);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 0);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(5), 3);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(8), 8);
	}
	it = barcode.begin();
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 0);
	BOOST_CHECK_EQUAL(it->death, -1);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 1);
	BOOST_CHECK_EQUAL(it->death, 5);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 2);
	BOOST_CHECK_EQUAL(it->death, 4);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 1);
	BOOST_CHECK_EQUAL(it->birth, 7);
	BOOST_CHECK_EQUAL(it->death, 8);
	++it;
	BOOST_CHECK_EQUAL(it->dim, 0);
	BOOST_CHECK_EQUAL(it->birth, 3);
	BOOST_CHECK_EQUAL(it->death, 6);
	++it;
	BOOST_CHECK(it == barcode.end());
}

