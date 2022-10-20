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
#define BOOST_MPL_LIMIT_LIST_SIZE 40
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

							Matrix<Default_options<Zp_field_element<5> > >,
							Matrix<Default_options<Zp_field_element<2> > >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::HEAP> >,

							Matrix<Multi_persistence_options<> >,
							Matrix<Multi_persistence_options<Column_types::LIST> >,
							Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> >,
							Matrix<Multi_persistence_options<Column_types::VECTOR> >,
							Matrix<Multi_persistence_options<Column_types::HEAP> >,

							Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<Cohomology_persistence_options<Zp_field_element<5> > >,
							Matrix<Cohomology_persistence_options<Zp_field_element<2> > >,

							Matrix<test_options1<> >,
							Matrix<test_options1<Column_types::LIST> >,
							Matrix<test_options1<Column_types::UNORDERED_SET> >,
							Matrix<test_options1<Column_types::VECTOR> >,
							Matrix<test_options1<Column_types::HEAP> >,

							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors, Matrix, list_of_matrix_types) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	//default constructor
	Matrix m;
	BOOST_CHECK_EQUAL(m.get_max_dimension(), -1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 0);

	//constructor from given boundary matrix
	Matrix mb(ordered_boundaries);
	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 2);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(6), 2);

	//constructor reserving column space
	Matrix mr(5);
	BOOST_CHECK_EQUAL(mr.get_max_dimension(), -1);
	BOOST_CHECK_EQUAL(mr.get_number_of_columns(), 0);

	//copy constructor
	Matrix mc1(mb);
	Matrix mc2 = mb;
	BOOST_CHECK_EQUAL(mc1.get_max_dimension(), 2);
	BOOST_CHECK_EQUAL(mc1.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mc1.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mc1.get_column_dimension(6), 2);
	BOOST_CHECK_EQUAL(mc2.get_max_dimension(), 2);
	BOOST_CHECK_EQUAL(mc2.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mc2.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mc2.get_column_dimension(6), 2);

	//move constructor
	Matrix mm(std::move(mb));
	BOOST_CHECK_EQUAL(mm.get_max_dimension(), 2);
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mm.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mm.get_column_dimension(6), 2);
	BOOST_CHECK_EQUAL(mb.get_max_dimension(), -1);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 0);

	//swap
	swap(mm, mb);
	BOOST_CHECK_EQUAL(mm.get_max_dimension(), -1);
	BOOST_CHECK_EQUAL(mm.get_number_of_columns(), 0);
	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 2);
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(mb.get_column_dimension(6), 2);
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

							Matrix<Default_options<Zp_field_element<5> > >,
							Matrix<Default_options<Zp_field_element<2> > >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::HEAP> >,

							Matrix<Multi_persistence_options<> >,
							Matrix<Multi_persistence_options<Column_types::LIST> >,
							Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> >,
							Matrix<Multi_persistence_options<Column_types::VECTOR> >,
							Matrix<Multi_persistence_options<Column_types::HEAP> >,

							Matrix<test_options1<> >,
							Matrix<test_options1<Column_types::LIST> >,
							Matrix<test_options1<Column_types::UNORDERED_SET> >,
							Matrix<test_options1<Column_types::VECTOR> >,
							Matrix<test_options1<Column_types::HEAP> >
						> list_of_boundary_matrix_types;

typedef boost::mpl::list<Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<Cohomology_persistence_options<Zp_field_element<5> > >,
							Matrix<Cohomology_persistence_options<Zp_field_element<2> > >,

							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_chain_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_methods, Matrix, list_of_boundary_matrix_types) {
	using boundary_type = typename Matrix::boundary_type;
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	boundary_type boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	boundary_type boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(!m.is_zero_column(3));

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(!m.is_zero_column(3));

	std::set<unsigned int> rowIndices;
	for (auto& cell : m.get_column(6)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 3);
	unsigned int i = 3;
	for (unsigned int r : rowIndices)
		BOOST_CHECK_EQUAL(r, i++);

	rowIndices.clear();
	for (auto& cell : m.get_column(5)){
		rowIndices.insert(cell.get_row_index());
	}

	BOOST_CHECK_EQUAL(m.get_pivot(0), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(1), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(2), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(3), 1);
	BOOST_CHECK_EQUAL(m.get_pivot(4), 2);
	BOOST_CHECK_EQUAL(m.get_pivot(6), 5);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_methods, Matrix, list_of_chain_matrix_types) {
	using boundary_type = typename Matrix::boundary_type;
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	boundary_type boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	boundary_type boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,2));
	BOOST_CHECK(!m.is_zero_cell(1,1));
	BOOST_CHECK(!m.is_zero_cell(3,3));
	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		BOOST_CHECK(!m.is_zero_column(i));
	}

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,2));
	BOOST_CHECK(!m.is_zero_cell(1,1));
	BOOST_CHECK(!m.is_zero_cell(3,3));
	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		BOOST_CHECK(!m.is_zero_column(i));
	}

	std::set<unsigned int> rowIndices;
	for (auto& cell : m.get_column(6)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 1);
	BOOST_CHECK_EQUAL(*rowIndices.begin(), 6);

	BOOST_CHECK_EQUAL(m.get_pivot(0), 0);
	BOOST_CHECK_EQUAL(m.get_pivot(1), 1);
	BOOST_CHECK_EQUAL(m.get_pivot(2), 2);
	BOOST_CHECK_EQUAL(m.get_pivot(3), 3);
	BOOST_CHECK_EQUAL(m.get_pivot(4), 4);
	BOOST_CHECK_EQUAL(m.get_pivot(5), 5);
	BOOST_CHECK_EQUAL(m.get_pivot(6), 6);
}

typedef boost::mpl::list<Matrix<Representative_cycles_options<Zp_field_element<2> > >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Representative_cycles_options<Zp_field_element<2>,Column_types::HEAP> >,

							Matrix<Default_options<Zp_field_element<2> > >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::VECTOR> >,
							Matrix<Default_options<Zp_field_element<2>,Column_types::HEAP> >,

							Matrix<Multi_persistence_options<> >,
							Matrix<Multi_persistence_options<Column_types::LIST> >,
							Matrix<Multi_persistence_options<Column_types::UNORDERED_SET> >,
							Matrix<Multi_persistence_options<Column_types::VECTOR> >,
							Matrix<Multi_persistence_options<Column_types::HEAP> >,

							Matrix<test_options1<> >,
							Matrix<test_options1<Column_types::LIST> >,
							Matrix<test_options1<Column_types::UNORDERED_SET> >,
							Matrix<test_options1<Column_types::VECTOR> >,
							Matrix<test_options1<Column_types::HEAP> >
						> list_of_Z2_boundary_matrix_types;
typedef boost::mpl::list<Matrix<Representative_cycles_options<Zp_field_element<5> > >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Representative_cycles_options<Zp_field_element<5>,Column_types::VECTOR> >,

							Matrix<Default_options<Zp_field_element<5> > >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::LIST> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::UNORDERED_SET> >,
							Matrix<Default_options<Zp_field_element<5>,Column_types::VECTOR> >
						> list_of_Zp_boundary_matrix_types;

typedef boost::mpl::list<Matrix<Zigzag_options<> >,
							Matrix<Zigzag_options<Column_types::LIST> >,

							Matrix<Cohomology_persistence_options<Zp_field_element<2> > >,

							Matrix<test_options2<> >,
							Matrix<test_options2<Column_types::LIST> >
						> list_of_Z2_chain_matrix_types;
typedef boost::mpl::list<Matrix<Cohomology_persistence_options<Zp_field_element<5> > >
						> list_of_Zp_chain_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_boundary_matrix_methods, Matrix, list_of_Z2_boundary_matrix_types) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::set<unsigned int> rowIndices;
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	unsigned int i = 1;
	for (unsigned int r : rowIndices)
		BOOST_CHECK_EQUAL(r, i++);

	m.add_to(3, 4);

	rowIndices.clear();
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	i = 0;
	for (unsigned int r : rowIndices){
		BOOST_CHECK_EQUAL(r, i);
		i = 2;
	}
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Zp_boundary_matrix_methods, Matrix, list_of_Zp_boundary_matrix_types) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	unsigned int i = 1;
	for (const std::pair<unsigned int,unsigned int>&r : rowIndices){
//		std::cout << "1: " << r.first << ", " << r.second << "\n";
		BOOST_CHECK_EQUAL(r.first, i++);
		BOOST_CHECK_EQUAL(r.second, 1);
	}

	m.add_to(3, 4);

	rowIndices.clear();
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 3);
	i = 0;
	for (const std::pair<unsigned int,unsigned int>& r : rowIndices){
//		std::cout << "2: " << r.first << ", " << r.second << "\n";
		if (i != 1) BOOST_CHECK_EQUAL(r.second, 1);
		else BOOST_CHECK_EQUAL(r.second, 2);
		BOOST_CHECK_EQUAL(r.first, i++);
	}
//	std::cout << "end\n";
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_matrix_methods, Matrix, list_of_Z2_chain_matrix_types) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::set<unsigned int> rowIndices;
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	unsigned int i = 3;
	for (unsigned int r : rowIndices){
//		std::cout << "1: " << r << "\n";
		BOOST_CHECK_EQUAL(r, i++);
	}

	m.add_to(3, 4);

	rowIndices.clear();
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 1);
	BOOST_CHECK_EQUAL(*rowIndices.begin(), 4);
//	std::cout << "end\n";
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Zp_chain_matrix_methods, Matrix, list_of_Zp_chain_matrix_types) {
	using boundary_matrix = typename Matrix::boundary_matrix;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	for (auto& cell : m.get_column(4)){
//		std::cout << "1: " << cell.get_row_index() << ", " << cell.get_element() << "\n";
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	unsigned int i = 3;
	for (const std::pair<unsigned int,unsigned int>&r : rowIndices){
		BOOST_CHECK_EQUAL(r.first, i++);
		BOOST_CHECK_EQUAL(r.second, 1);
	}

	m.add_to(3, 4);

	rowIndices.clear();
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
//		std::cout << "2: " << cell.get_row_index() << ", " << cell.get_element() << "\n";
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	i = 3;
	for (const std::pair<unsigned int,unsigned int>&r : rowIndices){
		if (i == 3) BOOST_CHECK_EQUAL(r.second, 2);
		else BOOST_CHECK_EQUAL(r.second, 1);
		BOOST_CHECK_EQUAL(r.first, i++);
	}
}

