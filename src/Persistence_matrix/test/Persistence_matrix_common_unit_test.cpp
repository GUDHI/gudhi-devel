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
// #include "gudhi/utilities/utilities.h"

using Gudhi::persistence_matrix::Z2_field_element;
using Gudhi::persistence_matrix::Zp_field_element;
using Gudhi::persistence_matrix::Matrix;
// using Gudhi::persistence_matrix::Representative_cycles_options;
using Gudhi::persistence_matrix::Default_options;
// using Gudhi::persistence_matrix::Zigzag_options;
// using Gudhi::persistence_matrix::Multi_persistence_options;
// using Gudhi::persistence_matrix::Cohomology_persistence_options;
using Gudhi::persistence_matrix::Column_types;

using Z5 = Zp_field_element<5>;
using Z2 = Z2_field_element;

template<class Field_type, Column_types column_type>
struct opt_b_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_p_nb : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_nb : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_r_p_nb : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_r_nb : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_i_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_i_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_i : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_i_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_b_ra_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_i_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_i : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_i_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_r_p : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_i_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
};

void build_boundary_matrix(std::vector<std::vector<unsigned int> >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back({0,1});
	boundaries.push_back({1,2});
	boundaries.push_back({0,2});
	boundaries.push_back({3,4,5});
}

template<typename Field_type>
void build_boundary_matrix(std::vector<std::vector<std::pair<unsigned int,Field_type> > >& boundaries)
{
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.emplace_back();
	boundaries.push_back({{0,Field_type(1)},{1,Field_type(4)}});
	boundaries.push_back({{1,Field_type(1)},{2,Field_type(4)}});
	boundaries.push_back({{0,Field_type(1)},{2,Field_type(4)}});
	boundaries.push_back({{3,Field_type(1)},{4,Field_type(1)},{5,Field_type(4)}});
}

typedef boost::mpl::list</*Matrix<opt_ra_i_p<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i_p<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i_p<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::LIST> >,
							Matrix<opt_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_p<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_p<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z2,Column_types::LIST> >,
							Matrix<opt_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_p<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_p<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_p<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z5,Column_types::LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::LIST> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::LIST> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra<Z2,Column_types::INTRUSIVE_SET> >*/
						> list_of_matrix_types_chain_with_row_access_no_removals;

typedef boost::mpl::list<Matrix<opt_b_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_matrix_types;

typedef boost::mpl::list<Matrix<opt_b_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra<Z5,Column_types::SET> >,
							Matrix<opt_b_ra<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >
						> list_of_matrix_types2;

template<class Matrix>
void test_constructors(){
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

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

BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors, Matrix, list_of_matrix_types) {
	test_constructors<Matrix>();
}
BOOST_AUTO_TEST_CASE_TEMPLATE(Matrix_constructors_bis, Matrix, list_of_matrix_types2) {
	test_constructors<Matrix>();
}

typedef boost::mpl::list<Matrix<opt_b_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra<Z5,Column_types::SET> >,
							Matrix<opt_b_ra<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_base_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_matrix_methods, Matrix, list_of_base_matrix_types) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	auto boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	auto boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 1);

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));
	BOOST_CHECK(!m.is_zero_column(5));
	BOOST_CHECK(!m.is_zero_column(6));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 2);

	BOOST_CHECK_EQUAL(m.get_column_dimension(0), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(1), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(3), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(4), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(5), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(6), 2);

	BOOST_CHECK(!m.is_zero_cell(3,1));
	BOOST_CHECK(!m.is_zero_column(3));
	m.zero_cell(3, 1);
	BOOST_CHECK(m.is_zero_cell(3,1));
	BOOST_CHECK(!m.is_zero_column(3));
	m.zero_column(3);
	BOOST_CHECK(m.is_zero_cell(3,1));
	BOOST_CHECK(m.is_zero_column(3));
}

typedef boost::mpl::list<Matrix<opt_b_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_boundary_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_matrix_methods, Matrix, list_of_boundary_matrix_types) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	auto boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	auto boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 1);

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,0));
	BOOST_CHECK(!m.is_zero_cell(3,0));
	BOOST_CHECK(m.is_zero_column(0));
	BOOST_CHECK(m.is_zero_column(1));
	BOOST_CHECK(m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));
	BOOST_CHECK(m.is_zero_column(5));
	BOOST_CHECK(!m.is_zero_column(6));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 2);

	BOOST_CHECK_EQUAL(m.get_column_dimension(0), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(1), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(3), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(4), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(5), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(6), 2);

	BOOST_CHECK_EQUAL(m.get_pivot(0), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(1), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(2), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(3), 1);
	BOOST_CHECK_EQUAL(m.get_pivot(4), 2);
	BOOST_CHECK_EQUAL(m.get_pivot(5), -1);
	BOOST_CHECK_EQUAL(m.get_pivot(6), 5);

	BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 3);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 4);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 6);
}

typedef boost::mpl::list<Matrix<opt_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >
						> list_of_chain_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_matrix_methods, Matrix, list_of_chain_matrix_types) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	auto boundary2 = ordered_boundaries.back();
	ordered_boundaries.pop_back();
	auto boundary1 = ordered_boundaries.back();
	ordered_boundaries.pop_back();

	Matrix m(ordered_boundaries);
	BOOST_CHECK(m.is_zero_cell(1,2));
	BOOST_CHECK(!m.is_zero_cell(1,1));
	BOOST_CHECK(!m.is_zero_cell(3,3));
	BOOST_CHECK(!m.is_zero_column(0));
	BOOST_CHECK(!m.is_zero_column(1));
	BOOST_CHECK(!m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 1);

	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 5);
	m.insert_boundary(boundary1);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 6);
	m.insert_boundary(boundary2);
	BOOST_CHECK_EQUAL(m.get_number_of_columns(), 7);
	BOOST_CHECK(m.is_zero_cell(1,2));
	BOOST_CHECK(!m.is_zero_cell(1,1));
	BOOST_CHECK(!m.is_zero_cell(3,3));
	BOOST_CHECK(!m.is_zero_column(0));
	BOOST_CHECK(!m.is_zero_column(1));
	BOOST_CHECK(!m.is_zero_column(2));
	BOOST_CHECK(!m.is_zero_column(3));
	BOOST_CHECK(!m.is_zero_column(4));
	BOOST_CHECK(!m.is_zero_column(5));
	BOOST_CHECK(!m.is_zero_column(6));

	BOOST_CHECK_EQUAL(m.get_max_dimension(), 2);

	BOOST_CHECK_EQUAL(m.get_column_dimension(0), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(1), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(2), 0);
	BOOST_CHECK_EQUAL(m.get_column_dimension(3), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(4), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(5), 1);
	BOOST_CHECK_EQUAL(m.get_column_dimension(6), 2);

	BOOST_CHECK_EQUAL(m.get_pivot(0), 0);
	BOOST_CHECK_EQUAL(m.get_pivot(1), 1);
	BOOST_CHECK_EQUAL(m.get_pivot(2), 2);
	BOOST_CHECK_EQUAL(m.get_pivot(3), 3);
	BOOST_CHECK_EQUAL(m.get_pivot(4), 4);
	BOOST_CHECK_EQUAL(m.get_pivot(5), 5);
	BOOST_CHECK_EQUAL(m.get_pivot(6), 6);

	BOOST_CHECK_EQUAL(m.get_column_with_pivot(0), 0);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(1), 1);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(2), 2);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(3), 3);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(4), 4);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(5), 5);
	BOOST_CHECK_EQUAL(m.get_column_with_pivot(6), 6);
}

typedef boost::mpl::list<Matrix<opt_b_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_Z2_base_matrix_types;
typedef boost::mpl::list<Matrix<opt_b_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_Z5_base_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_base_matrix_methods, Matrix, list_of_Z2_base_matrix_types) {
	std::vector<std::vector<unsigned int> > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	unsigned int i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		std::set<unsigned int> rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert(cell.get_row_index());
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(*itCol, *itOrdB);
			itCol++; itOrdB++;
		}
	}

	m.add_to(3, 4);

	std::set<unsigned int> rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	BOOST_CHECK_EQUAL(*rowIndices.begin(), 0);
	BOOST_CHECK_EQUAL(*rowIndices.rbegin(), 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z5_base_matrix_methods, Matrix, list_of_Z5_base_matrix_types) {
	std::vector<std::vector<std::pair<unsigned int,Z5> > > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	unsigned int i = 0;
	for (auto& b : ordered_boundaries){
		const auto& col = m.get_column(i++);
		std::set<std::pair<unsigned int,unsigned int> > rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert({cell.get_row_index(), cell.get_element()});
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(itCol->first, itOrdB->first);
			BOOST_CHECK_EQUAL(itCol->second, itOrdB->second);
			itCol++; itOrdB++;
		}
	}

	m.add_to(3, 4);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	BOOST_CHECK_EQUAL(rowIndices.begin()->first, 0);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->first, 2);
	BOOST_CHECK_EQUAL(rowIndices.begin()->second, 1);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->second, 4);
}

typedef boost::mpl::list<Matrix<opt_b_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_Z2_boundary_matrix_types;
typedef boost::mpl::list<Matrix<opt_b_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_Z5_boundary_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_boundary_matrix_methods, Matrix, list_of_Z2_boundary_matrix_types) {
	std::vector<std::vector<unsigned int> > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		if (i == 5) continue;

		const auto& b = ordered_boundaries[i];
		const auto& col = m.get_column(i);
		std::set<unsigned int> rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert(cell.get_row_index());
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(*itCol, *itOrdB);
			itCol++; itOrdB++;
		}
	}
	BOOST_CHECK(m.get_column(5).is_empty());

	m.add_to(3, 4);

	std::set<unsigned int> rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	BOOST_CHECK_EQUAL(*rowIndices.begin(), 0);
	BOOST_CHECK_EQUAL(*rowIndices.rbegin(), 2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z5_boundary_matrix_methods, Matrix, list_of_Z5_boundary_matrix_types) {
	std::vector<std::vector<std::pair<unsigned int,Z5> > > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		if (i == 5) continue;

		const auto& b = ordered_boundaries[i];
		const auto& col = m.get_column(i);
		std::set<std::pair<unsigned int,unsigned int> > rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert({cell.get_row_index(), cell.get_element()});
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(itCol->first, itOrdB->first);
			BOOST_CHECK_EQUAL(itCol->second, itOrdB->second);
			itCol++; itOrdB++;
		}
	}
	BOOST_CHECK(m.get_column(5).is_empty());

	m.add_to(3, 4);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	BOOST_CHECK_EQUAL(rowIndices.begin()->first, 0);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->first, 2);
	BOOST_CHECK_EQUAL(rowIndices.begin()->second, 1);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->second, 4);
}

typedef boost::mpl::list<Matrix<opt_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >
						> list_of_Z2_chain_matrix_types;
typedef boost::mpl::list<Matrix<opt_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >
						> list_of_Z5_chain_matrix_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_matrix_methods, Matrix, list_of_Z2_chain_matrix_types) {
	std::vector<std::vector<unsigned int> > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::vector<std::vector<unsigned int> > columns;
	columns.push_back({0});
	columns.push_back({0,1});
	columns.push_back({0,2});
	columns.push_back({3});
	columns.push_back({3,4});
	columns.push_back({3,4,5});
	columns.push_back({6});

	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		const auto& b = columns[i];
		const auto& col = m.get_column(i);
		std::set<unsigned int> rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert(cell.get_row_index());
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(*itCol, *itOrdB);
			itCol++; itOrdB++;
		}
	}

	m.add_to(3, 4);

	std::set<unsigned int> rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert(cell.get_row_index());
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 1);
	BOOST_CHECK_EQUAL(*rowIndices.begin(), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z5_chain_matrix_methods, Matrix, list_of_Z5_chain_matrix_types) {
	std::vector<std::vector<std::pair<unsigned int,Z5> > > ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);
	Matrix m(ordered_boundaries);

	std::vector<std::vector<std::pair<unsigned int,Z5> > > columns;
	columns.push_back({{0,Z5(1)}});
	columns.push_back({{0,Z5(1)},{1,Z5(4)}});
	columns.push_back({{0,Z5(1)},{2,Z5(4)}});
	columns.push_back({{3,Z5(1)}});
	columns.push_back({{3,Z5(1)},{4,Z5(1)}});
	columns.push_back({{3,Z5(1)},{4,Z5(1)},{5,Z5(4)}});
	columns.push_back({{6,Z5(1)}});

	for (unsigned int i = 0; i < ordered_boundaries.size(); ++i){
		const auto& b = columns[i];
		const auto& col = m.get_column(i);
		std::set<std::pair<unsigned int,unsigned int> > rowIndices;
		//not all columns are ordered
		for (auto& cell : col){
			rowIndices.insert({cell.get_row_index(), cell.get_element()});
		}
		auto itCol = rowIndices.begin();
		auto itOrdB = b.begin();
		while (itCol != rowIndices.end()) {
			BOOST_CHECK_EQUAL(itCol->first, itOrdB->first);
			BOOST_CHECK_EQUAL(itCol->second, itOrdB->second);
			itCol++; itOrdB++;
		}
	}

	m.add_to(3, 4);

	std::set<std::pair<unsigned int,unsigned int> > rowIndices;
	//not all columns are ordered
	for (auto& cell : m.get_column(4)){
		rowIndices.insert({cell.get_row_index(), cell.get_element()});
	}
	BOOST_CHECK_EQUAL(rowIndices.size(), 2);
	BOOST_CHECK_EQUAL(rowIndices.begin()->first, 3);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->first, 4);
	BOOST_CHECK_EQUAL(rowIndices.begin()->second, 2);
	BOOST_CHECK_EQUAL(rowIndices.rbegin()->second, 1);
}

typedef boost::mpl::list<Matrix<opt_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z5,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >
						> list_of_z5_chain_options_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_row_access_columns_options, Matrix, list_of_z5_chain_options_with_row_access) {
	using boundary_matrix = std::vector<std::vector<std::pair<unsigned int,Z5> > >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	rows.push_back({{0,Z5(1)},{1,Z5(1)},{2,Z5(1)}});
	rows.push_back({{1,Z5(4)}});
	rows.push_back({{2,Z5(4)}});
	rows.push_back({{3,Z5(1)},{4,Z5(1)},{5,Z5(1)}});
	rows.push_back({{4,Z5(1)},{5,Z5(1)}});
	rows.push_back({{5,Z5(4)}});
	rows.push_back({{6,Z5(1)}});

	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (auto& cell : mb.get_row(i)){
			ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->first, itChain->first);
			BOOST_CHECK_EQUAL(itRow->second, itChain->second);
			++itRow; ++itChain;
		}
	}
}

typedef boost::mpl::list<Matrix<opt_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >
						> list_of_z2_chain_options_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_row_access_columns_options, Matrix, list_of_z2_chain_options_with_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<unsigned int> > rows;
	rows.push_back({0,1,2});
	rows.push_back({1});
	rows.push_back({2});
	rows.push_back({3,4,5});
	rows.push_back({4,5});
	rows.push_back({5});
	rows.push_back({6});

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (auto& cell : mb.get_row(i)){
			ordered_rows[i].insert(cell.get_column_index());
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(*itRow, *itChain);
			++itRow; ++itChain;
		}
	}
}

typedef boost::mpl::list<Matrix<opt_b_ra_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >
						> list_of_z5_base_options_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_row_access_columns_options, Matrix, list_of_z5_base_options_with_row_access) {
	using boundary_matrix = std::vector<std::vector<std::pair<unsigned int,Z5> > >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	rows.push_back({{3,Z5(1)},{5,Z5(1)}});
	rows.push_back({{3,Z5(4)},{4,Z5(1)}});
	rows.push_back({{4,Z5(4)},{5,Z5(4)}});
	rows.push_back({{6,Z5(1)}});
	rows.push_back({{6,Z5(1)}});
	rows.push_back({{6,Z5(4)}});

	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (auto& cell : mb.get_row(i)){
			ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->first, itChain->first);
			BOOST_CHECK_EQUAL(itRow->second, itChain->second);
			++itRow; ++itChain;
		}
	}
}

typedef boost::mpl::list<Matrix<opt_b_ra_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_i_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_z2_base_options_with_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_base_row_access_columns_options, Matrix, list_of_z2_base_options_with_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	std::vector<std::vector<unsigned int> > rows;
	rows.push_back({3,5});
	rows.push_back({3,4});
	rows.push_back({4,5});
	rows.push_back({6});
	rows.push_back({6});
	rows.push_back({6});

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		for (auto& cell : mb.get_row(i)){
			ordered_rows[i].insert(cell.get_column_index());
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		BOOST_CHECK_EQUAL(row.size(), rows[i].size());
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(*itRow, *itChain);
			++itRow; ++itChain;
		}
	}
}

typedef boost::mpl::list<Matrix<opt_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_p_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_r_nb<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_r_p<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z5,Column_types::SET> >,
							Matrix<opt_b_ra_i_r_p<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_r_p<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_r<Z2,Column_types::SET> >,
							Matrix<opt_b_ra_i_r_p<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_b_ra_i_r<Z2,Column_types::INTRUSIVE_SET> >
						> list_of_options_with_remove;

BOOST_AUTO_TEST_CASE_TEMPLATE(Removable_columns_options, Matrix, list_of_options_with_remove) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 2);
	if constexpr (Matrix::Option_list::has_vine_update){
		const auto& barcode = mb.get_current_barcode();
		BOOST_CHECK_EQUAL(barcode.back().death, 6);
	}
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 7);
	mb.remove_maximal_simplex(6);
	BOOST_CHECK_EQUAL(mb.get_max_dimension(), 1);
	if constexpr (Matrix::Option_list::has_vine_update){
		const auto& barcode2 = mb.get_current_barcode();
		BOOST_CHECK_EQUAL(barcode2.back().death, -1);
	}
	BOOST_CHECK_EQUAL(mb.get_number_of_columns(), 6);
}

//TODO: test removal with row access

