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
#define BOOST_MPL_LIMIT_LIST_SIZE 50
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
struct opt_bar_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_bar_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_bar : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_bar_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
};

template<class Field_type, Column_types column_type>
struct opt_rep_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
};

template<class Field_type, Column_types column_type>
struct opt_rep_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
};

template<class Field_type, Column_types column_type>
struct opt_rep : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
};

template<class Field_type, Column_types column_type>
struct opt_rep_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine_b_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
};

template<class Field_type, Column_types column_type>
struct opt_vine_b_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
};

template<class Field_type, Column_types column_type>
struct opt_vine_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
};

template<class Field_type, Column_types column_type>
struct opt_vine_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_removable_columns = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
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
	boundaries.emplace_back();
	boundaries.push_back({1,7});
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
	boundaries.emplace_back();
	boundaries.push_back({{1,Field_type(1)},{7,Field_type(4)}});
}

typedef boost::mpl::list<Matrix<opt_bar_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_bar_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_b<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_b_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_bar_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_b<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_b_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_bar_b<Z2,Column_types::SET> >,
							Matrix<opt_bar_b<Z5,Column_types::SET> >,
							Matrix<opt_bar_b_r<Z2,Column_types::SET> >,
							Matrix<opt_bar_b_r<Z5,Column_types::SET> >,
							Matrix<opt_bar<Z2,Column_types::SET> >,
							Matrix<opt_bar<Z5,Column_types::SET> >,
							Matrix<opt_bar_r<Z2,Column_types::SET> >,
							Matrix<opt_bar_r<Z5,Column_types::SET> >,
							Matrix<opt_bar_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_bar_b<Z5,Column_types::VECTOR> >,
							Matrix<opt_bar_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_bar_b_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_bar<Z2,Column_types::VECTOR> >,
							Matrix<opt_bar<Z5,Column_types::VECTOR> >,
							Matrix<opt_bar_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_bar_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_bar_b<Z2,Column_types::LIST> >,
							Matrix<opt_bar_b<Z5,Column_types::LIST> >,
							Matrix<opt_bar_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_bar_b_r<Z5,Column_types::LIST> >,
							Matrix<opt_bar<Z2,Column_types::LIST> >,
							Matrix<opt_bar<Z5,Column_types::LIST> >,
							Matrix<opt_bar_r<Z2,Column_types::LIST> >,
							Matrix<opt_bar_r<Z5,Column_types::LIST> >,
							Matrix<opt_bar_b<Z2,Column_types::HEAP> >,
							Matrix<opt_bar_b_r<Z2,Column_types::HEAP> >
						> list_of_options_with_barcode_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Barcode_options, Matrix, list_of_options_with_barcode_access) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;
	
	struct BarComp {
		bool operator()(const std::tuple<int,int,int>& c1, const std::tuple<int,int,int>& c2) const
		{
			if (std::get<0>(c1) == std::get<0>(c2))
				return std::get<1>(c1) < std::get<1>(c2);
			return std::get<0>(c1) < std::get<0>(c2);
		}
	};

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	const auto& barcode = mb.get_current_barcode();
	std::set<std::tuple<int,int,int>,BarComp> bars;
	//bars are not ordered the same for base matrices
	for (auto it = barcode.begin(); it != barcode.end(); ++it){
		bars.emplace(it->dim, it->birth, it->death);
	}
	auto it = bars.begin();
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<2>(*it), -1);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 1);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 3);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 2);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 4);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 0);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 7);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 8);
	++it;
	BOOST_CHECK_EQUAL(std::get<0>(*it), 1);
	BOOST_CHECK_EQUAL(std::get<1>(*it), 5);
	BOOST_CHECK_EQUAL(std::get<2>(*it), 6);
	++it;
	BOOST_CHECK(it == bars.end());

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		BOOST_CHECK_EQUAL(mb.get_pivot(3), 1);
		BOOST_CHECK_EQUAL(mb.get_pivot(4), 2);
		BOOST_CHECK_EQUAL(mb.get_pivot(6), 5);
		BOOST_CHECK_EQUAL(mb.get_pivot(8), 7);
	} else {
		BOOST_CHECK_EQUAL(mb.get_pivot(0), 0);
		BOOST_CHECK_EQUAL(mb.get_pivot(1), 1);
		BOOST_CHECK_EQUAL(mb.get_pivot(2), 2);
		BOOST_CHECK_EQUAL(mb.get_pivot(3), 3);
		BOOST_CHECK_EQUAL(mb.get_pivot(4), 4);
		BOOST_CHECK_EQUAL(mb.get_pivot(5), 5);
		BOOST_CHECK_EQUAL(mb.get_pivot(6), 6);
		BOOST_CHECK_EQUAL(mb.get_pivot(7), 7);
		BOOST_CHECK_EQUAL(mb.get_pivot(8), 8);
	}
}

typedef boost::mpl::list<Matrix<opt_rep_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_rep_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_b<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_b_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_rep_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_b<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_b_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_rep_b<Z2,Column_types::SET> >,
							Matrix<opt_rep_b<Z5,Column_types::SET> >,
							Matrix<opt_rep_b_r<Z2,Column_types::SET> >,
							Matrix<opt_rep_b_r<Z5,Column_types::SET> >,
							Matrix<opt_rep<Z2,Column_types::SET> >,
							Matrix<opt_rep<Z5,Column_types::SET> >,
							Matrix<opt_rep_r<Z2,Column_types::SET> >,
							Matrix<opt_rep_r<Z5,Column_types::SET> >,
							Matrix<opt_rep_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_rep_b<Z5,Column_types::VECTOR> >,
							Matrix<opt_rep_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_rep_b_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_rep<Z2,Column_types::VECTOR> >,
							Matrix<opt_rep<Z5,Column_types::VECTOR> >,
							Matrix<opt_rep_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_rep_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_rep_b<Z2,Column_types::LIST> >,
							Matrix<opt_rep_b<Z5,Column_types::LIST> >,
							Matrix<opt_rep_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_rep_b_r<Z5,Column_types::LIST> >,
							Matrix<opt_rep<Z2,Column_types::LIST> >,
							Matrix<opt_rep<Z5,Column_types::LIST> >,
							Matrix<opt_rep_r<Z2,Column_types::LIST> >,
							Matrix<opt_rep_r<Z5,Column_types::LIST> >,
							Matrix<opt_rep_b<Z2,Column_types::HEAP> >,
							Matrix<opt_rep_b_r<Z2,Column_types::HEAP> >
						> list_of_options_with_rep_cycles;

BOOST_AUTO_TEST_CASE_TEMPLATE(Representative_cycle_options, Matrix, list_of_options_with_rep_cycles) {
	using boundary_matrix = typename std::conditional<
								Matrix::Option_list::is_z2,
								std::vector<std::vector<unsigned int> >,
								std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >
							>::type;

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

	if constexpr (Matrix::Option_list::is_of_boundary_type){
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
	} else {
		std::vector<unsigned int> tmp;
		tmp.push_back(0);
		BOOST_CHECK(cycles[0] == tmp);
		tmp.push_back(1);
		BOOST_CHECK(cycles[1] == tmp);
		tmp[1] = 2;
		BOOST_CHECK(cycles[2] == tmp);
		tmp[0] = 3;
		tmp[1] = 4;
		tmp.push_back(5);
		BOOST_CHECK(cycles[3] == tmp);
		tmp[0] = 0;
		tmp[1] = 7;
		tmp.pop_back();
		BOOST_CHECK(cycles[4] == tmp);
	}
}

typedef boost::mpl::list<Matrix<opt_vine_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_b<Z2,Column_types::SET> >,
							Matrix<opt_vine_b_r<Z2,Column_types::SET> >,
							Matrix<opt_vine_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_b<Z2,Column_types::LIST> >,
							Matrix<opt_vine_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_vine_b<Z2,Column_types::HEAP> >,
							Matrix<opt_vine_b_r<Z2,Column_types::HEAP> >,
							Matrix<opt_vine_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_ii<Z2,Column_types::SET> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::SET> >,
							Matrix<opt_vine_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_ii<Z2,Column_types::LIST> >,
							Matrix<opt_vine_r_ii<Z2,Column_types::LIST> >
						> list_of_options_with_vine_and_position_index;

BOOST_AUTO_TEST_CASE_TEMPLATE(Vine_option_with_position_indexing, Matrix, list_of_options_with_vine_and_position_index) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

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

typedef boost::mpl::list<Matrix<opt_vine_b_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::SET> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::SET> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::LIST> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::LIST> >,
							Matrix<opt_vine_b_ii<Z2,Column_types::HEAP> >,
							Matrix<opt_vine_b_r_ii<Z2,Column_types::HEAP> >,
							Matrix<opt_vine<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_vine<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_vine<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_vine<Z2,Column_types::SET> >,
							Matrix<opt_vine_r<Z2,Column_types::SET> >,
							Matrix<opt_vine<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_vine<Z2,Column_types::LIST> >,
							Matrix<opt_vine_r<Z2,Column_types::LIST> >
						> list_of_options_with_vine_and_id_index;

BOOST_AUTO_TEST_CASE_TEMPLATE(Vine_option_with_id_indexing, Matrix, list_of_options_with_vine_and_id_index) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

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
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(6), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(3));
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

	if constexpr (Matrix::Option_list::is_of_boundary_type){
		next = mb.vine_swap_with_z_eq_1_case(3, 8);		//use of simplex id
		BOOST_CHECK_EQUAL(next, 3);

		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(1), 5);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(2), 4);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(3), 8);
		BOOST_CHECK_EQUAL(mb.get_column_with_pivot(7), 6);
		BOOST_CHECK(mb.is_zero_column(0));
		BOOST_CHECK(mb.is_zero_column(1));
		BOOST_CHECK(mb.is_zero_column(2));
		BOOST_CHECK(mb.is_zero_column(7));
		BOOST_CHECK(mb.is_zero_column(3));
	} else {
		next = mb.vine_swap_with_z_eq_1_case(5, 8);		//use of internal chain id =/= initial simplex id
		BOOST_CHECK_EQUAL(next, 5);

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

template<class Field_type, Column_types column_type>
struct opt_ra_bar_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_bar_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_bar : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_bar_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_rep_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_rep_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_rep : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_rep_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_b_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_b_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_vine_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_bar_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_bar_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_bar : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_bar_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_rep_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_rep_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_rep : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_rep_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool can_retrieve_representative_cycles = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_b : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_b_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_r : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool is_indexed_by_position = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_b_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_b_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = true;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool is_indexed_by_position = false;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

template<class Field_type, Column_types column_type>
struct opt_ra_ni_vine_r_ii : Default_options<Field_type::get_characteristic() == 2, Field_type, column_type, false>{
	static const bool is_of_boundary_type = false;
	static const bool has_removable_columns = true;
	static const bool has_removable_rows = true;
	static const bool has_column_pairings = true;
	static const bool has_vine_update = true;
	static const bool has_row_access = true;
	static const bool has_intrusive_rows = false;
};

typedef boost::mpl::list<Matrix<opt_ra_bar_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_b<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar_b_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_bar_b<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_b_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_b<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_bar_b_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_bar<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_bar_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_bar_b<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_b_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_bar<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_b<Z5,Column_types::LIST> >,
							Matrix<opt_ra_bar_b_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_bar<Z5,Column_types::LIST> >,
							Matrix<opt_ra_bar_r<Z5,Column_types::LIST> >
						> list_of_options_with_barcode_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Barcode_and_row_access_options, Matrix, list_of_options_with_barcode_and_row_access) {
	using boundary_matrix = std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.get_current_barcode();

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		rows.push_back({{3,Z5(1)}});
		rows.push_back({{3,Z5(4)},{4,Z5(1)},{8,Z5(1)}});
		rows.push_back({{4,Z5(4)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{6,Z5(4)}});
		rows.push_back({});
		rows.push_back({{8,Z5(4)}});
	} else {
		rows.push_back({{0,Z5(1)},{1,Z5(1)},{2,Z5(1)},{7,Z5(1)}});
		rows.push_back({{1,Z5(4)}});
		rows.push_back({{2,Z5(4)}});
		rows.push_back({{3,Z5(1)},{4,Z5(1)},{5,Z5(1)},{8,Z5(1)}});
		rows.push_back({{4,Z5(1)},{5,Z5(1)}});
		rows.push_back({{5,Z5(4)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{7,Z5(4)}});
		rows.push_back({{8,Z5(1)}});
	}

	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}

typedef boost::mpl::list<Matrix<opt_ra_bar_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_bar_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_bar_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_bar_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_bar_b<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_bar_b_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_bar<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_bar_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_bar_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_bar<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_bar_b<Z2,Column_types::LIST> >,
							Matrix<opt_ra_bar_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_bar<Z2,Column_types::LIST> >,
							Matrix<opt_ra_bar_r<Z2,Column_types::LIST> >
						> list_of_z2_options_with_barcode_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Barcode_and_row_access_z2_options, Matrix, list_of_z2_options_with_barcode_and_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.get_current_barcode();

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		rows.push_back({3});
		rows.push_back({3,4,8});
		rows.push_back({4});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({});
		rows.push_back({8});
	} else {
		rows.push_back({0,1,2,7});
		rows.push_back({1});
		rows.push_back({2});
		rows.push_back({3,4,5,8});
		rows.push_back({4,5});
		rows.push_back({5});
		rows.push_back({6});
		rows.push_back({7});
		rows.push_back({8});
	}

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert(cell.get_column_index());
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}

typedef boost::mpl::list<Matrix<opt_ra_rep_b<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_b_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_r<Z5,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_b<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep_b_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep_r<Z5,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_rep_b<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_b_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_r<Z5,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_b<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_rep_b_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_rep<Z5,Column_types::SET> >,
							Matrix<opt_ra_ni_rep_r<Z5,Column_types::SET> >,
							Matrix<opt_ra_rep_b<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_b_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_rep<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_r<Z5,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_b<Z5,Column_types::LIST> >,
							Matrix<opt_ra_rep_b_r<Z5,Column_types::LIST> >,
							Matrix<opt_ra_rep<Z5,Column_types::LIST> >,
							Matrix<opt_ra_rep_r<Z5,Column_types::LIST> >
						> list_of_options_with_rep_cycles_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Representative_cycle_and_row_access_options, Matrix, list_of_options_with_rep_cycles_and_row_access) {
	using boundary_matrix = std::vector<std::vector<std::pair<unsigned int,typename Matrix::Field_type> > >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.update_representative_cycles();
	mb.get_representative_cycles();

	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		rows.push_back({{3,Z5(1)}});
		rows.push_back({{3,Z5(4)},{4,Z5(1)},{8,Z5(1)}});
		rows.push_back({{4,Z5(4)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{6,Z5(4)}});
		rows.push_back({});
		rows.push_back({{8,Z5(4)}});
	} else {
		rows.push_back({{0,Z5(1)},{1,Z5(1)},{2,Z5(1)},{7,Z5(1)}});
		rows.push_back({{1,Z5(4)}});
		rows.push_back({{2,Z5(4)}});
		rows.push_back({{3,Z5(1)},{4,Z5(1)},{5,Z5(1)},{8,Z5(1)}});
		rows.push_back({{4,Z5(1)},{5,Z5(1)}});
		rows.push_back({{5,Z5(4)}});
		rows.push_back({{6,Z5(1)}});
		rows.push_back({{7,Z5(4)}});
		rows.push_back({{8,Z5(1)}});
	}

	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}

typedef boost::mpl::list<Matrix<opt_ra_rep_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_rep_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_rep_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_rep_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_rep_b<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_rep_b_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_rep<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_rep_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_rep_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_rep<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_rep_b<Z2,Column_types::LIST> >,
							Matrix<opt_ra_rep_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_rep<Z2,Column_types::LIST> >,
							Matrix<opt_ra_rep_r<Z2,Column_types::LIST> >
						> list_of_z2_options_with_rep_cycles_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Representative_cycle_and_row_access_z2_options, Matrix, list_of_z2_options_with_rep_cycles_and_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.update_representative_cycles();
	mb.get_representative_cycles();

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		rows.push_back({3});
		rows.push_back({3,4,8});
		rows.push_back({4});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({});
		rows.push_back({8});
	} else {
		rows.push_back({0,1,2,7});
		rows.push_back({1});
		rows.push_back({2});
		rows.push_back({3,4,5,8});
		rows.push_back({4,5});
		rows.push_back({5});
		rows.push_back({6});
		rows.push_back({7});
		rows.push_back({8});
	}

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert(cell.get_column_index());
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}

typedef boost::mpl::list<Matrix<opt_ra_vine_b<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_b_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_b<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_vine_b_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_vine_b<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_b_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_b<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_vine_b_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_vine_b<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_b_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_b<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine_b_r<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_r_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_vine_r_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_vine_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_r_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_ii<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_vine_r_ii<Z2,Column_types::SET> >,
							Matrix<opt_ra_vine_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_r_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_ii<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine_r_ii<Z2,Column_types::LIST> >
						> list_of_options_with_vine_and_position_index_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Vine_option_with_position_indexing_and_row_access, Matrix, list_of_options_with_vine_and_position_index_and_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.vine_swap(6);
	mb.vine_swap(5);
	mb.vine_swap(4);
	mb.vine_swap(3);
	mb.vine_swap(7);
	mb.vine_swap(0);
	mb.vine_swap_with_z_eq_1_case(4);
	mb.vine_swap(5);
	mb.vine_swap_with_z_eq_1_case(6);

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
//		rows.push_back({5,6});
//		rows.push_back({4,5});
//		rows.push_back({4});
//		rows.push_back({6});
//		rows.push_back({8});
//		rows.push_back({8});
//		rows.push_back({});
//		rows.push_back({8});
//		rows.push_back({});
		rows.push_back({5,8});
		rows.push_back({4,5});
		rows.push_back({4});
		rows.push_back({8});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({});
		rows.push_back({6});
	} else {
		rows.push_back({0,1,2});
		rows.push_back({1,7});
		rows.push_back({2});
		rows.push_back({7});
		rows.push_back({3,4,5,8});
		rows.push_back({3,5,8});
		rows.push_back({8});
		rows.push_back({5});
		rows.push_back({6});
	}

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert(cell.get_column_index());
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}

typedef boost::mpl::list<Matrix<opt_ra_vine_b_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_b_r_ii<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_b_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_vine_b_r_ii<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_vine_b_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_b_r_ii<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_b_ii<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_vine_b_r_ii<Z2,Column_types::SET> >,
							Matrix<opt_ra_vine_b_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_b_r_ii<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_b_ii<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine_b_r_ii<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine_r<Z2,Column_types::INTRUSIVE_SET> >,
							Matrix<opt_ra_vine<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_vine_r<Z2,Column_types::INTRUSIVE_LIST> >,
							Matrix<opt_ra_ni_vine<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine_r<Z2,Column_types::UNORDERED_SET> >,
							Matrix<opt_ra_ni_vine<Z2,Column_types::SET> >,
							Matrix<opt_ra_ni_vine_r<Z2,Column_types::SET> >,
							Matrix<opt_ra_vine<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine_r<Z2,Column_types::VECTOR> >,
							Matrix<opt_ra_vine<Z2,Column_types::LIST> >,
							Matrix<opt_ra_vine_r<Z2,Column_types::LIST> >
						> list_of_options_with_vine_and_id_index_and_row_access;

BOOST_AUTO_TEST_CASE_TEMPLATE(Vine_option_with_id_indexing_and_row_access, Matrix, list_of_options_with_vine_and_id_index_and_row_access) {
	using boundary_matrix = std::vector<std::vector<unsigned int> >;

	boundary_matrix ordered_boundaries;
	build_boundary_matrix(ordered_boundaries);

	Matrix mb(ordered_boundaries);

	mb.vine_swap(6, 7);
	mb.vine_swap(5, 7);
	mb.vine_swap(4, 7);
	mb.vine_swap(3, 7);
	mb.vine_swap(6, 8);
	mb.vine_swap(0, 1);
	mb.vine_swap_with_z_eq_1_case(3, 4);
	mb.vine_swap(3, 5);

	std::vector<std::vector<unsigned int> > rows;
	if constexpr (Matrix::Option_list::is_of_boundary_type){
		mb.vine_swap_with_z_eq_1_case(3, 8);	//use of simplex id

		rows.push_back({4,5});
		rows.push_back({5,8});
		rows.push_back({4});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({6});
		rows.push_back({});
		rows.push_back({8});
	} else {
		mb.vine_swap_with_z_eq_1_case(5, 8);	//use of internal chain id =/= initial simplex id

		rows.push_back({1,7});
		rows.push_back({0,1,2});
		rows.push_back({2});
		rows.push_back({5});
		rows.push_back({3,4,5,8});
		rows.push_back({3,5,8});
		rows.push_back({6});
		rows.push_back({7});
		rows.push_back({8});
	}

	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(rows.size());
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
			for (auto& cell : mb.get_row(i)){
				ordered_rows[i].insert(cell.get_column_index());
			}
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		if (!Matrix::Option_list::is_of_boundary_type || !Matrix::Option_list::has_removable_rows || i != 6){
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
}


