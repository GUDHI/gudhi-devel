/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_TESTS_TYPE_LISTS_H
#define PM_COLUMN_TESTS_TYPE_LISTS_H

#include <boost/mp11.hpp>

#include "pm_column_tests_mastermatrix.h"
#include "pm_common_boost_type_lists.h"
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistence_matrix/columns/intrusive_list_column.h>
#include <gudhi/Persistence_matrix/columns/intrusive_set_column.h>
#include <gudhi/Persistence_matrix/columns/list_column.h>
#include <gudhi/Persistence_matrix/columns/set_column.h>
#include <gudhi/Persistence_matrix/columns/unordered_set_column.h>
#include <gudhi/Persistence_matrix/columns/vector_column.h>
#include <gudhi/Persistence_matrix/columns/naive_vector_column.h>
#include <gudhi/Persistence_matrix/columns/heap_column.h>

using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Intrusive_list_column;
using Gudhi::persistence_matrix::Intrusive_set_column;
using Gudhi::persistence_matrix::List_column;
using Gudhi::persistence_matrix::Set_column;
using Gudhi::persistence_matrix::Unordered_set_column;
using Gudhi::persistence_matrix::Vector_column;
using Gudhi::persistence_matrix::Naive_vector_column;
using Gudhi::persistence_matrix::Heap_column;

template<typename is_z2_only, typename col_type, typename has_row, typename rem_row, typename intr_row>
struct c_base_options{
	using type = Base_col_options<is_z2_only::t,col_type::t,has_row::t,rem_row::t,intr_row::t>;
};

template<typename is_z2_only, typename col_type, typename has_row, typename rem_row, typename intr_row>
struct c_boundary_options{
	using type = Boundary_col_options<is_z2_only::t,col_type::t,has_row::t,rem_row::t,intr_row::t>;
};

template<typename is_z2_only, typename col_type, typename has_row, typename rem_row, typename intr_row>
struct c_chain_options{
	using type = Chain_col_options<is_z2_only::t,col_type::t,has_row::t,rem_row::t,intr_row::t>;
};

template <class col_type>
class column_non_validity{
private:
	static constexpr bool is_non_valide(){
		if constexpr (col_type::Master::Option_list::column_type == Column_types::INTRUSIVE_LIST){
			return !std::is_same_v<col_type, Intrusive_list_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::INTRUSIVE_SET){
			return !std::is_same_v<col_type, Intrusive_set_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::LIST){
			return !std::is_same_v<col_type, List_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::SET){
			return !std::is_same_v<col_type, Set_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::UNORDERED_SET){
			return !std::is_same_v<col_type, Unordered_set_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::NAIVE_VECTOR){
			return !std::is_same_v<col_type, Naive_vector_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::VECTOR){
			return !std::is_same_v<col_type, Vector_column<typename col_type::Master> >;
		} else if constexpr (col_type::Master::Option_list::column_type == Column_types::HEAP){
			return !std::is_same_v<col_type, Heap_column<typename col_type::Master> >;
		} else {
			return true;	//we should not enter here, except if we want to ignore a column type.
		}
	}
public:
	static constexpr bool value = is_non_valide();
};

//if a new column type is implemented, create a `ct_*` structure for it and add it to this list...
using col_type_list = boost::mp11::mp_list<ct_intrusive_list, ct_intrusive_set, ct_list, ct_set, ct_heap, ct_unordered_set, ct_vector, ct_naive_vector>;
using row_col_type_list = boost::mp11::mp_list<ct_intrusive_list, ct_intrusive_set, ct_list, ct_set, ct_unordered_set, ct_vector, ct_naive_vector>;
//...and add the column name here.
using column_list = mp_list_q<Intrusive_list_column,Intrusive_set_column,List_column,Set_column,Unordered_set_column,Naive_vector_column,Vector_column,Heap_column>;
using c_matrix_type_list = mp_list_q<Column_mini_matrix>;

template<typename option_name_list, typename bool_is_z2, typename col_t, typename bool_has_row, typename bool_rem_row, typename bool_intr_row>
using option_template = boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, option_name_list, bool_is_z2, col_t, bool_has_row, bool_rem_row, bool_intr_row> >;

template<typename option_name_list>
using z2_no_ra_option_list = option_template<option_name_list, true_value_list, col_type_list, false_value_list, false_value_list, false_value_list>;
template<typename option_name_list>
using z2_only_ra_r_option_list = option_template<option_name_list, true_value_list, row_col_type_list, true_value_list, true_value_list, bool_value_list>;
template<typename option_name_list>
using z2_only_ra_option_list = option_template<option_name_list, true_value_list, row_col_type_list, true_value_list, false_value_list, bool_value_list>;
template<typename option_name_list>
using z5_no_ra_option_list = option_template<option_name_list, false_value_list, col_type_list, false_value_list, false_value_list, false_value_list>;
template<typename option_name_list>
using z5_only_ra_r_option_list = option_template<option_name_list, false_value_list, row_col_type_list, true_value_list, true_value_list, bool_value_list>;
template<typename option_name_list>
using z5_only_ra_option_list = option_template<option_name_list, false_value_list, row_col_type_list, true_value_list, false_value_list, bool_value_list>;

template<typename complete_option_list>
using c_matrices_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, c_matrix_type_list, complete_option_list>;

//the remove_if here is quite unefficiant as it has to remove a lot. But I did not found another way to do it without having to define something for each column type individually.
template<typename complete_option_list>
using columns_list = boost::mp11::mp_remove_if<boost::mp11::mp_product<boost::mp11::mp_invoke_q, column_list, c_matrices_list<complete_option_list> >, column_non_validity>;

#endif // PM_COLUMN_TESTS_TYPE_LISTS_H
