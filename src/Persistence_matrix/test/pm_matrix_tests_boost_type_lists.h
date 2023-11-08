/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_MATRIX_TESTS_TYPE_LISTS_H
#define PM_MATRIX_TESTS_TYPE_LISTS_H

#include <boost/mp11.hpp>
#include <type_traits>

#include "pm_matrix_tests_options.h"
#include "pm_common_boost_type_lists.h"
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/matrix.h>

using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Matrix;

template<bool ra_v, bool rem_row_v, bool intr_row_v>
struct ra_value{
	static constexpr const bool ra = ra_v;
	static constexpr const bool rem_row = rem_row_v;
	static constexpr const bool intr_row = intr_row_v;
};

template<bool z2_v, bool vine_v, bool rep_v>
struct ru_opt_values{
	static constexpr const bool vine = vine_v;
	static constexpr const bool rep = rep_v;
	static constexpr const bool z2 = z2_v;
};

template<typename z2_v, typename barcode_v, typename vine_v, typename rep_v>
struct chain_opt_values{
	static constexpr const bool vine = vine_v::t;
	static constexpr const bool rep = rep_v::t;
	static constexpr const bool barcode = barcode_v::t;
	static constexpr const bool z2 = z2_v::t;
};

template<typename is_z2_only, typename col_type, typename ra_value, typename rem_col, typename swaps>
struct m_base_options{
	using type = typename std::conditional<
								ra_value::ra, 
								Base_options_with_row_access<is_z2_only::t,col_type::t,ra_value::rem_row,ra_value::intr_row,rem_col::t,swaps::t>, 
								Base_options<is_z2_only::t,col_type::t,rem_col::t,swaps::t> 
							>::type;
};

template<typename is_z2_only, typename col_type, typename ra_value>
struct m_col_comp_options{
	using type = typename std::conditional<
								ra_value::ra, 
								Column_compression_options_with_row_access<is_z2_only::t,col_type::t,ra_value::rem_row,ra_value::intr_row>, 
								Column_compression_options<is_z2_only::t,col_type::t> 
							>::type;
};

template<typename is_z2_only, typename col_type, typename ra_value, typename rem_col, typename swaps, typename pos_idx>
struct m_boundary_options{
	using type = typename std::conditional<
								ra_value::ra, 
								Boundary_options_with_row_access<is_z2_only::t,col_type::t,ra_value::rem_row,ra_value::intr_row,rem_col::t,swaps::t,pos_idx::t>, 
								Boundary_options<is_z2_only::t,col_type::t,rem_col::t,swaps::t,pos_idx::t> 
							>::type;
};

template<typename col_type, typename ru_opt_values, typename ra_value, typename rem_col, typename pos_idx, typename dim, typename barcode>
struct m_ru_options{
	using type = typename std::conditional<
								ra_value::ra, 
								typename std::conditional<
									ru_opt_values::vine, 
									RU_vine_options_with_row_access<col_type::t,ru_opt_values::rep,ra_value::rem_row,ra_value::intr_row,rem_col::t,pos_idx::t,dim::t,barcode::t>, 
									RU_rep_options_with_row_access<ru_opt_values::z2,col_type::t,ra_value::rem_row,ra_value::intr_row,rem_col::t,pos_idx::t,dim::t,barcode::t> 
								>::type, 
								typename std::conditional<
									ru_opt_values::vine, 
									RU_vine_options<col_type::t,ru_opt_values::rep,rem_col::t,pos_idx::t,dim::t,barcode::t>, 
									RU_rep_options<ru_opt_values::z2,col_type::t,rem_col::t,pos_idx::t,dim::t,barcode::t> 
								>::type
							>::type;
};

template<typename col_type, typename chain_opt_values, typename ra_value, typename rem_col, typename pos_idx, typename dim>
struct m_chain_options{
	using type = typename std::conditional<
								ra_value::ra, 
								typename std::conditional<
									chain_opt_values::vine, 
									Chain_vine_options_with_row_access<col_type::t,chain_opt_values::rep,chain_opt_values::barcode,ra_value::rem_row,ra_value::intr_row,rem_col::t,pos_idx::t,dim::t>, 
									typename std::conditional<
										chain_opt_values::rep, 
										Chain_rep_options_with_row_access<chain_opt_values::z2,col_type::t,chain_opt_values::barcode,ra_value::rem_row,ra_value::intr_row,rem_col::t,pos_idx::t,dim::t>, 
										Chain_barcode_options_with_row_access<chain_opt_values::z2,col_type::t,ra_value::rem_row,ra_value::intr_row,rem_col::t,pos_idx::t,dim::t> 
									>::type
								>::type, 
								typename std::conditional<
									chain_opt_values::vine, 
									Chain_vine_options<col_type::t,chain_opt_values::rep,chain_opt_values::barcode,rem_col::t,pos_idx::t,dim::t>,
									typename std::conditional<
										chain_opt_values::rep, 
										Chain_rep_options<chain_opt_values::z2,col_type::t,chain_opt_values::barcode,rem_col::t,pos_idx::t,dim::t>, 
										Chain_barcode_options<chain_opt_values::z2,col_type::t,rem_col::t,pos_idx::t,dim::t> 
									>::type
								>::type
							>::type;
};

template <class option>
class matrix_non_validity{
private:
	static constexpr bool is_non_valide(){
		return option::has_row_access && option::column_type == Column_types::HEAP;
	}
public:
	static constexpr bool value = is_non_valide();
};

template <class option>
class chain_opt_non_validity{
private:
	static constexpr bool is_non_valide(){
		return (!option::vine && !option::rep && !option::barcode) || (option::vine && !option::z2);
	}
public:
	static constexpr bool value = is_non_valide();
};

using matrix_type_list = mp_list_q<Matrix>;
using col_type_list = boost::mp11::mp_list<ct_intrusive_list>;	//to avoid long compilation time and high memory usage, the tests are restricted to one column type by default for users. But in case real changes were made to this module, it would be better to test at least once with all column types as below:
// using col_type_list = boost::mp11::mp_list<ct_intrusive_list, ct_intrusive_set, ct_list, ct_set, ct_heap, ct_unordered_set, ct_vector, ct_naive_vector>;

using all_ra_values_list = boost::mp11::mp_list<ra_value<false,false,false>,ra_value<true,false,false>,ra_value<true,true,true>,ra_value<true,true,false>,ra_value<true,false,true> >;
using ra_values_list = boost::mp11::mp_list<ra_value<true,false,false>,ra_value<true,true,true>,ra_value<true,true,false>,ra_value<true,false,true> >;
using ra_r_values_list = boost::mp11::mp_list<ra_value<true,true,true>,ra_value<true,true,false> >;

using z2_ru_vine_values_list = boost::mp11::mp_list<ru_opt_values<true,true,false>,ru_opt_values<true,true,true> >;
using z2_ru_vine_rep_values_list = boost::mp11::mp_list<ru_opt_values<true,true,true> >;
using z2_ru_rep_values_list = boost::mp11::mp_list<ru_opt_values<true,false,true> >;
using zp_ru_rep_values_list = boost::mp11::mp_list<ru_opt_values<false,false,true> >;

// using z2_chain_vine_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, bool_value_list, true_value_list, bool_value_list>;
using z2_chain_vine_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, true_value_list, true_value_list>;
using z2_chain_vine_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, true_value_list, bool_value_list>;
using z2_chain_vine_no_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, false_value_list, true_value_list, bool_value_list>;
using z2_chain_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, bool_value_list, false_value_list, true_value_list>;
using z2_chain_rep_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, false_value_list, true_value_list>;
using z2_chain_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, false_value_list, false_value_list>;
using zp_chain_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, bool_value_list, false_value_list, true_value_list>;
using zp_chain_rep_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, true_value_list, false_value_list, true_value_list>;
using zp_chain_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, true_value_list, false_value_list, false_value_list>;

using z2_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, true_value_list, col_type_list, all_ra_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, true_value_list, col_type_list, all_ra_values_list> >, matrix_non_validity>;
using z2_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, true_value_list, col_type_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, false_value_list, col_type_list, all_ra_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, false_value_list, col_type_list, all_ra_values_list> >, matrix_non_validity>;
using zp_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, false_value_list, col_type_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

//TODO: test if the compilation is faster when using mp_remove_if on the full_*_option_list instead of recalculating them with other templates.
// or computing all smaller partitions first to concatenate them for the bigger ones (mp_append)

using z2_ra_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, true_value_list, col_type_list, ra_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, true_value_list, col_type_list, ra_values_list> >, matrix_non_validity>;
using z2_ra_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, true_value_list, col_type_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, false_value_list, col_type_list, ra_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, false_value_list, col_type_list, ra_values_list> >, matrix_non_validity>;
using zp_ra_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, false_value_list, col_type_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_ra_r_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, true_value_list, col_type_list, ra_r_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, true_value_list, col_type_list, ra_r_values_list> >, matrix_non_validity>;
using z2_ra_r_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, true_value_list, col_type_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_r_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, false_value_list, col_type_list, ra_r_values_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_r_col_comp_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, false_value_list, col_type_list, ra_r_values_list> >, matrix_non_validity>;
using zp_ra_r_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, false_value_list, col_type_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_r_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, true_value_list, col_type_list, all_ra_values_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, true_value_list, col_type_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_r_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, false_value_list, col_type_list, all_ra_values_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using zp_r_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, false_value_list, col_type_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_dim_boundary_option_list = z2_boundary_option_list;
using zp_dim_boundary_option_list = zp_boundary_option_list;

using z2_barcode_boundary_option_list = z2_boundary_option_list;
using zp_barcode_boundary_option_list = zp_boundary_option_list;

using z2_swap_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, true_value_list, col_type_list, all_ra_values_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using z2_swap_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, true_value_list, col_type_list, all_ra_values_list, bool_value_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using zp_swap_base_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, false_value_list, col_type_list, all_ra_values_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using zp_swap_boundary_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, false_value_list, col_type_list, all_ra_values_list, bool_value_list, true_value_list, bool_value_list> >, matrix_non_validity>;

using z2_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_dim_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using z2_barcode_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using z2_rep_ru_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_vine_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_r_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_r_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_dim_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using zp_dim_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list, bool_value_list> >, matrix_non_validity>;
using z2_barcode_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, z2_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using zp_barcode_ru_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, zp_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;

// using z2_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_no_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_barcode_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_barcode_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_barcode_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_dim_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using z2_barcode_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_no_barcode_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_no_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_rep_chain_vine_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_vine_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_r_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_r_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_dim_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using zp_dim_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using z2_barcode_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_rep_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_barcode_chain_rep_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_rep_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;

using z2_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_barcode_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_barcode_values_list, ra_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_ra_r_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_barcode_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_ra_r_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_barcode_values_list, ra_r_values_list, bool_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_r_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_barcode_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using zp_r_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_barcode_values_list, all_ra_values_list, true_value_list, bool_value_list, bool_value_list> >, matrix_non_validity>;
using z2_dim_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, z2_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;
using zp_dim_chain_bar_option_list = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, zp_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_value_list, true_value_list> >, matrix_non_validity>;

template<typename complete_option_list>
using matrices_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, matrix_type_list, complete_option_list>;

#endif // PM_MATRIX_TESTS_TYPE_LISTS_H
