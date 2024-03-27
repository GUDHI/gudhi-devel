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
		return (option::has_row_access || option::has_column_compression) && option::column_type == Column_types::HEAP;
	}
public:
	static constexpr bool value = is_non_valide();
};

//to avoid long compilation time and high memory usage, the tests are restricted to one column type by default for users. But in case real changes were made to this module, it would be better to test at least once with all column types as below:
// using col_type_list = boost::mp11::mp_list<ct_intrusive_list, ct_intrusive_set, ct_list, ct_set, ct_heap, ct_unordered_set, ct_vector, ct_naive_vector>;
#ifdef PM_TEST_INTR_LIST
using col_type_list = boost::mp11::mp_list<ct_intrusive_list>;
#else
#ifdef PM_TEST_INTR_SET
using col_type_list = boost::mp11::mp_list<ct_intrusive_set>;
#else
#ifdef PM_TEST_LIST
using col_type_list = boost::mp11::mp_list<ct_list>;
#else
#ifdef PM_TEST_SET
using col_type_list = boost::mp11::mp_list<ct_set>;
#else
#ifdef PM_TEST_HEAP
using col_type_list = boost::mp11::mp_list<ct_heap>;	//WARNING: unit tests involving row access will not compile (they template list will be empty), so they have to be commented to test heap columns alone
#else
#ifdef PM_TEST_UNORD_SET
using col_type_list = boost::mp11::mp_list<ct_unordered_set>;
#else
#ifdef PM_TEST_NAIVE_VECTOR
using col_type_list = boost::mp11::mp_list<ct_naive_vector>;
#else
using col_type_list = boost::mp11::mp_list<ct_vector>;
#endif
#endif
#endif
#endif
#endif
#endif
#endif

using matrix_type_list = mp_list_q<Matrix>;

using all_ra_values_list = boost::mp11::mp_list<ra_value<false,false,false>,ra_value<true,false,false>,ra_value<true,true,true>,ra_value<true,true,false>,ra_value<true,false,true> >;
using ra_values_list = boost::mp11::mp_list<ra_value<true,false,false>,ra_value<true,true,true>,ra_value<true,true,false>,ra_value<true,false,true> >;
using ra_r_values_list = boost::mp11::mp_list<ra_value<true,true,true>,ra_value<true,true,false> >;

// Base matrices

template<typename bool_is_z2, typename has_row_t, typename bool_rem_col, typename bool_swap>
using base_option_template = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_base_options>, bool_is_z2, col_type_list, has_row_t, bool_rem_col, bool_swap> >, matrix_non_validity>;

using opt_base_z2 = base_option_template<true_value_list, all_ra_values_list, bool_value_list, bool_value_list>;
using opt_base_z2_ra = base_option_template<true_value_list, ra_values_list, bool_value_list, bool_value_list>;
using opt_base_z2_ra_r = base_option_template<true_value_list, ra_r_values_list, bool_value_list, bool_value_list>;
using opt_base_z2_r = base_option_template<true_value_list, all_ra_values_list, true_value_list, bool_value_list>;
using opt_base_z2_swap = base_option_template<true_value_list, all_ra_values_list, bool_value_list, true_value_list>;

using opt_base_zp = base_option_template<false_value_list, all_ra_values_list, bool_value_list, bool_value_list>;
using opt_base_zp_ra = base_option_template<false_value_list, ra_values_list, bool_value_list, bool_value_list>;
using opt_base_zp_ra_r = base_option_template<false_value_list, ra_r_values_list, bool_value_list, bool_value_list>;
using opt_base_zp_r = base_option_template<false_value_list, all_ra_values_list, true_value_list, bool_value_list>;
using opt_base_zp_swap = base_option_template<false_value_list, all_ra_values_list, bool_value_list, true_value_list>;

// Compression matrices

template<typename bool_is_z2, typename has_row_t>
using compression_option_template = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_col_comp_options>, bool_is_z2, col_type_list, has_row_t> >, matrix_non_validity>;

using opt_col_comp_z2 = compression_option_template<true_value_list, all_ra_values_list>;
using opt_col_comp_z2_ra = compression_option_template<true_value_list, ra_values_list>;
using opt_col_comp_z2_ra_r = compression_option_template<true_value_list, ra_r_values_list>;

using opt_col_comp_zp = compression_option_template<false_value_list, all_ra_values_list>;
using opt_col_comp_zp_ra = compression_option_template<false_value_list, ra_values_list>;
using opt_col_comp_zp_ra_r = compression_option_template<false_value_list, ra_r_values_list>;

// Boundary matrices

template<typename bool_is_z2, typename has_row_t, typename bool_rem_col, typename bool_swap, typename bool_pos_idx>
using boundary_option_template = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_boundary_options>, bool_is_z2, col_type_list, has_row_t, bool_rem_col, bool_swap, bool_pos_idx> >, matrix_non_validity>;

template<typename bool_pos_idx>
using opt_boundary_z2 = boundary_option_template<true_value_list, all_ra_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_ra = boundary_option_template<true_value_list, ra_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_ra_r = boundary_option_template<true_value_list, ra_r_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_r = boundary_option_template<true_value_list, all_ra_values_list, true_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_dim = opt_boundary_z2<bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_barcode = opt_boundary_z2<bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_z2_swap = boundary_option_template<true_value_list, all_ra_values_list, bool_value_list, true_value_list, bool_pos_idx>;

template<typename bool_pos_idx>
using opt_boundary_zp = boundary_option_template<false_value_list, all_ra_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_ra = boundary_option_template<false_value_list, ra_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_ra_r = boundary_option_template<false_value_list, ra_r_values_list, bool_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_r = boundary_option_template<false_value_list, all_ra_values_list, true_value_list, bool_value_list, bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_dim = opt_boundary_zp<bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_barcode = opt_boundary_zp<bool_pos_idx>;
template<typename bool_pos_idx>
using opt_boundary_zp_swap = boundary_option_template<false_value_list, all_ra_values_list, bool_value_list, true_value_list, bool_pos_idx>;

// RU matrices

using z2_ru_vine_values_list = boost::mp11::mp_list<ru_opt_values<true,true,false>,ru_opt_values<true,true,true> >;
using z2_ru_vine_rep_values_list = boost::mp11::mp_list<ru_opt_values<true,true,true> >;
using z2_ru_rep_values_list = boost::mp11::mp_list<ru_opt_values<true,false,true> >;
using zp_ru_rep_values_list = boost::mp11::mp_list<ru_opt_values<false,false,true> >;

template<typename ru_opt_t, typename has_row_t, typename bool_rem_col, typename bool_pos_idx, typename bool_dim, typename bool_barcode>
using ru_option_template = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_ru_options>, col_type_list, ru_opt_t, has_row_t, bool_rem_col, bool_pos_idx, bool_dim, bool_barcode> >, matrix_non_validity>;

template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_vine_z2 = ru_option_template<z2_ru_vine_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_vine_z2_ra = ru_option_template<z2_ru_vine_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_vine_z2_ra_r = ru_option_template<z2_ru_vine_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_vine_z2_r = ru_option_template<z2_ru_vine_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_vine_z2_rep = ru_option_template<z2_ru_vine_rep_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;

template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_z2 = ru_option_template<z2_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_z2_ra = ru_option_template<z2_ru_rep_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_z2_ra_r = ru_option_template<z2_ru_rep_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_z2_r = ru_option_template<z2_ru_rep_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim, bool_barcode>;

template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_zp = ru_option_template<zp_ru_rep_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_zp_ra = ru_option_template<zp_ru_rep_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_zp_ra_r = ru_option_template<zp_ru_rep_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim, bool_barcode>;
template<typename bool_pos_idx, typename bool_barcode, typename bool_dim>
using opt_ru_rep_zp_r = ru_option_template<zp_ru_rep_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim, bool_barcode>;

// Chain matrices

using z2_chain_vine_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, true_value_list, true_value_list>;
using z2_chain_vine_rep_no_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, false_value_list, true_value_list, true_value_list>;
using z2_chain_vine_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, true_value_list, bool_value_list>;
using z2_chain_vine_no_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, false_value_list, true_value_list, bool_value_list>;
using z2_chain_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, bool_value_list, false_value_list, true_value_list>;
using z2_chain_rep_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, false_value_list, true_value_list>;
using zp_chain_rep_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, bool_value_list, false_value_list, true_value_list>;
using zp_chain_rep_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, true_value_list, false_value_list, true_value_list>;
using z2_chain_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, true_value_list, true_value_list, false_value_list, false_value_list>;
using zp_chain_barcode_values_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<chain_opt_values>, false_value_list, true_value_list, false_value_list, false_value_list>;

template<typename chain_opt_t, typename has_row_t, typename bool_rem_col, typename bool_pos_idx, typename bool_dim>
using chain_option_template = boost::mp11::mp_remove_if<boost::mp11::mp_transform<get_type, boost::mp11::mp_product<boost::mp11::mp_invoke_q, mp_list_q<m_chain_options>, col_type_list, chain_opt_t, has_row_t, bool_rem_col, bool_pos_idx, bool_dim> >, matrix_non_validity>;

template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_ra = chain_option_template<z2_chain_vine_barcode_values_list, ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_ra_no_barcode = chain_option_template<z2_chain_vine_no_barcode_values_list, ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_ra_r = chain_option_template<z2_chain_vine_barcode_values_list, ra_r_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_ra_r_no_barcode = chain_option_template<z2_chain_vine_no_barcode_values_list, ra_r_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_barcode = chain_option_template<z2_chain_vine_barcode_values_list, all_ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_no_barcode = chain_option_template<z2_chain_vine_no_barcode_values_list, all_ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_rep = chain_option_template<z2_chain_vine_rep_values_list, all_ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_rem_col, typename bool_dim>
using opt_chain_vine_z2_rep_no_barcode = chain_option_template<z2_chain_vine_rep_no_barcode_values_list, all_ra_values_list, bool_rem_col, bool_pos_idx, bool_dim>;

template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_z2 = chain_option_template<z2_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_z2_ra = chain_option_template<z2_chain_rep_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_z2_ra_r = chain_option_template<z2_chain_rep_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_z2_r = chain_option_template<z2_chain_rep_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_z2_barcode = chain_option_template<z2_chain_rep_barcode_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;

template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_zp = chain_option_template<zp_chain_rep_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_zp_ra = chain_option_template<zp_chain_rep_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_zp_ra_r = chain_option_template<zp_chain_rep_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_zp_r = chain_option_template<zp_chain_rep_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_rep_zp_barcode = chain_option_template<zp_chain_rep_barcode_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;

template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_z2 = chain_option_template<z2_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_z2_ra = chain_option_template<z2_chain_barcode_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_z2_ra_r = chain_option_template<z2_chain_barcode_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_z2_r = chain_option_template<z2_chain_barcode_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim>;

template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_zp = chain_option_template<zp_chain_barcode_values_list, all_ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_zp_ra = chain_option_template<zp_chain_barcode_values_list, ra_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_zp_ra_r = chain_option_template<zp_chain_barcode_values_list, ra_r_values_list, bool_value_list, bool_pos_idx, bool_dim>;
template<typename bool_pos_idx, typename bool_dim>
using opt_chain_bar_zp_r = chain_option_template<zp_chain_barcode_values_list, all_ra_values_list, true_value_list, bool_pos_idx, bool_dim>;

// Final template

template<typename complete_option_list>
using matrices_list = boost::mp11::mp_product<boost::mp11::mp_invoke_q, matrix_type_list, complete_option_list >;

#endif // PM_MATRIX_TESTS_TYPE_LISTS_H
