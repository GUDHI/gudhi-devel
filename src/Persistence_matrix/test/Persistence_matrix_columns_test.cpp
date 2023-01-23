/*	This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
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
#include "gudhi/column_types/row_access.h"
#include "gudhi/column_types/cell.h"

#include "gudhi/column_types/list_column.h"
#include "gudhi/column_types/set_column.h"
#include "gudhi/column_types/unordered_set_column.h"
#include "gudhi/column_types/vector_column.h"
#include "gudhi/column_types/z2_heap_column.h"
#include "gudhi/column_types/z2_list_column.h"
#include "gudhi/column_types/z2_set_column.h"
#include "gudhi/column_types/z2_unordered_set_column.h"
#include "gudhi/column_types/z2_vector_column.h"
#include "gudhi/column_types/intrusive_list_column.h"
#include "gudhi/column_types/intrusive_set_column.h"
#include "gudhi/column_types/z2_intrusive_list_column.h"
#include "gudhi/column_types/z2_intrusive_set_column.h"

#include "gudhi/column_types/boundary_columns/boundary_list_column.h"
#include "gudhi/column_types/boundary_columns/boundary_set_column.h"
#include "gudhi/column_types/boundary_columns/boundary_unordered_set_column.h"
#include "gudhi/column_types/boundary_columns/boundary_vector_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_heap_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_list_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_set_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_unordered_set_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_vector_column.h"
#include "gudhi/column_types/boundary_columns/boundary_intrusive_list_column.h"
#include "gudhi/column_types/boundary_columns/boundary_intrusive_set_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_intrusive_list_column.h"
#include "gudhi/column_types/boundary_columns/z2_boundary_intrusive_set_column.h"

#include "gudhi/column_types/chain_columns/chain_list_column.h"
#include "gudhi/column_types/chain_columns/chain_set_column.h"
#include "gudhi/column_types/chain_columns/chain_unordered_set_column.h"
#include "gudhi/column_types/chain_columns/chain_vector_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_heap_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_list_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_set_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_unordered_set_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_vector_column.h"
#include "gudhi/column_types/chain_columns/chain_intrusive_list_column.h"
#include "gudhi/column_types/chain_columns/chain_intrusive_set_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_intrusive_list_column.h"
#include "gudhi/column_types/chain_columns/z2_chain_intrusive_set_column.h"

using Gudhi::persistence_matrix::Zp_field_element;

using Gudhi::persistence_matrix::List_column;
using Gudhi::persistence_matrix::Set_column;
using Gudhi::persistence_matrix::Unordered_set_column;
using Gudhi::persistence_matrix::Vector_column;
using Gudhi::persistence_matrix::Z2_heap_column;
using Gudhi::persistence_matrix::Z2_list_column;
using Gudhi::persistence_matrix::Z2_set_column;
using Gudhi::persistence_matrix::Z2_unordered_set_column;
using Gudhi::persistence_matrix::Z2_vector_column;
using Gudhi::persistence_matrix::Intrusive_list_column;
using Gudhi::persistence_matrix::Intrusive_set_column;
using Gudhi::persistence_matrix::Z2_intrusive_list_column;
using Gudhi::persistence_matrix::Z2_intrusive_set_column;

using Gudhi::persistence_matrix::List_boundary_column;
using Gudhi::persistence_matrix::Set_boundary_column;
using Gudhi::persistence_matrix::Unordered_set_boundary_column;
using Gudhi::persistence_matrix::Vector_boundary_column;
using Gudhi::persistence_matrix::Z2_heap_boundary_column;
using Gudhi::persistence_matrix::Z2_list_boundary_column;
using Gudhi::persistence_matrix::Z2_set_boundary_column;
using Gudhi::persistence_matrix::Z2_unordered_set_boundary_column;
using Gudhi::persistence_matrix::Z2_vector_boundary_column;
using Gudhi::persistence_matrix::Intrusive_list_boundary_column;
using Gudhi::persistence_matrix::Intrusive_set_boundary_column;
using Gudhi::persistence_matrix::Z2_intrusive_list_boundary_column;
using Gudhi::persistence_matrix::Z2_intrusive_set_boundary_column;

using Gudhi::persistence_matrix::List_chain_column;
using Gudhi::persistence_matrix::Set_chain_column;
using Gudhi::persistence_matrix::Unordered_set_chain_column;
using Gudhi::persistence_matrix::Vector_chain_column;
using Gudhi::persistence_matrix::Z2_heap_chain_column;
using Gudhi::persistence_matrix::Z2_list_chain_column;
using Gudhi::persistence_matrix::Z2_set_chain_column;
using Gudhi::persistence_matrix::Z2_unordered_set_chain_column;
using Gudhi::persistence_matrix::Z2_vector_chain_column;
using Gudhi::persistence_matrix::Intrusive_list_chain_column;
using Gudhi::persistence_matrix::Intrusive_set_chain_column;
using Gudhi::persistence_matrix::Z2_intrusive_list_chain_column;
using Gudhi::persistence_matrix::Z2_intrusive_set_chain_column;

using Gudhi::persistence_matrix::Base_cell;
using Gudhi::persistence_matrix::Z2_base_cell;
using Gudhi::persistence_matrix::Row_cell;
using Gudhi::persistence_matrix::Z2_row_cell;
using Gudhi::persistence_matrix::Intrusive_row_cell;
using Gudhi::persistence_matrix::Z2_intrusive_row_cell;
using Gudhi::persistence_matrix::Intrusive_list_cell;
using Gudhi::persistence_matrix::Intrusive_set_cell;
using Gudhi::persistence_matrix::Intrusive_list_row_cell;
using Gudhi::persistence_matrix::Intrusive_set_row_cell;
using Gudhi::persistence_matrix::Z2_intrusive_list_cell;
using Gudhi::persistence_matrix::Z2_intrusive_set_cell;
using Gudhi::persistence_matrix::Z2_intrusive_list_row_cell;
using Gudhi::persistence_matrix::Z2_intrusive_set_row_cell;
using Gudhi::persistence_matrix::base_hook_matrix_row;

using Gudhi::persistence_matrix::Row_access;
using Gudhi::persistence_matrix::dimension_type;

using Z5 = Zp_field_element<5>;
using Z2 = Zp_field_element<2>;
using dict_type = std::vector<unsigned int>;

template<class Cell_type>
struct RowCellComp : std::binary_function<const Cell_type&, const Cell_type&, bool> {
	bool operator()(const Cell_type& c1, const Cell_type& c2) const
	{
		if (c1.get_row_index() == c2.get_row_index())
			return c1.get_column_index() < c2.get_column_index();
		return c1.get_row_index() < c2.get_row_index();
	}
};

template<class Cell_type>
using Row_type = boost::intrusive::list <
					Cell_type
				  , boost::intrusive::constant_time_size<false>
				  , boost::intrusive::base_hook< base_hook_matrix_row >
				>;
template<class Cell_type>
using row_container_type = std::vector<Row_type<Cell_type> >;
template<class Cell_type>
using row_container_set = std::vector<std::set<Cell_type,RowCellComp<Cell_type> > >;
template<class Cell_type>
using RA = Row_access<row_container_type<Cell_type>,Cell_type,true,false>;
template<class Cell_type>
using RA_set = Row_access<row_container_set<Cell_type>,Cell_type,false,false>;

struct Dummy_column_pairing{
	Dummy_column_pairing& operator=([[maybe_unused]] Dummy_column_pairing other){return *this;}
	friend void swap([[maybe_unused]] Dummy_column_pairing& d1, [[maybe_unused]] Dummy_column_pairing& d2){}

	Dummy_column_pairing(){}
	Dummy_column_pairing([[maybe_unused]] const Dummy_column_pairing& toCopy){}
	Dummy_column_pairing([[maybe_unused]] Dummy_column_pairing&& other) noexcept{}

	static constexpr bool isActive_ = false;
};

struct Dummy_row_access{
	friend void swap([[maybe_unused]] Dummy_row_access& d1, [[maybe_unused]] Dummy_row_access& d2){}

	Dummy_row_access(){}
	Dummy_row_access([[maybe_unused]] Dummy_row_access&& other) noexcept{}

	static constexpr bool isActive_ = false;
};

////////////// BASICS //////////////

typedef boost::mpl::list<List_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Set_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Unordered_set_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Vector_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Intrusive_list_column<Z5, Intrusive_list_cell<Z5>, Dummy_row_access>,
							Intrusive_set_column<Z5, Intrusive_set_cell<Z5>, Dummy_row_access>,
							List_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Set_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Unordered_set_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Vector_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Intrusive_list_boundary_column<Z5, Intrusive_list_cell<Z5>, Dummy_row_access>,
							Intrusive_set_boundary_column<Z5, Intrusive_set_cell<Z5>, Dummy_row_access>
						> list_of_5_columns;

typedef boost::mpl::list<Z2_heap_column,
							Z2_list_column<Z2_base_cell, Dummy_row_access>,
							Z2_set_column<Z2_base_cell, Dummy_row_access>,
							Z2_unordered_set_column<Z2_base_cell, Dummy_row_access>,
							Z2_vector_column<Z2_base_cell, Dummy_row_access>,
							Z2_intrusive_list_column<Z2_intrusive_list_cell, Dummy_row_access>,
							Z2_intrusive_set_column<Z2_intrusive_set_cell, Dummy_row_access>,
							Z2_heap_boundary_column,
							Z2_list_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_set_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_unordered_set_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_vector_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_intrusive_list_boundary_column<Z2_intrusive_list_cell, Dummy_row_access>,
							Z2_intrusive_set_boundary_column<Z2_intrusive_set_cell, Dummy_row_access>
						> list_of_2_columns;

typedef boost::mpl::list<List_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Set_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Unordered_set_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Vector_boundary_column<Z5, Base_cell<Z5>, Dummy_row_access>,
							Intrusive_list_boundary_column<Z5, Intrusive_list_cell<Z5>, Dummy_row_access>,
							Intrusive_set_boundary_column<Z5, Intrusive_set_cell<Z5>, Dummy_row_access>
						> list_of_5_boundary_columns;

typedef boost::mpl::list<Z2_heap_boundary_column,
							Z2_list_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_set_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_unordered_set_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_vector_boundary_column<Z2_base_cell, Dummy_row_access>,
							Z2_intrusive_list_boundary_column<Z2_intrusive_list_cell, Dummy_row_access>,
							Z2_intrusive_set_boundary_column<Z2_intrusive_set_cell, Dummy_row_access>
						> list_of_2_boundary_columns;

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Set_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Unordered_set_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Vector_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_cell<Z5>, Dummy_row_access>,
							Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_cell<Z5>, Dummy_row_access>
						> list_of_5_chain_columns;

typedef boost::mpl::list<Z2_heap_chain_column<dict_type>,
							Z2_list_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_set_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_unordered_set_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_vector_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_cell, Dummy_row_access>,
							Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_cell, Dummy_row_access>
						> list_of_2_chain_columns;

////////////// ROWS //////////////

typedef boost::mpl::list<List_column<Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >,
							Vector_column<Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >,
							List_boundary_column<Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >,
							Vector_boundary_column<Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >
						> list_of_5_columns_with_row_non_intr;

typedef boost::mpl::list<List_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Set_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Unordered_set_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Vector_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							List_boundary_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Set_boundary_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Unordered_set_boundary_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Vector_boundary_column<Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >
						> list_of_5_columns_with_row_set_non_intr;

typedef boost::mpl::list<Intrusive_list_column<Z5, Intrusive_list_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > >,
							Intrusive_list_boundary_column<Z5, Intrusive_list_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_5_columns_with_row_intr_list;

typedef boost::mpl::list<Intrusive_list_column<Z5, Intrusive_list_row_cell<Row_cell<Z5> >, RA_set<Intrusive_list_row_cell<Row_cell<Z5> > > >,
							Intrusive_list_boundary_column<Z5, Intrusive_list_row_cell<Row_cell<Z5> >, RA_set<Intrusive_list_row_cell<Row_cell<Z5> > > >
						> list_of_5_columns_with_row_set_intr_list;

typedef boost::mpl::list<Intrusive_set_column<Z5, Intrusive_set_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > >,
							Intrusive_set_boundary_column<Z5, Intrusive_set_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_5_columns_with_row_intr_set;

typedef boost::mpl::list<Intrusive_set_column<Z5, Intrusive_set_row_cell<Row_cell<Z5> >, RA_set<Intrusive_set_row_cell<Row_cell<Z5> > > >,
							Intrusive_set_boundary_column<Z5, Intrusive_set_row_cell<Row_cell<Z5> >, RA_set<Intrusive_set_row_cell<Row_cell<Z5> > > >
						> list_of_5_columns_with_row_set_intr_set;

typedef boost::mpl::list<Z2_list_column<Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >,
							Z2_vector_column<Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >,
							Z2_list_boundary_column<Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >,
							Z2_vector_boundary_column<Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >
						> list_of_2_columns_with_row_non_intr;

typedef boost::mpl::list<Z2_list_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_set_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_unordered_set_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_vector_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_list_boundary_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_set_boundary_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_unordered_set_boundary_column<Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_vector_boundary_column<Z2_row_cell, RA_set<Z2_row_cell> >
						> list_of_2_columns_with_row_set_non_intr;

typedef boost::mpl::list<Z2_intrusive_list_column<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > >,
							Z2_intrusive_list_boundary_column<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > >
						> list_of_2_columns_with_row_intr_list;

typedef boost::mpl::list<Z2_intrusive_list_column<Z2_intrusive_list_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_list_row_cell<Z2_row_cell> > >,
							Z2_intrusive_list_boundary_column<Z2_intrusive_list_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_list_row_cell<Z2_row_cell> > >
						> list_of_2_columns_with_row_set_intr_list;

typedef boost::mpl::list<Z2_intrusive_set_column<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > >,
							Z2_intrusive_set_boundary_column<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > >
						> list_of_2_columns_with_row_intr_set;

typedef boost::mpl::list<Z2_intrusive_set_column<Z2_intrusive_set_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_set_row_cell<Z2_row_cell> > >,
							Z2_intrusive_set_boundary_column<Z2_intrusive_set_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_set_row_cell<Z2_row_cell> > >
						> list_of_2_columns_with_row_set_intr_set;

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >,
							Vector_chain_column<dict_type, Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >
						> list_of_5_chain_columns_with_row_non_intr;

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Set_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Unordered_set_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Vector_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >
						> list_of_5_chain_columns_with_row_set_non_intr;

typedef boost::mpl::list<Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_5_chain_columns_with_row_intr_list;

typedef boost::mpl::list<Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_row_cell<Row_cell<Z5> >, RA_set<Intrusive_list_row_cell<Row_cell<Z5> > > >
						> list_of_5_chain_columns_with_row_set_intr_list;

typedef boost::mpl::list<Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_5_chain_columns_with_row_intr_set;

typedef boost::mpl::list<Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_row_cell<Row_cell<Z5> >, RA_set<Intrusive_set_row_cell<Row_cell<Z5> > > >
						> list_of_5_chain_columns_with_row_set_intr_set;

typedef boost::mpl::list<Z2_list_chain_column<dict_type, Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >,
							Z2_vector_chain_column<dict_type, Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >
						> list_of_2_chain_columns_with_row_non_intr;

typedef boost::mpl::list<Z2_list_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_set_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_unordered_set_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_vector_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >
						> list_of_2_chain_columns_with_row_set_non_intr;

typedef boost::mpl::list<Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > >
						> list_of_2_chain_columns_with_row_intr_list;

typedef boost::mpl::list<Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_list_row_cell<Z2_row_cell> > >
						> list_of_2_chain_columns_with_row_set_intr_list;

typedef boost::mpl::list<Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > >
						> list_of_2_chain_columns_with_row_intr_set;

typedef boost::mpl::list<Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_set_row_cell<Z2_row_cell> > >
						> list_of_2_chain_columns_with_row_set_intr_set;

////////////// PAIRING //////////////

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Set_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Unordered_set_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Vector_chain_column<dict_type, Z5, Base_cell<Z5>, Dummy_row_access>,
							Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_cell<Z5>, Dummy_row_access>,
							Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_cell<Z5>, Dummy_row_access>,
							Z2_heap_chain_column<dict_type>,
							Z2_list_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_set_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_unordered_set_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_vector_chain_column<dict_type, Z2_base_cell, Dummy_row_access>,
							Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_cell, Dummy_row_access>,
							Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_cell, Dummy_row_access>
						> list_of_chain_pairing_columns;

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >,
							Vector_chain_column<dict_type, Z5, Intrusive_row_cell<Z5>, RA<Intrusive_row_cell<Z5> > >
						> list_of_pairing_chain_columns_with_row_non_intr;

typedef boost::mpl::list<List_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Set_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Unordered_set_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >,
							Vector_chain_column<dict_type, Z5, Row_cell<Z5>, RA_set<Row_cell<Z5> > >
						> list_of_pairing_chain_columns_with_set_row_non_intr;

typedef boost::mpl::list<Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_pairing_chain_columns_with_row_intr_list;

typedef boost::mpl::list<Intrusive_list_chain_column<dict_type, Z5, Intrusive_list_row_cell<Row_cell<Z5> >, RA_set<Intrusive_list_row_cell<Row_cell<Z5> > > >
						> list_of_pairing_chain_columns_with_set_row_intr_list;

typedef boost::mpl::list<Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_row_cell<Intrusive_row_cell<Z5> >, RA<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > >
						> list_of_pairing_chain_columns_with_row_intr_set;

typedef boost::mpl::list<Intrusive_set_chain_column<dict_type, Z5, Intrusive_set_row_cell<Row_cell<Z5> >, RA_set<Intrusive_set_row_cell<Row_cell<Z5> > > >
						> list_of_pairing_chain_columns_with_set_row_intr_set;

typedef boost::mpl::list<Z2_list_chain_column<dict_type, Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >,
							Z2_vector_chain_column<dict_type, Z2_intrusive_row_cell, RA<Z2_intrusive_row_cell> >
						> list_of_pairing_chain_columns_with_row_z2_non_intr;

typedef boost::mpl::list<Z2_list_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_set_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_unordered_set_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >,
							Z2_vector_chain_column<dict_type, Z2_row_cell, RA_set<Z2_row_cell> >
						> list_of_pairing_chain_columns_with_set_row_z2_non_intr;

typedef boost::mpl::list<Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > >
						> list_of_pairing_chain_columns_with_row_z2_intr_list;

typedef boost::mpl::list<Z2_intrusive_list_chain_column<dict_type, Z2_intrusive_list_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_list_row_cell<Z2_row_cell> > >
						> list_of_pairing_chain_columns_with_set_row_z2_intr_list;

typedef boost::mpl::list<Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_row_cell<Z2_intrusive_row_cell>, RA<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > >
						> list_of_pairing_chain_columns_with_row_z2_intr_set;

typedef boost::mpl::list<Z2_intrusive_set_chain_column<dict_type, Z2_intrusive_set_row_cell<Z2_row_cell>, RA_set<Z2_intrusive_set_row_cell<Z2_row_cell> > >
						> list_of_pairing_chain_columns_with_set_row_z2_intr_set;

template<class Column>
void common_5_test(std::vector<Column> &matrix)
{
	std::set<std::pair<unsigned int,unsigned int> > rows;
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	auto it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[1]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[2]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 0);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 5);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(!matrix[1].is_empty());
	BOOST_CHECK(!matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 3);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].begin() == matrix[1].end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(!matrix[1].is_non_zero(2));

	matrix[0] *= 4;
	matrix[1] *= 2;
	matrix[2] *= 0;

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert({c.get_row_index(), c.get_element()});
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(it->first, 1);
	BOOST_CHECK_EQUAL(it->second, 1);
	++it;
	BOOST_CHECK_EQUAL(it->first, 2);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK_EQUAL(it->first, 3);
	BOOST_CHECK_EQUAL(it->second, 2);
	++it;
	BOOST_CHECK_EQUAL(it->first, 6);
	BOOST_CHECK_EQUAL(it->second, 4);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].begin() == matrix[1].end());
	BOOST_CHECK(matrix[2].begin() == matrix[2].end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(!matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	swap(matrix[3], matrix[4]);
	auto res = matrix[4].get_content(7);
	for (auto& f : res) BOOST_CHECK_EQUAL(f, 0u);
	res = matrix[3].get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 4u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	Column move(std::move(matrix[5]));
	res = move.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 1u);
	BOOST_CHECK_EQUAL(res[6], 1u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_5, Column, list_of_5_columns) {
	std::vector<Column> matrix;

	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4));
	matrix.push_back(Column());
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}));
	matrix.push_back(Column(matrix[1]));

	Column add = matrix[4] + matrix[5];
	auto res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0u);
	BOOST_CHECK_EQUAL(res[1], 4u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 0u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	Column mul1 = matrix[4] * 2;
	res = mul1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 2u);
	BOOST_CHECK_EQUAL(res[1], 4u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 1u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 3u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	Column mul2 = 2 * matrix[4];
	BOOST_CHECK(res == mul2.get_content(7));

	common_5_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_common_5, Column, list_of_5_chain_columns) {
	std::vector<Column> matrix;
	dict_type dict(7);

	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, dict));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4, dict));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4, dict));
	matrix.push_back(Column(dict));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, dict));
	matrix.push_back(Column(matrix[1]));

	Column add = matrix[4] + matrix[5];
	auto res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0u);
	BOOST_CHECK_EQUAL(res[1], 4u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 0u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	Column mul1 = matrix[4] * 2;
	res = mul1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 2u);
	BOOST_CHECK_EQUAL(res[1], 4u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 1u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 3u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	Column mul2 = 2 * matrix[4];
	BOOST_CHECK(res == mul2.get_content(7));

	common_5_test(matrix);
}

template<class Column, class Rows>
void common_5_test_with_rows(Rows &rows){
	std::vector<Column> matrix;

	matrix.push_back(Column(0, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, rows));
	matrix.push_back(Column(1, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4, rows));
	matrix.push_back(Column(2, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4, rows));
	matrix.push_back(Column(3, rows));
	matrix.push_back(Column(4, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, rows));
	matrix.push_back(Column(matrix[1], 5u));

	common_5_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_5_non_intr, Column, list_of_5_columns_with_row_non_intr) {
	row_container_type<Intrusive_row_cell<Z5> > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_5_non_intr, Column, list_of_5_columns_with_row_set_non_intr) {
	row_container_set<Row_cell<Z5> > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_5_intr_list, Column, list_of_5_columns_with_row_intr_list) {
	row_container_type<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_5_intr_list, Column, list_of_5_columns_with_row_set_intr_list) {
	row_container_set<Intrusive_list_row_cell<Row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_5_intr_set, Column, list_of_5_columns_with_row_intr_set) {
	row_container_type<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_5_intr_set, Column, list_of_5_columns_with_row_set_intr_set) {
	row_container_set<Intrusive_set_row_cell<Row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_test_with_rows<Column>(rows);
}

template<class Column, class Rows>
void common_5_chain_test_with_rows(Rows &rows){
	std::vector<Column> matrix;
	dict_type dict(7);

	matrix.push_back(Column(0, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, rows, dict));
	matrix.push_back(Column(1, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4, rows, dict));
	matrix.push_back(Column(2, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(3)},{2,Z5(4)},{5,Z5(4)},{6,Z5(4)}}, 4, rows, dict));
	matrix.push_back(Column(3, rows, dict));
	matrix.push_back(Column(4, std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}}, 4, rows, dict));
	matrix.push_back(Column(matrix[1], 5u));

	common_5_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_5_non_intr, Column, list_of_5_chain_columns_with_row_non_intr) {
	row_container_type<Intrusive_row_cell<Z5> > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_5_non_intr, Column, list_of_5_chain_columns_with_row_set_non_intr) {
	row_container_set<Row_cell<Z5> > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_5_intr_list, Column, list_of_5_chain_columns_with_row_intr_list) {
	row_container_type<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_5_intr_list, Column, list_of_5_chain_columns_with_row_set_intr_list) {
	row_container_set<Intrusive_list_row_cell<Row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_5_intr_set, Column, list_of_5_chain_columns_with_row_intr_set) {
	row_container_type<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_5_intr_set, Column, list_of_5_chain_columns_with_row_set_intr_set) {
	row_container_set<Intrusive_set_row_cell<Row_cell<Z5> > > rows;	//defined here to ensure it is destroyed after the columns
	common_5_chain_test_with_rows<Column>(rows);
}

template<class Column>
void common_2_test(std::vector<Column> &matrix)
{
	std::set<unsigned int> rows;
	for (auto c : matrix[0]){
		rows.insert(c.get_row_index());
	}
	auto it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 3);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[1]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	rows.clear();
	for (auto c : matrix[2]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 0);
	++it;
	BOOST_CHECK_EQUAL(*it, 1);
	++it;
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 5);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(!matrix[1].is_empty());
	BOOST_CHECK(!matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(1));
	BOOST_CHECK(matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	matrix[0] += matrix[1];
	matrix[1] += matrix[2];

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 3);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[1].is_empty());		//forces heap column to prune.
	BOOST_CHECK(matrix[1].begin() == matrix[1].end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(2));
	BOOST_CHECK(!matrix[1].is_non_zero(2));

	matrix[0] *= 5;
	matrix[2] *= 8;

	rows.clear();
	for (auto c : matrix[0]){
		rows.insert(c.get_row_index());
	}
	it = rows.begin();
	BOOST_CHECK_EQUAL(*it, 2);
	++it;
	BOOST_CHECK_EQUAL(*it, 3);
	++it;
	BOOST_CHECK_EQUAL(*it, 6);
	++it;
	BOOST_CHECK(it == rows.end());

	BOOST_CHECK(matrix[2].begin() == matrix[2].end());

	BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
	BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
	BOOST_CHECK(!matrix[0].is_empty());
	BOOST_CHECK(matrix[1].is_empty());
	BOOST_CHECK(matrix[2].is_empty());
	BOOST_CHECK(matrix[0].is_non_zero(2));
	BOOST_CHECK(!matrix[1].is_non_zero(2));
	BOOST_CHECK(!matrix[2].is_non_zero(3));

	swap(matrix[3], matrix[4]);
	auto res = matrix[4].get_content(7);
	for (auto f : res) BOOST_CHECK_EQUAL(f, 0);
	res = matrix[3].get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	Column move(std::move(matrix[5]));
	res = move.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_common_2, Column, list_of_2_columns) {
	std::vector<Column> matrix;

	matrix.push_back(Column(std::vector<unsigned int>{0,1,3,5}, 4));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4));
	matrix.push_back(Column());
	matrix.push_back(Column(std::vector<unsigned int>{0,1,3,5}));
	matrix.push_back(Column(matrix[1]));

	Column add = matrix[4] + matrix[5];
	auto res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0);
	BOOST_CHECK_EQUAL(res[1], 0);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 0);
	BOOST_CHECK_EQUAL(res[6], 1);

	Column mul1 = matrix[4] * 3;
	res = mul1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	Column mul2 = 4 * matrix[4];
	BOOST_CHECK(mul2.is_empty());
	for (auto f : mul2.get_content(7))
		BOOST_CHECK_EQUAL(f, 0);

	common_2_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_common_2, Column, list_of_2_chain_columns) {
	std::vector<Column> matrix;
	dict_type dict(7);

	matrix.push_back(Column(std::vector<unsigned int>{0,1,3,5}, 4, dict));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4, dict));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,2,5,6}, 4, dict));
	matrix.push_back(Column(dict));
	matrix.push_back(Column(std::vector<unsigned int>{0,1,3,5}, 4, dict));
	matrix.push_back(Column(matrix[1]));

	Column add = matrix[4] + matrix[5];
	auto res = add.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 0);
	BOOST_CHECK_EQUAL(res[1], 0);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 0);
	BOOST_CHECK_EQUAL(res[6], 1);

	Column mul1 = matrix[4] * 3;
	res = mul1.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	Column mul2 = 4 * matrix[4];
	BOOST_CHECK(mul2.is_empty());
	for (auto f : mul2.get_content(7))
		BOOST_CHECK_EQUAL(f, 0);

	common_2_test(matrix);
}

template<class Column, class Rows>
void common_2_test_with_rows(Rows &rows){
	std::vector<Column> matrix;

	matrix.push_back(Column(0, std::vector<unsigned int>{0,1,3,5}, 4, rows));
	matrix.push_back(Column(1, std::vector<unsigned int>{0,1,2,5,6}, 4, rows));
	matrix.push_back(Column(2, std::vector<unsigned int>{0,1,2,5,6}, 4, rows));
	matrix.push_back(Column(3, rows));
	matrix.push_back(Column(4, std::vector<unsigned int>{0,1,3,5}, rows));
	matrix.push_back(Column(matrix[1], 5u));

	common_2_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_2_non_intr, Column, list_of_2_columns_with_row_non_intr) {
	row_container_type<Z2_intrusive_row_cell> rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_2_non_intr, Column, list_of_2_columns_with_row_set_non_intr) {
	row_container_set<Z2_row_cell> rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_2_intr_list, Column, list_of_2_columns_with_row_intr_list) {
	row_container_type<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_2_intr_list, Column, list_of_2_columns_with_row_set_intr_list) {
	row_container_set<Z2_intrusive_list_row_cell<Z2_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_common_2_intr_set, Column, list_of_2_columns_with_row_intr_set) {
	row_container_type<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_set_row_common_2_intr_set, Column, list_of_2_columns_with_row_set_intr_set) {
	row_container_set<Z2_intrusive_set_row_cell<Z2_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_test_with_rows<Column>(rows);
}

template<class Column, class Rows>
void common_2_chain_test_with_rows(Rows &rows){
	std::vector<Column> matrix;
	dict_type dict(7);

	matrix.push_back(Column(0, std::vector<unsigned int>{0,1,3,5}, 4, rows, dict));
	matrix.push_back(Column(1, std::vector<unsigned int>{0,1,2,5,6}, 4, rows, dict));
	matrix.push_back(Column(2, std::vector<unsigned int>{0,1,2,5,6}, 4, rows, dict));
	matrix.push_back(Column(3, rows, dict));
	matrix.push_back(Column(4, std::vector<unsigned int>{0,1,3,5}, 4, rows, dict));
	matrix.push_back(Column(matrix[1], 5u));

	common_2_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_2_non_intr, Column, list_of_2_chain_columns_with_row_non_intr) {
	row_container_type<Z2_intrusive_row_cell> rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_2_non_intr, Column, list_of_2_chain_columns_with_row_set_non_intr) {
	row_container_set<Z2_row_cell> rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_2_intr_list, Column, list_of_2_chain_columns_with_row_intr_list) {
	row_container_type<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_2_intr_list, Column, list_of_2_chain_columns_with_row_set_intr_list) {
	row_container_set<Z2_intrusive_list_row_cell<Z2_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_common_2_intr_set, Column, list_of_2_chain_columns_with_row_intr_set) {
	row_container_type<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_common_2_intr_set, Column, list_of_2_chain_columns_with_row_set_intr_set) {
	row_container_set<Z2_intrusive_set_row_cell<Z2_row_cell> > rows;	//defined here to ensure it is destroyed after the columns
	common_2_chain_test_with_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Boundary_column_types_methods, Column, list_of_5_boundary_columns) {
	Column col1;
	Column col2(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)},{1,Z5(2)},{3,Z5(3)},{5,Z5(4)}});
	Column col3(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)},{1,Z5(2)},{2,Z5(1)},{5,Z5(1)},{6,Z5(1)}}, 4);
	Column col4(col3);

	BOOST_CHECK_EQUAL(col2.get_dimension(), 3);
	BOOST_CHECK_EQUAL(col3.get_dimension(), 4);

	auto res = col1.get_content(7);
	for (auto& f : res) BOOST_CHECK_EQUAL(f, 0u);

	res = col2.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 0u);
	BOOST_CHECK_EQUAL(res[3], 3u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 4u);
	BOOST_CHECK_EQUAL(res[6], 0u);

	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 1u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	BOOST_CHECK(res == col4.get_content(7));

	BOOST_CHECK_EQUAL(col1.get_pivot(), -1);
	BOOST_CHECK_EQUAL(col2.get_pivot(), 5);
	BOOST_CHECK_EQUAL(col3.get_pivot(), 6);
	BOOST_CHECK_EQUAL(col4.get_pivot(), 6);
	BOOST_CHECK_EQUAL(col1.get_pivot_value(), 0u);
	BOOST_CHECK_EQUAL(col2.get_pivot_value(), 4u);
	BOOST_CHECK_EQUAL(col3.get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(col4.get_pivot_value(), 1u);

	col4.clear(5);
	res = col4.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 2u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 0u);
	BOOST_CHECK_EQUAL(res[5], 0u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	BOOST_CHECK_EQUAL(col4.get_pivot(), 6);
	BOOST_CHECK_EQUAL(col4.get_pivot_value(), 1u);

	col4.clear();
	BOOST_CHECK(col4.is_empty());
	BOOST_CHECK(col1.get_content(7) == col4.get_content(7));

	BOOST_CHECK_EQUAL(col4.get_pivot(), -1);
	BOOST_CHECK_EQUAL(col4.get_pivot_value(), 0u);

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col3.reorder(permutation);
	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 4u);
	BOOST_CHECK_EQUAL(res[1], 0u);
	BOOST_CHECK_EQUAL(res[2], 1u);
	BOOST_CHECK_EQUAL(res[3], 0u);
	BOOST_CHECK_EQUAL(res[4], 1u);
	BOOST_CHECK_EQUAL(res[5], 2u);
	BOOST_CHECK_EQUAL(res[6], 1u);

	BOOST_CHECK_EQUAL(col3.get_pivot(), 6);
	BOOST_CHECK_EQUAL(col3.get_pivot_value(), 1u);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_boundary_column_types_methods, Column, list_of_2_boundary_columns) {
	Column col1;
	Column col2(std::vector<unsigned int>{0,1,3,5});
	Column col3(std::vector<unsigned int>{0,1,2,5,6}, 4);
	Column col4(col3);

	BOOST_CHECK_EQUAL(col2.get_dimension(), 3);
	BOOST_CHECK_EQUAL(col3.get_dimension(), 4);

	auto res = col1.get_content(7);
	for (auto f : res) BOOST_CHECK_EQUAL(f, 0);

	res = col2.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 0);
	BOOST_CHECK_EQUAL(res[3], 1);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 0);

	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);

	BOOST_CHECK(res == col4.get_content(7));

	BOOST_CHECK_EQUAL(col1.get_pivot(), -1);
	BOOST_CHECK_EQUAL(col2.get_pivot(), 5);
	BOOST_CHECK_EQUAL(col3.get_pivot(), 6);
	BOOST_CHECK_EQUAL(col4.get_pivot(), 6);

	col4.clear(5);
	res = col4.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 1);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 0);
	BOOST_CHECK_EQUAL(res[5], 0);
	BOOST_CHECK_EQUAL(res[6], 1);

	BOOST_CHECK_EQUAL(col4.get_pivot(), 6);

	col4.clear();
	BOOST_CHECK(col4.is_empty());
	BOOST_CHECK(col1.get_content(7) == col4.get_content(7));

	BOOST_CHECK_EQUAL(col4.get_pivot(), -1);

	std::vector<unsigned int> permutation{0,5,2,3,1,4,6};
	col3.reorder(permutation);
	res = col3.get_content(7);
	BOOST_CHECK_EQUAL(res[0], 1);
	BOOST_CHECK_EQUAL(res[1], 0);
	BOOST_CHECK_EQUAL(res[2], 1);
	BOOST_CHECK_EQUAL(res[3], 0);
	BOOST_CHECK_EQUAL(res[4], 1);
	BOOST_CHECK_EQUAL(res[5], 1);
	BOOST_CHECK_EQUAL(res[6], 1);

	BOOST_CHECK_EQUAL(col3.get_pivot(), 6);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_methods, Column, list_of_5_chain_columns) {
	std::vector<Column> matrix;
	dict_type pivotToColumnIndex{0,1,2,3};

	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(1)}}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(4)}, {1,Z5(3)}}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{1,Z5(2)},{2,Z5(4)}}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<std::pair<unsigned int,Z5> >{{0,Z5(2)},{2,Z5(3)},{3,Z5(1)}}, 4, pivotToColumnIndex));
	matrix.push_back(Column(pivotToColumnIndex));
	matrix.push_back(matrix[0]);

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 0);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 1);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 2);
	BOOST_CHECK_EQUAL(matrix[3].get_pivot(), 3);
	BOOST_CHECK_EQUAL(matrix[4].get_pivot(), -1);
	BOOST_CHECK_EQUAL(matrix[5].get_pivot(), 0);
	BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 3u);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot_value(), 4u);
	BOOST_CHECK_EQUAL(matrix[3].get_pivot_value(), 1u);
	BOOST_CHECK_EQUAL(matrix[4].get_pivot_value(), 0u);
	BOOST_CHECK_EQUAL(matrix[5].get_pivot_value(), 1u);

	BOOST_CHECK_EQUAL(pivotToColumnIndex[0], 0);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[1], 1);

	matrix[0] += matrix[1];

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 0);
	BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 3u);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 4u);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[0], 1);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[1], 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_methods, Column, list_of_2_chain_columns) {
	std::vector<Column> matrix;
	dict_type pivotToColumnIndex{0,1,2,3};

	matrix.push_back(Column(std::vector<unsigned int>{0}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<unsigned int>{0,1}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<unsigned int>{1,2}, 4, pivotToColumnIndex));
	matrix.push_back(Column(std::vector<unsigned int>{0,2,3}, 4, pivotToColumnIndex));
	matrix.push_back(Column(pivotToColumnIndex));
	matrix.push_back(matrix[0]);

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 0);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 1);
	BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 2);
	BOOST_CHECK_EQUAL(matrix[3].get_pivot(), 3);
	BOOST_CHECK_EQUAL(matrix[4].get_pivot(), -1);
	BOOST_CHECK_EQUAL(matrix[5].get_pivot(), 0);

	BOOST_CHECK_EQUAL(pivotToColumnIndex[0], 0);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[1], 1);

	matrix[0] += matrix[1];

	BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 0);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[0], 1);
	BOOST_CHECK_EQUAL(pivotToColumnIndex[1], 0);
}

template<class Column, class Row_container_type>
void row_methods_test_com(std::vector<Column>& matrix,
						  std::vector<std::vector<std::pair<unsigned int,Z5> > >& chains,
						  Row_container_type& row_cont)
{
	std::vector<std::vector<std::pair<unsigned int,Z5> > > rows;
	rows.push_back({{0,Z5(1)},{3,Z5(2)},{5,Z5(1)}});
	rows.push_back({{1,Z5(3)},{2,Z5(2)}});
	rows.push_back({{2,Z5(4)},{3,Z5(3)}});
	rows.push_back({});
	rows.push_back({{3,Z5(1)}});

	//unordered sets are unordered...
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_columns(matrix.size());
	for (unsigned int i = 0; i < matrix.size(); ++i){
		for (auto& cell : matrix[i]){
			ordered_columns[i].insert({cell.get_row_index(), cell.get_element()});
		}
	}
	for (unsigned int i = 0; i < matrix.size(); ++i){
		auto itCol = ordered_columns[i].begin();
		auto itChain = chains[i].begin();
		while (itChain != chains[i].end()){
			BOOST_CHECK_EQUAL(itCol->first, itChain->first);
			BOOST_CHECK_EQUAL(itCol->second, itChain->second);
			++itCol; ++itChain;
		}
	}

	BOOST_CHECK_EQUAL(row_cont.size(), rows.size());
	//rows are unordered
	std::vector<std::set<std::pair<unsigned int,Z5> > > ordered_rows(row_cont.size());
	for (unsigned int i = 0; i < row_cont.size(); ++i){
		for (auto& cell : row_cont[i]){
			ordered_rows[i].insert({cell.get_column_index(), cell.get_element()});
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(itRow->first, itChain->first);
			BOOST_CHECK_EQUAL(itRow->second, itChain->second);
			++itRow; ++itChain;
		}
	}
}

template<class Column, class Row_container_type>
void row_methods_test(Row_container_type& row_cont) {
	std::vector<Column> matrix;

	std::vector<std::vector<std::pair<unsigned int,Z5> > > chains;
	chains.push_back({{0,Z5(1)}});
	chains.push_back({{1,Z5(3)}});
	chains.push_back({{1,Z5(2)},{2,Z5(4)}});
	chains.push_back({{0,Z5(2)},{2,Z5(3)},{4,Z5(1)}});
	chains.emplace_back();
	chains.push_back({{0,Z5(1)}});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, row_cont));
	}
	matrix.push_back(Column(4, row_cont));
	matrix.push_back(Column(matrix[0], 5u));

	row_methods_test_com(matrix, chains, row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_methods_non_intr, Column, list_of_5_columns_with_row_non_intr) {
	row_container_type<Intrusive_row_cell<Z5> > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_set_methods_non_intr, Column, list_of_5_columns_with_row_set_non_intr) {
	row_container_set<Row_cell<Z5> > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_methods_intr_list, Column, list_of_5_columns_with_row_intr_list) {
	row_container_type<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_set_methods_intr_list, Column, list_of_5_columns_with_row_set_intr_list) {
	row_container_set<Intrusive_list_row_cell<Row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_methods_intr_set, Column, list_of_5_columns_with_row_intr_set) {
	row_container_type<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Column_types_with_row_set_methods_intr_set, Column, list_of_5_columns_with_row_set_intr_set) {
	row_container_set<Intrusive_set_row_cell<Row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	row_methods_test<Column>(row_cont);
}

template<class Column, class Row_container_type>
void chain_row_methods_test(Row_container_type& row_cont) {
	std::vector<Column> matrix;
	dict_type dict(7);

	std::vector<std::vector<std::pair<unsigned int,Z5> > > chains;
	chains.push_back({{0,Z5(1)}});
	chains.push_back({{1,Z5(3)}});
	chains.push_back({{1,Z5(2)},{2,Z5(4)}});
	chains.push_back({{0,Z5(2)},{2,Z5(3)},{4,Z5(1)}});
	chains.emplace_back();
	chains.push_back({{0,Z5(1)}});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, row_cont, dict));
	}
	matrix.push_back(Column(4, row_cont, dict));
	matrix.push_back(Column(matrix[0], 5u));

	row_methods_test_com(matrix, chains, row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_methods_non_intr, Column, list_of_5_chain_columns_with_row_non_intr) {
	row_container_type<Intrusive_row_cell<Z5> > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_set_methods_non_intr, Column, list_of_5_chain_columns_with_row_set_non_intr) {
	row_container_set<Row_cell<Z5> > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_methods_intr_list, Column, list_of_5_chain_columns_with_row_intr_list) {
	row_container_type<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_set_methods_intr_list, Column, list_of_5_chain_columns_with_row_set_intr_list) {
	row_container_set<Intrusive_list_row_cell<Row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_methods_intr_set, Column, list_of_5_chain_columns_with_row_intr_set) {
	row_container_type<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_set_methods_intr_set, Column, list_of_5_chain_columns_with_row_set_intr_set) {
	row_container_set<Intrusive_set_row_cell<Row_cell<Z5> > > row_cont;	//defined here to ensure it is destroyed after the columns
	chain_row_methods_test<Column>(row_cont);
}

template<class Column, class Row_container_type>
void z2_row_methods_test_com(std::vector<Column>& matrix,
							 std::vector<std::vector<unsigned int> >& chains,
							 Row_container_type& row_cont)
{
	std::vector<std::vector<unsigned int> > rows;
	rows.push_back({0,3,5});
	rows.push_back({1,2});
	rows.push_back({2,3});
	rows.push_back({});
	rows.push_back({3});

	//unordered sets and heaps are unordered
	std::vector<std::set<unsigned int> > ordered_columns(matrix.size());
	for (unsigned int i = 0; i < matrix.size(); ++i){
		for (auto& cell : matrix[i]){
			ordered_columns[i].insert(cell.get_row_index());
		}
	}
	for (unsigned int i = 0; i < matrix.size(); ++i){
		auto itCol = ordered_columns[i].begin();
		auto itChain = chains[i].begin();
		while (itChain != chains[i].end()){
			BOOST_CHECK_EQUAL(*itCol, *itChain);
			++itCol; ++itChain;
		}
	}

	BOOST_CHECK_EQUAL(row_cont.size(), rows.size());
	//rows are unordered
	std::vector<std::set<unsigned int> > ordered_rows(row_cont.size());
	for (unsigned int i = 0; i < row_cont.size(); ++i){
		for (auto& cell : row_cont[i]){
			ordered_rows[i].insert(cell.get_column_index());
		}
	}
	for (unsigned int i = 0; i < rows.size(); ++i){
		auto& row = ordered_rows[i];
		auto itRow = row.begin();
		auto itChain = rows[i].begin();
		while (itChain != rows[i].end()){
			BOOST_CHECK_EQUAL(*itRow, *itChain);
			++itRow; ++itChain;
		}
	}
}

template<class Column, class Row_container_type>
void z2_row_methods_test(Row_container_type& row_cont) {
	std::vector<Column> matrix;

	std::vector<std::vector<unsigned int> > chains;
	chains.push_back({0});
	chains.push_back({1});
	chains.push_back({1,2});
	chains.push_back({0,2,4});
	chains.emplace_back();
	chains.push_back({0});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, row_cont));
	}
	matrix.push_back(Column(4, row_cont));
	matrix.push_back(Column(matrix[0], 5u));

	z2_row_methods_test_com(matrix, chains, row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_methods_non_intr, Column, list_of_2_columns_with_row_non_intr) {
	row_container_type<Z2_intrusive_row_cell> row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_set_methods_non_intr, Column, list_of_2_columns_with_row_set_non_intr) {
	row_container_set<Z2_row_cell> row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_methods_intr_list, Column, list_of_2_columns_with_row_intr_list) {
	row_container_type<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_set_methods_intr_list, Column, list_of_2_columns_with_row_set_intr_list) {
	row_container_set<Z2_intrusive_list_row_cell<Z2_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_methods_intr_set, Column, list_of_2_columns_with_row_intr_set) {
	row_container_type<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_column_types_with_row_set_methods_intr_set, Column, list_of_2_columns_with_row_set_intr_set) {
	row_container_set<Z2_intrusive_set_row_cell<Z2_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_row_methods_test<Column>(row_cont);
}

template<class Column, class Row_container_type>
void z2_chain_row_methods_test(Row_container_type& row_cont) {
	std::vector<Column> matrix;
	dict_type dict(7);

	std::vector<std::vector<unsigned int> > chains;
	chains.push_back({0});
	chains.push_back({1});
	chains.push_back({1,2});
	chains.push_back({0,2,4});
	chains.emplace_back();
	chains.push_back({0});

	for (unsigned int i = 0; i < 4; ++i){
		matrix.push_back(Column(i, chains[i], 4, row_cont, dict));
	}
	matrix.push_back(Column(4, row_cont, dict));
	matrix.push_back(Column(matrix[0], 5u));

	z2_row_methods_test_com(matrix, chains, row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_methods_non_intr, Column, list_of_2_chain_columns_with_row_non_intr) {
	row_container_type<Z2_intrusive_row_cell> row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_set_methods_non_intr, Column, list_of_2_chain_columns_with_row_set_non_intr) {
	row_container_set<Z2_row_cell> row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_methods_intr_list, Column, list_of_2_chain_columns_with_row_intr_list) {
	row_container_type<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_set_methods_intr_list, Column, list_of_2_chain_columns_with_row_set_intr_list) {
	row_container_set<Z2_intrusive_list_row_cell<Z2_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_methods_intr_set, Column, list_of_2_chain_columns_with_row_intr_set) {
	row_container_type<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Z2_chain_column_types_with_row_set_methods_intr_set, Column, list_of_2_chain_columns_with_row_set_intr_set) {
	row_container_set<Z2_intrusive_set_row_cell<Z2_row_cell> > row_cont;	//defined here to ensure it is destroyed after the columns
	z2_chain_row_methods_test<Column>(row_cont);
}

template<class Column>
void pairing_test(std::vector<Column> &matrix)
{
	matrix[0].assign_paired_chain(1);
	matrix[1].assign_paired_chain(0);
	matrix[3].assign_paired_chain(2);

	BOOST_CHECK(matrix[0].is_paired());
	BOOST_CHECK(matrix[1].is_paired());
	BOOST_CHECK(!matrix[2].is_paired());
	BOOST_CHECK(matrix[3].is_paired());
	BOOST_CHECK_EQUAL(matrix[0].get_paired_chain_index(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_paired_chain_index(), 0);
	BOOST_CHECK_EQUAL(matrix[2].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[3].get_paired_chain_index(), 2);

	matrix[3].unassign_paired_chain();
	matrix[1].unassign_paired_chain();

	BOOST_CHECK(matrix[0].is_paired());
	BOOST_CHECK(!matrix[1].is_paired());
	BOOST_CHECK(!matrix[2].is_paired());
	BOOST_CHECK(!matrix[3].is_paired());
	BOOST_CHECK_EQUAL(matrix[0].get_paired_chain_index(), 1);
	BOOST_CHECK_EQUAL(matrix[1].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[2].get_paired_chain_index(), -1);
	BOOST_CHECK_EQUAL(matrix[3].get_paired_chain_index(), -1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_pairing_option, Column, list_of_chain_pairing_columns) {
	std::vector<Column> matrix;
	dict_type pivotToColumnIndex{0,0,0,0};
	for (int i = 0; i < 4; ++i)
		matrix.push_back(Column(pivotToColumnIndex));

	pairing_test(matrix);
}

template<class Column, class Rows>
void pairing_test_with_chain_rows(Rows &rows){
	std::vector<Column> matrix;
	dict_type dict(7);

	for (int i = 0; i < 4; ++i)
		matrix.push_back(Column(i, rows, dict));

	pairing_test(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_non_intr, Column, list_of_pairing_chain_columns_with_row_non_intr) {
	row_container_type<Intrusive_row_cell<Z5> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_non_intr, Column, list_of_pairing_chain_columns_with_set_row_non_intr) {
	row_container_set<Row_cell<Z5> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_intr_list, Column, list_of_pairing_chain_columns_with_row_intr_list) {
	row_container_type<Intrusive_list_row_cell<Intrusive_row_cell<Z5> > > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_intr_list, Column, list_of_pairing_chain_columns_with_set_row_intr_list) {
	row_container_set<Intrusive_list_row_cell<Row_cell<Z5> > > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_intr_set, Column, list_of_pairing_chain_columns_with_row_intr_set) {
	row_container_type<Intrusive_set_row_cell<Intrusive_row_cell<Z5> > > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_intr_set, Column, list_of_pairing_chain_columns_with_set_row_intr_set) {
	row_container_set<Intrusive_set_row_cell<Row_cell<Z5> > > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_z2_non_intr, Column, list_of_pairing_chain_columns_with_row_z2_non_intr) {
	row_container_type<Z2_intrusive_row_cell> rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_z2_non_intr, Column, list_of_pairing_chain_columns_with_set_row_z2_non_intr) {
	row_container_set<Z2_row_cell> rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_z2_intr_list, Column, list_of_pairing_chain_columns_with_row_z2_intr_list) {
	row_container_type<Z2_intrusive_list_row_cell<Z2_intrusive_row_cell> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_z2_intr_list, Column, list_of_pairing_chain_columns_with_set_row_z2_intr_list) {
	row_container_set<Z2_intrusive_list_row_cell<Z2_row_cell> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_row_pairing_option_z2_intr_set, Column, list_of_pairing_chain_columns_with_row_z2_intr_set) {
	row_container_type<Z2_intrusive_set_row_cell<Z2_intrusive_row_cell> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Chain_column_types_with_set_row_pairing_option_z2_intr_set, Column, list_of_pairing_chain_columns_with_set_row_z2_intr_set) {
	row_container_set<Z2_intrusive_set_row_cell<Z2_row_cell> > rows;
	pairing_test_with_chain_rows<Column>(rows);
}


