/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_TESTS_MASTER_H
#define PM_COLUMN_TESTS_MASTER_H

#include <map>
#include <set>
#include <vector>

#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>

#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistence_matrix/columns/column_dimension_holder.h>
#include <gudhi/Persistence_matrix/columns/chain_column_extra_properties.h>
#include <gudhi/Persistence_matrix/columns/cell_types.h>
#include <gudhi/Persistence_matrix/columns/row_access.h>
#include <gudhi/Persistence_matrix/allocators/cell_constructors.h>
#include <gudhi/Fields/Z2_field_operators.h>
#include <gudhi/Fields/Zp_field_operators.h>

using Gudhi::persistence_matrix::Cell;
using Gudhi::persistence_matrix::Cell_column_index;
using Gudhi::persistence_matrix::Cell_field_element;
using Gudhi::persistence_matrix::Chain_column_extra_properties;
using Gudhi::persistence_matrix::Column_dimension_holder;
using Gudhi::persistence_matrix::Column_types;
using Gudhi::persistence_matrix::Dummy_cell_column_index_mixin;
using Gudhi::persistence_matrix::Dummy_cell_field_element_mixin;
using Gudhi::persistence_matrix::Dummy_chain_properties;
using Gudhi::persistence_matrix::Dummy_dimension_holder;
using Gudhi::persistence_matrix::Dummy_row_access;
using Gudhi::persistence_matrix::New_cell_constructor;
using Gudhi::persistence_matrix::Pool_cell_constructor;
using Gudhi::persistence_matrix::Row_access;

using Zp = Gudhi::persistence_fields::Zp_field_operators<>;
using Z2 = Gudhi::persistence_fields::Z2_field_operators;

template <class Options>
struct Column_mini_matrix {
  using Option_list = Options;
  using Field_operators = typename Options::Field_coeff_operators;
  using Index = typename Options::Index;
  using ID_index = typename Options::Index;
  using Pos_index = typename Options::Index;
  using Dimension = typename Options::Dimension;
  using Element = typename std::conditional<Options::is_z2, bool, typename Field_operators::Element>::type;
  using Characteristic =
      typename std::conditional<Options::is_z2, unsigned int, typename Field_operators::Characteristic>::type;

  struct Matrix_row_tag;
  struct Matrix_column_tag;

  using Base_hook_matrix_row =
      boost::intrusive::list_base_hook<boost::intrusive::tag<Matrix_row_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::auto_unlink> >;
  using Base_hook_matrix_list_column =
      boost::intrusive::list_base_hook<boost::intrusive::tag<Matrix_column_tag>,
                                       boost::intrusive::link_mode<boost::intrusive::safe_link> >;
  using Base_hook_matrix_set_column =
      boost::intrusive::set_base_hook<boost::intrusive::tag<Matrix_column_tag>,
                                      boost::intrusive::link_mode<boost::intrusive::safe_link> >;
  struct Dummy_row_hook {};
  struct Dummy_column_hook {};

  using Row_hook = typename std::conditional<Options::has_row_access && Options::has_intrusive_rows,
                                                  Base_hook_matrix_row,
                                                  Dummy_row_hook
                                                 >::type;
  using Column_hook =
      typename std::conditional<Options::column_type == Column_types::INTRUSIVE_LIST,
                                Base_hook_matrix_list_column,
                                typename std::conditional<Options::column_type == Column_types::INTRUSIVE_SET,
                                                          Base_hook_matrix_set_column,
                                                          Dummy_column_hook
                                                         >::type
                               >::type;

  using Cell_column_index_option =
      typename std::conditional<Options::has_row_access, Cell_column_index<Index>, Dummy_cell_column_index_mixin>::type;
  using Cell_field_element_option = typename std::conditional<Options::is_z2,
                                                              Dummy_cell_field_element_mixin,
                                                              Cell_field_element<Element>
                                                             >::type;
  using Matrix_cell = Cell<Column_mini_matrix<Options> >;

  inline static New_cell_constructor<Matrix_cell> defaultCellConstructor;
  using Cell_constructor = New_cell_constructor<Matrix_cell>;

  struct Column_z2_settings {
    Column_z2_settings() : cellConstructor() {}
    Column_z2_settings([[maybe_unused]] Characteristic characteristic) : cellConstructor() {}

    Cell_constructor cellConstructor;
  };

  struct Column_zp_settings {
    Column_zp_settings() : operators(), cellConstructor() {}
    Column_zp_settings(Characteristic characteristic) : operators(characteristic), cellConstructor() {}

    Field_operators operators;
    Cell_constructor cellConstructor;
  };

  using Column_settings = typename std::conditional<Options::is_z2, Column_z2_settings, Column_zp_settings>::type;

  template <class Matrix_cell>
  struct RowCellComp {
    bool operator()(const Matrix_cell& c1, const Matrix_cell& c2) const {
      return c1.get_column_index() < c2.get_column_index();
    }
  };

  using Row =
      typename std::conditional<Options::has_intrusive_rows,
                                boost::intrusive::list<Matrix_cell,
                                                       boost::intrusive::constant_time_size<false>,
                                                       boost::intrusive::base_hook<Base_hook_matrix_row> >,
                                std::set<Matrix_cell, RowCellComp<Matrix_cell> >
                               >::type;

  using Row_container = typename std::conditional<Options::has_removable_rows,
                                                       std::map<ID_index, Row>,
                                                       std::vector<Row>
                                                      >::type;

  using Row_access_option = typename std::conditional<Options::has_row_access,
                                                      Row_access<Column_mini_matrix<Options> >,
                                                      Dummy_row_access
                                                     >::type;

  static const bool isNonBasic = !Options::is_basic;

  using Column_dimension_option =
      typename std::conditional<isNonBasic,
                                Column_dimension_holder<Column_mini_matrix<Options> >,
                                Dummy_dimension_holder
                               >::type;

  using Chain_column_option = typename std::conditional<isNonBasic && !Options::is_of_boundary_type,
                                                        Chain_column_extra_properties<Column_mini_matrix<Options> >,
                                                        Dummy_chain_properties
                                                       >::type;

  using Boundary = typename std::conditional<Options::is_z2,
                                                  std::initializer_list<ID_index>,
                                                  std::initializer_list<std::pair<ID_index, Element> >
                                                 >::type;
};

template <bool is_z2_only, Column_types col_type, bool has_row, bool rem_row, bool intr_row>
struct Base_col_options {
  using Field_coeff_operators = typename std::conditional<is_z2_only, Z2, Zp>::type;
  using Index = unsigned int;
  using Dimension = int;  // needs to be signed.

  static const bool is_basic = true;  // exists just for the tests
  static const bool is_of_boundary_type = true;

  static const bool is_z2 = is_z2_only;
  static const Column_types column_type = col_type;

  static const bool has_row_access = has_row;
  static const bool has_removable_rows = rem_row;   // ignored if has_row_access == false
  static const bool has_intrusive_rows = intr_row;  // ignored if has_row_access == false
};

template <bool is_z2_only, Column_types col_type, bool has_row, bool rem_row, bool intr_row>
struct Boundary_col_options {
  using Field_coeff_operators = typename std::conditional<is_z2_only, Z2, Zp>::type;
  using Index = unsigned int;
  using Dimension = int;  // needs to be signed.

  static const bool is_basic = false;  // exists just for the tests
  static const bool is_of_boundary_type = true;

  static const bool is_z2 = is_z2_only;
  static const Column_types column_type = col_type;

  static const bool has_row_access = has_row;
  static const bool has_removable_rows = rem_row;   // ignored if has_row_access == false
  static const bool has_intrusive_rows = intr_row;  // ignored if has_row_access == false
};

template <bool is_z2_only, Column_types col_type, bool has_row, bool rem_row, bool intr_row>
struct Chain_col_options {
  using Field_coeff_operators = typename std::conditional<is_z2_only, Z2, Zp>::type;
  using Index = unsigned int;
  using Dimension = int;  // needs to be signed.

  static const bool is_basic = false;  // exists just for the tests
  static const bool is_of_boundary_type = false;

  static const bool is_z2 = is_z2_only;
  static const Column_types column_type = col_type;

  static const bool has_row_access = has_row;
  static const bool has_removable_rows = rem_row;   // ignored if has_row_access == false
  static const bool has_intrusive_rows = intr_row;  // ignored if has_row_access == false
};

#endif  // PM_COLUMN_TESTS_MASTER_H
