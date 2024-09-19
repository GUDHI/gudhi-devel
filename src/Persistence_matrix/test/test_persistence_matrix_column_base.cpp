/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistence_matrix"
#include <boost/test/unit_test.hpp>

#include "pm_column_tests.h"
#include "pm_column_tests_boost_type_lists.h"

using option_name_list = mp_list_q<c_base_options>;

#ifdef PM_TEST_Z2
#ifdef PM_TEST_NO_ROW

using z2_no_row_access_columns = columns_list<z2_no_ra_option_list<option_name_list> >;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_constructors, Column, z2_no_row_access_columns) {
  column_test_common_constructors<Column>();
  column_test_base_boundary_constructors<Column>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_content_access, Column, z2_no_row_access_columns) {
  typename Column::Column_settings settings;
  std::vector<Column> matrix = build_column_matrix<Column>(settings);
  column_test_common_z2_content_access(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_operators, Column, z2_no_row_access_columns) {
  typename Column::Column_settings settings;
  std::vector<Column> matrix = build_column_matrix<Column>(settings);
  column_test_common_z2_operators(matrix);

  matrix = build_column_matrix<Column>(settings);
  column_test_base_z2_operators(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_other, Column, z2_no_row_access_columns) {
  column_test_base_boundary_z2_methods<Column>();
}

#else

#ifdef PM_TEST_REM_ROW
using z2_only_row_access_columns = columns_list<z2_only_ra_r_option_list<option_name_list> >;
#else
using z2_only_row_access_columns = columns_list<z2_only_ra_option_list<option_name_list> >;
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_with_row_access_content_access, Column, z2_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings;
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);
  column_test_common_z2_content_access(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_with_row_access_operators, Column, z2_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings;
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);
  column_test_common_z2_operators(matrix);

  matrix.clear();
  rows.clear();
  matrix = build_column_matrix<Column>(rows, settings);
  column_test_base_z2_operators(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_row_access_constructors, Column, z2_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings;
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);

  column_test_row_access_constructors(matrix, rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z2_column_other, Column, z2_only_row_access_columns) {
  column_test_base_boundary_z2_methods<Column>();
}

#endif
#endif

#ifdef PM_TEST_Z5
#ifdef PM_TEST_NO_ROW

using z5_no_row_access_columns = columns_list<z5_no_ra_option_list<option_name_list> >;

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_constructors, Column, z5_no_row_access_columns) {
  column_test_common_constructors<Column>();
  column_test_base_boundary_constructors<Column>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_content_access, Column, z5_no_row_access_columns) {
  typename Column::Column_settings settings(5);
  std::vector<Column> matrix = build_column_matrix<Column>(settings);
  column_test_common_z5_content_access(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_operators, Column, z5_no_row_access_columns) {
  typename Column::Column_settings settings(5);
  std::vector<Column> matrix = build_column_matrix<Column>(settings);
  column_test_common_z5_operators(matrix);

  matrix = build_column_matrix<Column>(settings);
  column_test_base_z5_operators(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_other, Column, z5_no_row_access_columns) {
  column_test_base_boundary_z5_methods<Column>();
}

#else

#ifdef PM_TEST_REM_ROW
using z5_only_row_access_columns = columns_list<z5_only_ra_r_option_list<option_name_list> >;
#else
using z5_only_row_access_columns = columns_list<z5_only_ra_option_list<option_name_list> >;
#endif

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_with_row_access_content_access, Column, z5_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings(5);
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);
  column_test_common_z5_content_access(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_with_row_access_operators, Column, z5_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings(5);
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);
  column_test_common_z5_operators(matrix);

  matrix.clear();
  rows.clear();
  matrix = build_column_matrix<Column>(rows, settings);
  column_test_base_z5_operators(matrix);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_column_row_access_constructors, Column, z5_only_row_access_columns) {
  typename Column::Master::Row_container rows;  // do not destroy before matrix
  typename Column::Column_settings settings(5);
  std::vector<Column> matrix = build_column_matrix<Column>(rows, settings);

  column_test_row_access_constructors(matrix, rows);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Base_z5_column_other, Column, z5_only_row_access_columns) {
  column_test_base_boundary_z5_methods<Column>();
}

#endif
#endif
