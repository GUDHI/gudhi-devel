/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PM_COLUMN_TESTS_H
#define PM_COLUMN_TESTS_H

#include <vector>
#include <set>

#include <boost/test/test_tools.hpp>

#include "pm_test_utilities.h"
#include "pm_column_tests_mastermatrix.h"

// assumes column was not modified since construction, ie no duplicated / erased values in heap or vector column
template <class Column>
std::vector<column_content<Column> > get_ordered_column_contents(std::vector<Column>& matrix) {
  std::vector<std::set<typename std::conditional<is_z2<Column>(),
                                                 unsigned int,
                                                 std::pair<unsigned int, typename Column::Field_element>
                                                >::type> > ordCol(matrix.size());
  for (unsigned int i = 0; i < matrix.size(); ++i) {
    Column& col = matrix[i];
    for (auto& entry : col) {
      if constexpr (is_z2<Column>()) {
        ordCol[i].insert(entry.get_row_index());
      } else {
        ordCol[i].insert({entry.get_row_index(), entry.get_element()});
      }
    }
  }
  return ordCol;
}

// assumes column was not modified since construction, ie no duplicated / erased values in heap or vector column
template <class Column>
std::vector<column_content<Column> > get_ordered_rows(std::vector<Column>& matrix) {
  std::vector<std::set<typename std::conditional<is_z2<Column>(),
                                                 unsigned int,
                                                 std::pair<unsigned int, typename Column::Field_element>
                                                >::type> > rows;
  for (Column& col : matrix) {
    for (auto& entry : col) {
      if (entry.get_row_index() >= rows.size()) rows.resize(entry.get_row_index() + 1);
      if constexpr (is_z2<Column>()) {
        rows[entry.get_row_index()].insert(entry.get_column_index());
      } else {
        rows[entry.get_row_index()].insert({entry.get_column_index(), entry.get_element()});
      }
    }
  }
  return rows;
}

template <class Column>
std::vector<Column> build_column_matrix(typename Column::Column_settings& settings) {
  std::vector<Column> matrix;

  if constexpr (is_z2<Column>()) {
    using cont = std::vector<unsigned int>;

    matrix.emplace_back(cont{0, 1, 3, 5}, 4, &settings);
    matrix.emplace_back(cont{0, 1, 2, 5, 6}, 4, &settings);
    matrix.emplace_back(cont{0, 1, 2, 5, 6}, 4, &settings);
    matrix.emplace_back(cont{}, 4, &settings);
    matrix.emplace_back(cont{0, 1, 3, 5}, 4, &settings);
    matrix.emplace_back(matrix[1]);
  } else {
    using cont = std::vector<std::pair<unsigned int, typename Column::Field_element> >;
    matrix.emplace_back(cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, 4, &settings);
    matrix.emplace_back(cont{{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}}, 4, &settings);
    matrix.emplace_back(cont{{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}}, 4, &settings);
    matrix.emplace_back(cont{}, 4, &settings);
    matrix.emplace_back(cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, 4, &settings);
    matrix.emplace_back(matrix[1]);
  }

  return matrix;
}

template <class Column, class Rows>
std::vector<Column> build_column_matrix(Rows& rows, typename Column::Column_settings& settings) {
  std::vector<Column> matrix;

  if constexpr (is_z2<Column>()) {
    using cont = std::vector<unsigned int>;
    matrix.emplace_back(0, cont{0, 1, 3, 5}, 4, &rows, &settings);
    matrix.emplace_back(1, cont{0, 1, 2, 5, 6}, 4, &rows, &settings);
    matrix.emplace_back(2, cont{0, 1, 2, 5, 6}, 4, &rows, &settings);
    matrix.emplace_back(3, cont{}, 4, &rows, &settings);
    matrix.emplace_back(4, cont{0, 1, 3, 5}, 4, &rows, &settings);
    matrix.emplace_back(matrix[1], 5, &rows);
  } else {
    using cont = std::vector<std::pair<unsigned int, typename Column::Field_element> >;
    matrix.emplace_back(0, cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, 4, &rows, &settings);
    matrix.emplace_back(1, cont{{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}}, 4, &rows, &settings);
    matrix.emplace_back(2, cont{{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}}, 4, &rows, &settings);
    matrix.emplace_back(3, cont{}, 4, &rows, &settings);
    matrix.emplace_back(4, cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, 4, &rows, &settings);
    matrix.emplace_back(matrix[1], 5, &rows);
  }

  return matrix;
}

template <class Column, class Rows>
std::vector<Column> build_base_boundary_column_matrix(Rows& rows, typename Column::Column_settings& settings) {
  std::vector<Column> matrix;

  if constexpr (is_z2<Column>()) {
    using cont = std::vector<unsigned int>;
    matrix.emplace_back(0, cont{0, 1, 3, 5}, &rows, &settings);
    matrix.emplace_back(1, cont{0, 1, 2, 5, 6}, &rows, &settings);
    matrix.emplace_back(2, cont{0, 1, 2, 5, 6}, &rows, &settings);
    matrix.emplace_back(3, cont{}, &rows, &settings);
    matrix.emplace_back(4, cont{0, 1, 3, 5}, &rows, &settings);
    matrix.emplace_back(matrix[1], 5, &rows);
  } else {
    using cont = std::vector<std::pair<unsigned int, typename Column::Field_element> >;
    matrix.emplace_back(0, cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, &rows, &settings);
    matrix.emplace_back(1, cont{{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}}, &rows, &settings);
    matrix.emplace_back(2, cont{{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}}, &rows, &settings);
    matrix.emplace_back(3, cont{}, &rows, &settings);
    matrix.emplace_back(4, cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, &rows, &settings);
    matrix.emplace_back(matrix[1], 5, &rows);
  }

  return matrix;
}

template <class Column>
std::vector<Column> build_base_boundary_column_matrix(typename Column::Column_settings& settings) {
  std::vector<Column> matrix;

  if constexpr (is_z2<Column>()) {
    using cont = std::vector<unsigned int>;
    matrix.emplace_back(cont{0, 1, 3, 5}, &settings);
    matrix.emplace_back(cont{0, 1, 2, 5, 6}, &settings);
    matrix.emplace_back(cont{0, 1, 2, 5, 6}, &settings);
    matrix.emplace_back(cont{}, &settings);
    matrix.emplace_back(cont{0, 1, 3, 5}, &settings);
    matrix.emplace_back(matrix[1]);
  } else {
    using cont = std::vector<std::pair<unsigned int, typename Column::Field_element> >;
    matrix.emplace_back(cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, &settings);
    matrix.emplace_back(cont{{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}}, &settings);
    matrix.emplace_back(cont{{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}}, &settings);
    matrix.emplace_back(cont{}, &settings);
    matrix.emplace_back(cont{{0, 1}, {1, 2}, {3, 3}, {5, 4}}, &settings);
    matrix.emplace_back(matrix[1]);
  }

  return matrix;
}

template <class Column>
std::vector<std::vector<typename std::conditional<
    is_z2<Column>(),
    unsigned int,
    std::pair<unsigned int, typename Column::Field_element>
   >::type> > build_rows() {
  using entry_type = typename std::conditional<is_z2<Column>(),
                                               unsigned int,
                                               std::pair<unsigned int, typename Column::Field_element>
                                              >::type;
  using Container = std::vector<entry_type>;
  std::vector<Container> rows(7);

  if constexpr (is_z2<Column>()) {
    rows[0] = {0, 1, 2, 4, 5};
    rows[1] = {0, 1, 2, 4, 5};
    rows[2] = {1, 2, 5};
    rows[3] = {0, 4};
    rows[4] = {};
    rows[5] = {0, 1, 2, 4, 5};
    rows[6] = {1, 2, 5};
  } else {
    rows[0] = {{0, 1}, {1, 4}, {2, 1}, {4, 1}, {5, 4}};
    rows[1] = {{0, 2}, {1, 2}, {2, 3}, {4, 2}, {5, 2}};
    rows[2] = {{1, 1}, {2, 4}, {5, 1}};
    rows[3] = {{0, 3}, {4, 3}};
    rows[4] = {};
    rows[5] = {{0, 4}, {1, 1}, {2, 4}, {4, 4}, {5, 1}};
    rows[6] = {{1, 1}, {2, 4}, {5, 1}};
  }

  return rows;
}

template <class Column>
std::vector<std::vector<typename std::conditional<
    is_z2<Column>(),
    unsigned int,
    std::pair<unsigned int, typename Column::Field_element>
   >::type> > build_column_values() {
  using entry_type = typename std::conditional<is_z2<Column>(),
                                               unsigned int,
                                               std::pair<unsigned int, typename Column::Field_element>
                                              >::type;
  using Container = std::vector<entry_type>;
  std::vector<Container> columns(6);

  if constexpr (is_z2<Column>()) {
    columns[0] = {0, 1, 3, 5};
    columns[1] = {0, 1, 2, 5, 6};
    columns[2] = {0, 1, 2, 5, 6};
    columns[3] = {};
    columns[4] = {0, 1, 3, 5};
    columns[5] = {0, 1, 2, 5, 6};
  } else {
    columns[0] = {{0, 1}, {1, 2}, {3, 3}, {5, 4}};
    columns[1] = {{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}};
    columns[2] = {{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}};
    columns[3] = {};
    columns[4] = {{0, 1}, {1, 2}, {3, 3}, {5, 4}};
    columns[5] = {{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}};
  }

  return columns;
}

template <class Column>
void column_test_common_constructors() {
  using entry_type = typename std::conditional<is_z2<Column>(),
                                               unsigned int,
                                               std::pair<unsigned int, typename Column::Field_element>
                                              >::type;
  using Container = std::vector<entry_type>;

  Container cont1, cont2;
  typename Column::Column_settings settings(5);

  if constexpr (is_z2<Column>()) {
    cont1 = {0, 2, 4};
    cont2 = {0, 5, 6};
  } else {
    cont1 = {{0, 1u}, {2, 2u}, {4, 3u}};
    cont2 = {{0, 1u}, {5, 2u}, {6, 3u}};
  }

  Column emptyCol(&settings);
  BOOST_CHECK_EQUAL(emptyCol.size(), 0);

  Column col(cont1, 2, &settings);
  BOOST_CHECK_EQUAL(col.size(), 3);

  std::vector<int> rows;  // type doesn't matter, as the row option is not enabled.
  Column rowCol(2, cont2, 2, &rows, &settings);
  BOOST_CHECK_EQUAL(rowCol.size(), 3);

  Column copyCol(col);
  BOOST_CHECK_EQUAL(copyCol.size(), 3);

  Column rowCopyCol(rowCol, 4, &rows);
  BOOST_CHECK_EQUAL(rowCopyCol.size(), 3);

  Column moveCol(std::move(col));
  BOOST_CHECK_EQUAL(moveCol.size(), 3);
  BOOST_CHECK_EQUAL(col.size(), 0);
  BOOST_CHECK(copyCol.get_content() == moveCol.get_content());

  Column assignCol = copyCol;
  BOOST_CHECK_EQUAL(assignCol.size(), 3);

  swap(copyCol, rowCol);
  BOOST_CHECK_EQUAL(copyCol.size(), 3);
  BOOST_CHECK_EQUAL(rowCol.size(), 3);
  BOOST_CHECK(copyCol.get_content() == rowCopyCol.get_content());
  BOOST_CHECK(rowCol.get_content() == assignCol.get_content());

  BOOST_CHECK_EQUAL(rows.size(), 0);
}

template <class Column, typename entry_type>
void column_test_common_content_access(Column& col,
                                       const std::set<entry_type>& setcont,
                                       const std::vector<typename Column::Field_element>& veccont) {
  BOOST_CHECK(get_column_content_via_iterators(col) == setcont);
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(col.size(), setcont.size());
  }

  if (veccont.empty())
    BOOST_CHECK(col.is_empty());
  else
    BOOST_CHECK(!col.is_empty());

  for (unsigned int i = 0; i < veccont.size(); ++i) {
    if (veccont[i] == 0u)
      BOOST_CHECK(!col.is_non_zero(i));
    else
      BOOST_CHECK(col.is_non_zero(i));
  }
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_common_z5_content_access(std::vector<Column>& matrix) {
  std::set<std::pair<unsigned int, typename Column::Field_element> > setcont;
  std::vector<typename Column::Field_element> veccont;

  setcont = {{0, 1}, {1, 2}, {3, 3}, {5, 4}};
  veccont = {1, 2, 0, 3, 0, 4};
  column_test_common_content_access(matrix[0], setcont, veccont);

  setcont = {{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}};
  veccont = {4, 2, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[1], setcont, veccont);

  setcont = {{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}};
  veccont = {1, 3, 4, 0, 0, 4, 4};
  column_test_common_content_access(matrix[2], setcont, veccont);

  setcont = {};
  veccont = {};
  column_test_common_content_access(matrix[3], setcont, veccont);

  setcont = {{0, 1}, {1, 2}, {3, 3}, {5, 4}};
  veccont = {1, 2, 0, 3, 0, 4};
  column_test_common_content_access(matrix[4], setcont, veccont);

  setcont = {{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}};
  veccont = {4, 2, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[5], setcont, veccont);
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_common_z2_content_access(std::vector<Column>& matrix) {
  std::set<unsigned int> setcont;
  std::vector<typename Column::Field_element> veccont;

  setcont = {0, 1, 3, 5};
  veccont = {1, 1, 0, 1, 0, 1};
  column_test_common_content_access(matrix[0], setcont, veccont);

  setcont = {0, 1, 2, 5, 6};
  veccont = {1, 1, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[1], setcont, veccont);

  setcont = {0, 1, 2, 5, 6};
  veccont = {1, 1, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[2], setcont, veccont);

  setcont = {};
  veccont = {};
  column_test_common_content_access(matrix[3], setcont, veccont);

  setcont = {0, 1, 3, 5};
  veccont = {1, 1, 0, 1, 0, 1};
  column_test_common_content_access(matrix[4], setcont, veccont);

  setcont = {0, 1, 2, 5, 6};
  veccont = {1, 1, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[5], setcont, veccont);
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_common_z5_operators(std::vector<Column>& matrix) {
  std::set<std::pair<unsigned int, typename Column::Field_element> > setcont;
  std::vector<typename Column::Field_element> veccont;

  matrix[0] += matrix[1];
  matrix[1] += matrix[2];

  setcont = {{1, 4}, {2, 1}, {3, 3}, {6, 1}};
  veccont = {0, 4, 1, 3, 0, 0, 1};
  column_test_common_content_access(matrix[0], setcont, veccont);
  BOOST_CHECK(matrix[1] == matrix[1]);
  BOOST_CHECK(matrix[1] == matrix[3]);
  BOOST_CHECK(matrix[1] < matrix[0]);
  if constexpr (Column::Master::Option_list::column_type == Column_types::HEAP) {
    BOOST_CHECK(matrix[0] < matrix[2]);  // order different for heap columns
  } else {
    BOOST_CHECK(matrix[2] < matrix[0]);
  }
  BOOST_CHECK(!(matrix[0] < matrix[0]));

  matrix[0] *= 4;
  matrix[1] *= 2;
  matrix[2] *= 1;

  setcont = {{1, 1}, {2, 4}, {3, 2}, {6, 4}};
  veccont = {0, 1, 4, 2, 0, 0, 4};
  column_test_common_content_access(matrix[0], setcont, veccont);
  setcont = {{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}};
  veccont = {1, 3, 4, 0, 0, 4, 4};
  column_test_common_content_access(matrix[2], setcont, veccont);
  BOOST_CHECK(matrix[1] == matrix[1]);
  BOOST_CHECK(matrix[1] == matrix[3]);
  BOOST_CHECK(matrix[1] < matrix[0]);
  if constexpr (Column::Master::Option_list::column_type == Column_types::HEAP) {
    BOOST_CHECK(matrix[0] < matrix[2]);  // order different for heap columns
  } else {
    BOOST_CHECK(matrix[2] < matrix[0]);
  }
  BOOST_CHECK(!(matrix[0] < matrix[0]));

  // this = v * this + column
  matrix[4].multiply_target_and_add(4, matrix[5]);
  setcont = {{0, 3}, {2, 1}, {3, 2}, {5, 2}, {6, 1}};
  veccont = {3, 0, 1, 2, 0, 2, 1};
  column_test_common_content_access(matrix[4], setcont, veccont);
  // this = this + column * v
  matrix[5].multiply_source_and_add(matrix[3], 3);
  setcont = {{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}};
  veccont = {4, 2, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[5], setcont, veccont);
  // this = this + column * v
  matrix[5].multiply_source_and_add(matrix[4], 3);
  setcont = {{0, 3}, {1, 2}, {2, 4}, {3, 1}, {5, 2}, {6, 4}};
  veccont = {3, 2, 4, 1, 0, 2, 4};
  column_test_common_content_access(matrix[5], setcont, veccont);
  // this = v * this + column
  matrix[3].multiply_target_and_add(4, matrix[5]);
  setcont = {{0, 3}, {1, 2}, {2, 4}, {3, 1}, {5, 2}, {6, 4}};
  veccont = {3, 2, 4, 1, 0, 2, 4};
  column_test_common_content_access(matrix[3], setcont, veccont);
  BOOST_CHECK(matrix[3] == matrix[5]);
  BOOST_CHECK(!(matrix[5] == matrix[4]));
  BOOST_CHECK(!(matrix[5] < matrix[3]));
  BOOST_CHECK(!(matrix[3] < matrix[5]));
  if constexpr (Column::Master::Option_list::column_type == Column_types::HEAP) {
    BOOST_CHECK(matrix[4] < matrix[5]);
    BOOST_CHECK(matrix[4] < matrix[3]);
  } else {
    BOOST_CHECK(matrix[5] < matrix[4]);
    BOOST_CHECK(matrix[3] < matrix[4]);
  }
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_common_z2_operators(std::vector<Column>& matrix) {
  std::set<unsigned int> setcont;
  std::vector<typename Column::Field_element> veccont;

  matrix[0] += matrix[1];
  matrix[1] += matrix[2];

  setcont = {2, 3, 6};
  veccont = {0, 0, 1, 1, 0, 0, 1};
  column_test_common_content_access(matrix[0], setcont, veccont);
  BOOST_CHECK(matrix[1] == matrix[1]);
  BOOST_CHECK(matrix[1] == matrix[3]);
  BOOST_CHECK(matrix[1] < matrix[0]);
  if constexpr (Column::Master::Option_list::column_type == Column_types::HEAP) {
    BOOST_CHECK(matrix[0] < matrix[2]);  // order different for heap columns
  } else {
    BOOST_CHECK(matrix[2] < matrix[0]);
  }
  BOOST_CHECK(!(matrix[0] < matrix[0]));

  matrix[0] *= 5;
  matrix[2] *= 1;

  setcont = {2, 3, 6};
  veccont = {0, 0, 1, 1, 0, 0, 1};
  column_test_common_content_access(matrix[0], setcont, veccont);
  setcont = {0, 1, 2, 5, 6};
  veccont = {1, 1, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[2], setcont, veccont);
  BOOST_CHECK(!(matrix[0] == matrix[2]));
  BOOST_CHECK(matrix[2] == matrix[2]);
  BOOST_CHECK(matrix[1] < matrix[0]);
  if constexpr (Column::Master::Option_list::column_type == Column_types::HEAP) {
    BOOST_CHECK(matrix[0] < matrix[2]);  // order different for heap columns
  } else {
    BOOST_CHECK(matrix[2] < matrix[0]);
  }
  BOOST_CHECK(!(matrix[0] < matrix[0]));

  // this = v * this + column
  matrix[4].multiply_target_and_add(3, matrix[5]);
  setcont = {2, 3, 6};
  veccont = {0, 0, 1, 1, 0, 0, 1};
  column_test_common_content_access(matrix[4], setcont, veccont);
  BOOST_CHECK(matrix[4] == matrix[0]);
  // this = this + column * v
  matrix[5].multiply_source_and_add(matrix[3], 3);
  setcont = {0, 1, 2, 5, 6};
  veccont = {1, 1, 1, 0, 0, 1, 1};
  column_test_common_content_access(matrix[5], setcont, veccont);
  BOOST_CHECK(matrix[5] == matrix[2]);
  // this = this + column * v
  matrix[5].multiply_source_and_add(matrix[4], 3);
  setcont = {0, 1, 3, 5};
  veccont = {1, 1, 0, 1, 0, 1, 0};
  column_test_common_content_access(matrix[5], setcont, veccont);
  // this = v * this + column
  matrix[3].multiply_target_and_add(3, matrix[5]);
  setcont = {0, 1, 3, 5};
  veccont = {1, 1, 0, 1, 0, 1, 0};
  column_test_common_content_access(matrix[3], setcont, veccont);
  BOOST_CHECK(matrix[3] == matrix[5]);
  BOOST_CHECK(!(matrix[5] < matrix[3]));
  BOOST_CHECK(!(matrix[3] < matrix[5]));
  BOOST_CHECK(matrix[5] < matrix[4]);
}

// common: assumes that matrix was build with build_column_matrix(rows) and not modified since.
// base/boundary special: assumes that matrix was build with build_base_boundary_column_matrix(rows) and not modified
// since.
template <class Column, class Row_container>
void column_test_row_access_constructors(std::vector<Column>& matrix, Row_container& rows) {
  test_matrix_equality<Column>(build_rows<Column>(), get_ordered_rows(matrix));
  test_matrix_equality<Column>(build_column_values<Column>(), get_ordered_column_contents(matrix));

  Column col(matrix[0], 6, &rows);
  for (auto& r : rows) {
    if constexpr (Column::Master::Option_list::has_removable_rows) {
      if (!r.second.empty()) {
        auto& entry = *r.second.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    } else {
      if (!r.empty()) {
        auto& entry = *r.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    }
  }

  Column moveCol(std::move(col));
  BOOST_CHECK_EQUAL(moveCol.size(), 4);
  BOOST_CHECK_EQUAL(col.size(), 0);
  BOOST_CHECK(matrix[0] == moveCol);
  for (auto& r : rows) {
    if constexpr (Column::Master::Option_list::has_removable_rows) {
      if (!r.second.empty()) {
        auto& entry = *r.second.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    } else {
      if (!r.empty()) {
        auto& entry = *r.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    }
  }

  swap(col, moveCol);
  BOOST_CHECK_EQUAL(moveCol.size(), 0);
  BOOST_CHECK_EQUAL(col.size(), 4);
  BOOST_CHECK(matrix[0] == col);
  for (auto& r : rows) {
    if constexpr (Column::Master::Option_list::has_removable_rows) {
      if (!r.second.empty()) {
        auto& entry = *r.second.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    } else {
      if (!r.empty()) {
        auto& entry = *r.rbegin();
        if (entry.get_row_index() == 0 || entry.get_row_index() == 1 || entry.get_row_index() == 3 ||
            entry.get_row_index() == 5) {
          BOOST_CHECK_EQUAL(entry.get_column_index(), 6);
        }
      }
    }
  }
}

template <class Column>
void column_test_base_boundary_constructors() {
  using entry_type = typename std::conditional<is_z2<Column>(),
                                               unsigned int,
                                               std::pair<unsigned int, typename Column::Field_element>
                                              >::type;
  using Container = std::vector<entry_type>;

  typename Column::Column_settings settings(5);

  Container cont1, cont2;

  if constexpr (is_z2<Column>()) {
    cont1 = {0, 2, 4};
    cont2 = {0, 5, 6};
  } else {
    cont1 = {{0, 1u}, {2, 2u}, {4, 3u}};
    cont2 = {{0, 1u}, {5, 2u}, {6, 3u}};
  }

  Column col(cont1, &settings);
  BOOST_CHECK_EQUAL(col.size(), 3);

  std::vector<int> rows;  // type doesn't matter, as the row option is not enabled.
  Column rowCol(2, cont2, &rows, &settings);
  BOOST_CHECK_EQUAL(rowCol.size(), 3);

  BOOST_CHECK(!(col == rowCol));
}

template <class Column>
void column_test_base_boundary_z5_methods() {
  std::vector<typename Column::Field_element> veccont;
  typename Column::Column_settings settings(5);

  Column col(
      std::vector<std::pair<unsigned int, typename Column::Field_element> >{
          {0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}},
      &settings);
  veccont = {4, 2, 1, 0, 0, 1, 1};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 5);

  std::vector<unsigned int> permutation{0, 5, 2, 3, 1, 4, 6};
  col.reorder(permutation);
  veccont = {4, 0, 1, 0, 1, 2, 1};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 5);

  col.clear(2);
  col.clear(6);
  veccont = {4, 0, 0, 0, 1, 2, 0};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 3);

  col.clear();
  veccont = {};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 0);
}

template <class Column>
void column_test_base_boundary_z2_methods() {
  std::vector<typename Column::Field_element> veccont;
  typename Column::Column_settings settings;

  Column col(std::vector<unsigned int>{0, 1, 2, 5, 6}, &settings);
  veccont = {1, 1, 1, 0, 0, 1, 1};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 5);

  std::vector<unsigned int> permutation{0, 5, 2, 3, 1, 4, 6};
  col.reorder(permutation);
  veccont = {1, 0, 1, 0, 1, 1, 1};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 5);

  col.clear(2);
  col.clear(6);
  veccont = {1, 0, 0, 0, 1, 1, 0};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 3);

  col.clear();
  veccont = {};
  BOOST_CHECK(col.get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(col.size(), 0);
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_base_z5_operators(std::vector<Column>& matrix) {
  using Entry = typename Column::Entry;
  std::set<Entry> setcont;
  std::vector<typename Column::Field_element> veccont;

  Entry entry(0);
  entry.set_element(4);
  setcont.insert(entry);
  entry = Entry(1);
  entry.set_element(2);
  setcont.insert(entry);
  entry = Entry(2);
  entry.set_element(1);
  setcont.insert(entry);
  entry = Entry(5);
  entry.set_element(1);
  setcont.insert(entry);
  entry = Entry(6);
  entry.set_element(1);
  setcont.insert(entry);
  matrix[0] += setcont;

  veccont = {0, 4, 1, 3, 0, 0, 1};
  BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[0].size(), 4);
  }

  setcont.clear();
  entry = Entry(0);
  entry.set_element(1);
  setcont.insert(entry);
  entry = Entry(1);
  entry.set_element(3);
  setcont.insert(entry);
  entry = Entry(2);
  entry.set_element(4);
  setcont.insert(entry);
  entry = Entry(5);
  entry.set_element(4);
  setcont.insert(entry);
  entry = Entry(6);
  entry.set_element(4);
  setcont.insert(entry);
  matrix[1] += setcont;

  veccont = {};
  BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[1].size(), 0);
  }

  // this = v * this + column
  matrix[4].multiply_target_and_add(4, setcont);
  veccont = {0, 1, 4, 2, 0, 0, 4};
  BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[4].size(), 4);
  }
  // this = this + column * v
  setcont = {};
  matrix[5].multiply_source_and_add(setcont, 3);
  veccont = {4, 2, 1, 0, 0, 1, 1};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 5);
  }
  // this = this + column * v
  setcont.clear();
  entry = Entry(0);
  entry.set_element(3);
  setcont.insert(entry);
  entry = Entry(2);
  entry.set_element(1);
  setcont.insert(entry);
  entry = Entry(3);
  entry.set_element(2);
  setcont.insert(entry);
  entry = Entry(5);
  entry.set_element(2);
  setcont.insert(entry);
  entry = Entry(6);
  entry.set_element(1);
  setcont.insert(entry);
  matrix[5].multiply_source_and_add(setcont, 3);
  veccont = {3, 2, 4, 1, 0, 2, 4};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 6);
  }
  // this = v * this + column
  setcont.clear();
  entry = Entry(0);
  entry.set_element(3);
  setcont.insert(entry);
  entry = Entry(1);
  entry.set_element(2);
  setcont.insert(entry);
  entry = Entry(2);
  entry.set_element(4);
  setcont.insert(entry);
  entry = Entry(3);
  entry.set_element(1);
  setcont.insert(entry);
  entry = Entry(5);
  entry.set_element(2);
  setcont.insert(entry);
  entry = Entry(6);
  entry.set_element(4);
  setcont.insert(entry);
  matrix[3].multiply_target_and_add(4, setcont);
  veccont = {3, 2, 4, 1, 0, 2, 4};
  BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[3].size(), 6);
  }
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_base_z2_operators(std::vector<Column>& matrix) {
  std::set<typename Column::Entry> setcont;
  std::vector<bool> veccont;

  setcont = {0, 1, 2, 5, 6};
  matrix[0] += setcont;

  veccont = {0, 0, 1, 1, 0, 0, 1};
  BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[0].size(), 3);
  }

  setcont = {0, 1, 2, 5, 6};
  matrix[1] += setcont;

  veccont = {};
  BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[1].size(), 0);
  }

  // this = v * this + column
  matrix[4].multiply_target_and_add(1, setcont);
  veccont = {0, 0, 1, 1, 0, 0, 1};
  BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[4].size(), 3);
  }
  // this = this + column * v
  setcont = {};
  matrix[5].multiply_source_and_add(setcont, 1);
  veccont = {1, 1, 1, 0, 0, 1, 1};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 5);
  }
  // this = this + column * v
  setcont = {2, 3, 6};
  matrix[5].multiply_source_and_add(setcont, 1);
  veccont = {1, 1, 0, 1, 0, 1, 0};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 4);
  }
  // this = v * this + column
  setcont = {0, 1, 3, 5};
  matrix[3].multiply_target_and_add(0, setcont);
  veccont = {1, 1, 0, 1, 0, 1, 0};
  BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[3].size(), 4);
  }
}

template <class Column>
void column_test_boundary_chain_methods(std::vector<Column>& matrix) {
  BOOST_CHECK_EQUAL(matrix[0].get_pivot(), 5);
  BOOST_CHECK_EQUAL(matrix[1].get_pivot(), 6);
  BOOST_CHECK_EQUAL(matrix[2].get_pivot(), 6);
  BOOST_CHECK_EQUAL(matrix[3].get_pivot(), -1);
  BOOST_CHECK_EQUAL(matrix[4].get_pivot(), 5);
  BOOST_CHECK_EQUAL(matrix[5].get_pivot(), 6);

  if constexpr (!is_z2<Column>()) {
    BOOST_CHECK_EQUAL(matrix[0].get_pivot_value(), 4u);
    BOOST_CHECK_EQUAL(matrix[1].get_pivot_value(), 1u);
    BOOST_CHECK_EQUAL(matrix[2].get_pivot_value(), 4u);
    BOOST_CHECK_EQUAL(matrix[3].get_pivot_value(), 0u);
    BOOST_CHECK_EQUAL(matrix[4].get_pivot_value(), 4u);
    BOOST_CHECK_EQUAL(matrix[5].get_pivot_value(), 1u);
  }

  BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[3].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[4].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[5].get_dimension(), 4);
}

template <class Column>
void column_test_boundary_methods(std::vector<Column>& matrix) {
  BOOST_CHECK_EQUAL(matrix[0].get_dimension(), 3);
  BOOST_CHECK_EQUAL(matrix[1].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[2].get_dimension(), 4);
  BOOST_CHECK_EQUAL(matrix[3].get_dimension(), 0);
  BOOST_CHECK_EQUAL(matrix[4].get_dimension(), 3);
  BOOST_CHECK_EQUAL(matrix[5].get_dimension(), 4);
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_boundary_z5_operators(std::vector<Column>& matrix) {
  using cont = std::vector<std::pair<unsigned int, typename Column::Field_element> >;

  std::vector<typename Column::Field_element> veccont;
  typename Column::Column_settings settings(5);

  const Column col0(cont{{0, 4}, {1, 2}, {2, 1}, {5, 1}, {6, 1}}, &settings);
  matrix[0] += col0;

  veccont = {0, 4, 1, 3, 0, 0, 1};
  BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[0].size(), 4);
  }

  const Column col1(cont{{0, 1}, {1, 3}, {2, 4}, {5, 4}, {6, 4}}, &settings);
  matrix[1] += col1;

  veccont = {};
  BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[1].size(), 0);
  }

  // this = v * this + column
  matrix[4].multiply_target_and_add(4, col1);
  veccont = {0, 1, 4, 2, 0, 0, 4};
  BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[4].size(), 4);
  }
  // this = this + column * v
  const Column col2(cont{}, &settings);
  matrix[5].multiply_source_and_add(col2, 3);
  veccont = {4, 2, 1, 0, 0, 1, 1};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  BOOST_CHECK_EQUAL(matrix[5].size(), 5);
  // this = this + column * v
  const Column col3(cont{{0, 3}, {2, 1}, {3, 2}, {5, 2}, {6, 1}}, &settings);
  matrix[5].multiply_source_and_add(col3, 3);
  veccont = {3, 2, 4, 1, 0, 2, 4};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 6);
  }
  // this = v * this + column
  const Column col4(cont{{0, 3}, {1, 2}, {2, 4}, {3, 1}, {5, 2}, {6, 4}}, &settings);
  matrix[3].multiply_target_and_add(4, col4);
  veccont = {3, 2, 4, 1, 0, 2, 4};
  BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[3].size(), 6);
  }
}

// assumes that matrix was build with build_column_matrix and was not modified since.
template <class Column>
void column_test_boundary_z2_operators(std::vector<Column>& matrix) {
  using cont = std::vector<unsigned int>;

  std::vector<bool> veccont;
  typename Column::Column_settings settings;

  const Column col0(cont{0, 1, 2, 5, 6}, &settings);
  matrix[0] += col0;

  veccont = {0, 0, 1, 1, 0, 0, 1};
  BOOST_CHECK(matrix[0].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[0].size(), 3);
  }

  const Column col1(cont{0, 1, 2, 5, 6}, &settings);
  matrix[1] += col1;

  veccont = {};
  BOOST_CHECK(matrix[1].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[1].size(), 0);
  }

  // this = v * this + column
  matrix[4].multiply_target_and_add(3, col1);
  veccont = {0, 0, 1, 1, 0, 0, 1};
  BOOST_CHECK(matrix[4].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[4].size(), 3);
  }
  // this = this + column * v
  const Column col2(cont{}, &settings);
  matrix[5].multiply_source_and_add(col2, 3);
  veccont = {1, 1, 1, 0, 0, 1, 1};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 5);
  }
  // this = this + column * v
  const Column col3(cont{2, 3, 6}, &settings);
  matrix[5].multiply_source_and_add(col3, 3);
  veccont = {1, 1, 0, 1, 0, 1, 0};
  BOOST_CHECK(matrix[5].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[5].size(), 4);
  }
  // this = v * this + column
  const Column col4(cont{0, 1, 3, 5}, &settings);
  matrix[3].multiply_target_and_add(4, col4);
  veccont = {1, 1, 0, 1, 0, 1, 0};
  BOOST_CHECK(matrix[3].get_content(veccont.size()) == veccont);
  if constexpr (Column::Master::Option_list::column_type != Column_types::HEAP) {
    BOOST_CHECK_EQUAL(matrix[3].size(), 4);
  }
}

template <class Column>
void column_test_chain_methods() {
  typename Column::Column_settings settings(5);

  Column col(&settings);

  BOOST_CHECK(!col.is_paired());
  BOOST_CHECK(col.get_paired_chain_index() == Column::Master::template get_null_value<typename Column::Index>());

  col.unassign_paired_chain();
  BOOST_CHECK(!col.is_paired());
  BOOST_CHECK(col.get_paired_chain_index() == Column::Master::template get_null_value<typename Column::Index>());

  col.assign_paired_chain(2);
  BOOST_CHECK(col.is_paired());
  BOOST_CHECK(col.get_paired_chain_index() == 2);

  col.unassign_paired_chain();
  BOOST_CHECK(!col.is_paired());
  BOOST_CHECK(col.get_paired_chain_index() == Column::Master::template get_null_value<typename Column::Index>());
}

#endif  // PM_COLUMN_TESTS_H
