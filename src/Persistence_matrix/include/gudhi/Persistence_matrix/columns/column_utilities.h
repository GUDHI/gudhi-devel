/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file column_utilities.h
 * @author Hannah Schreiber
 * @brief Contains helper methods for column addition and column hasher.
 */

#ifndef PM_COLUMN_UTILITIES_H
#define PM_COLUMN_UTILITIES_H

#include <cstddef>
#include <stdexcept>

#include <gudhi/persistence_matrix_options.h>

namespace Gudhi {
namespace persistence_matrix {

template <class Cell, typename Cell_iterator>
Cell* _get_cell(const Cell_iterator& itTarget)
{
  if constexpr (Cell::Master::Option_list::column_type == Column_types::INTRUSIVE_LIST ||
                Cell::Master::Option_list::column_type == Column_types::INTRUSIVE_SET) {
    return &*itTarget;
  } else {
    return *itTarget;
  }
}

// works only for ordered columns
template <class Column, class Cell_iterator, typename F1, typename F2, typename F3, typename F4>
void _generic_merge_cell_to_column(Column& targetColumn,
                                   Cell_iterator& itSource,
                                   typename Column::Column_support::iterator& itTarget,
                                   F1&& process_target,
                                   F2&& process_source,
                                   F3&& update_target1,
                                   F4&& update_target2,
                                   bool& pivotIsZeroed)
{
  typename Column::Cell* targetCell = _get_cell<typename Column::Cell>(itTarget);

  if (targetCell->get_row_index() < itSource->get_row_index()) {
    process_target(targetCell);
    ++itTarget;
  } else if (targetCell->get_row_index() > itSource->get_row_index()) {
    process_source(itSource, itTarget);
    ++itSource;
  } else {
    if constexpr (Column::Master::Option_list::is_z2) {
      //_multiply_*_and_add never enters here so not treated
      if constexpr (Column::Master::isNonBasic && !Column::Master::Option_list::is_of_boundary_type) {
        if (targetCell->get_row_index() == targetColumn.get_pivot()) pivotIsZeroed = true;
      }
      targetColumn._delete_cell(itTarget);
    } else {
      update_target1(targetCell->get_element(), itSource);
      if (targetCell->get_element() == Column::Field_operators::get_additive_identity()) {
        if constexpr (Column::Master::isNonBasic && !Column::Master::Option_list::is_of_boundary_type) {
          if (targetCell->get_row_index() == targetColumn.get_pivot()) pivotIsZeroed = true;
        }
        targetColumn._delete_cell(itTarget);
      } else {
        update_target2(targetCell);
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_cell(*targetCell);
        ++itTarget;
      }
    }
    ++itSource;
  }
}

// works only for ordered columns
template <class Column, class Cell_range, typename F1, typename F2, typename F3, typename F4, typename F5>
bool _generic_add_to_column(const Cell_range& source,
                            Column& targetColumn,
                            F1&& process_target,
                            F2&& process_source,
                            F3&& update_target1,
                            F4&& update_target2,
                            F5&& finish_target)
{
  bool pivotIsZeroed = false;

  auto& target = targetColumn.column_;
  auto itTarget = target.begin();
  auto itSource = source.begin();
  while (itTarget != target.end() && itSource != source.end()) {
    _generic_merge_cell_to_column(targetColumn, itSource, itTarget,
                                  process_target, process_source, update_target1, update_target2,
                                  pivotIsZeroed);
  }

  finish_target(itTarget);

  while (itSource != source.end()) {
    process_source(itSource, target.end());
    ++itSource;
  }

  return pivotIsZeroed;
}

template <class Column, class Cell_range>
bool _add_to_column(const Cell_range& source, Column& targetColumn)
{
  return _generic_add_to_column(
      source,
      targetColumn,
      [&]([[maybe_unused]] typename Column::Cell* cellTarget) {},
      [&](typename Cell_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        if constexpr (Column::Master::Option_list::is_z2) {
          targetColumn._insert_cell(itSource->get_row_index(), itTarget);
        } else {
          targetColumn._insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
        }
      },
      [&](typename Column::Field_element& targetElement, typename Cell_range::const_iterator& itSource) {
        if constexpr (!Column::Master::Option_list::is_z2)
          targetColumn.operators_->add_inplace(targetElement, itSource->get_element());
      },
      [&]([[maybe_unused]] typename Column::Cell* cellTarget) {},
      [&]([[maybe_unused]] typename Column::Column_support::iterator& itTarget) {}
    );
}

template <class Column, class Cell_range>
bool _multiply_target_and_add_to_column(const typename Column::Field_element& val,
                                        const Cell_range& source,
                                        Column& targetColumn)
{
  if (val == 0u) {
    if constexpr (Column::Master::isNonBasic && !Column::Master::Option_list::is_of_boundary_type) {
      throw std::invalid_argument("A chain column should not be multiplied by 0.");
      // this would not only mess up the base, but also the pivots stored.
    } else {
      targetColumn.clear();
    }
  }

  return _generic_add_to_column(
      source,
      targetColumn,
      [&](typename Column::Cell* cellTarget) {
        targetColumn.operators_->multiply_inplace(cellTarget->get_element(), val);
        // targetColumn.RA_opt::update_cell(*itTarget) produces an internal compiler error
        // even though it works in _generic_add_to_column... Probably because of the lambda.
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_cell(*cellTarget);
      },
      [&](typename Cell_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        targetColumn._insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
      },
      [&](typename Column::Field_element& targetElement, typename Cell_range::const_iterator& itSource) {
        targetColumn.operators_->multiply_and_add_inplace_front(targetElement, val, itSource->get_element());
      },
      [&]([[maybe_unused]] typename Column::Cell* cellTarget) {},
      [&](typename Column::Column_support::iterator& itTarget) {
        while (itTarget != targetColumn.column_.end()) {
          typename Column::Cell* targetCell = _get_cell<typename Column::Cell>(itTarget);
          targetColumn.operators_->multiply_inplace(targetCell->get_element(), val);
          if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_cell(*targetCell);
          itTarget++;
        }
      }
    );
}

template <class Column, class Cell_range>
bool _multiply_source_and_add_to_column(const typename Column::Field_element& val,
                                        const Cell_range& source,
                                        Column& targetColumn)
{
  if (val == 0u) {
    return false;
  }

  return _generic_add_to_column(
      source,
      targetColumn,
      []([[maybe_unused]] typename Column::Cell* cellTarget) {},
      [&](typename Cell_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        typename Column::Cell* cell =
            targetColumn._insert_cell(itSource->get_element(), itSource->get_row_index(), itTarget);
        targetColumn.operators_->multiply_inplace(cell->get_element(), val);
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_cell(*cell);
      },
      [&](typename Column::Field_element& targetElement, typename Cell_range::const_iterator& itSource) {
        targetColumn.operators_->multiply_and_add_inplace_back(itSource->get_element(), val, targetElement);
      },
      [&]([[maybe_unused]] typename Column::Cell* cellTarget) {},
      []([[maybe_unused]] typename Column::Column_support::iterator& itTarget) {});
}

// column has to be ordered (ie. not suited for unordered_map and heap) and contain the exact values
// (ie. not suited for vector and heap). A same column but ordered differently will have another hash value.
template <class Column>
std::size_t hash_column(const Column& column)
{
  std::size_t seed = 0;
  for (auto& cell : column) {
    seed ^= std::hash<unsigned int>()(cell.get_row_index() * static_cast<unsigned int>(cell.get_element())) +
            0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_COLUMN_UTILITIES_H
