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

template <class Entry, typename Entry_iterator>
Entry* _get_entry(const Entry_iterator& itTarget)
{
  if constexpr (Entry::Master::Option_list::column_type == Column_types::INTRUSIVE_LIST ||
                Entry::Master::Option_list::column_type == Column_types::INTRUSIVE_SET) {
    return &*itTarget;
  } else {
    return *itTarget;
  }
}

// works only for ordered columns
template <class Column, class Entry_iterator, typename F1, typename F2, typename F3, typename F4>
void _generic_merge_entry_to_column(Column& targetColumn,
                                    Entry_iterator& itSource,
                                    typename Column::Column_support::iterator& itTarget,
                                    F1&& process_target,
                                    F2&& process_source,
                                    F3&& update_target1,
                                    F4&& update_target2,
                                    bool& pivotIsZeroed)
{
  typename Column::Entry* targetEntry = _get_entry<typename Column::Entry>(itTarget);

  if (targetEntry->get_row_index() < itSource->get_row_index()) {
    process_target(targetEntry);
    ++itTarget;
  } else if (targetEntry->get_row_index() > itSource->get_row_index()) {
    process_source(itSource, itTarget);
    ++itSource;
  } else {
    if constexpr (Column::Master::Option_list::is_z2) {
      //_multiply_*_and_add never enters here so not treated
      if constexpr (Column::Master::isNonBasic && !Column::Master::Option_list::is_of_boundary_type) {
        if (targetEntry->get_row_index() == targetColumn.get_pivot()) pivotIsZeroed = true;
      }
      targetColumn._delete_entry(itTarget);
    } else {
      update_target1(targetEntry->get_element(), itSource);
      if (targetEntry->get_element() == Column::Field_operators::get_additive_identity()) {
        if constexpr (Column::Master::isNonBasic && !Column::Master::Option_list::is_of_boundary_type) {
          if (targetEntry->get_row_index() == targetColumn.get_pivot()) pivotIsZeroed = true;
        }
        targetColumn._delete_entry(itTarget);
      } else {
        update_target2(targetEntry);
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_entry(*targetEntry);
        ++itTarget;
      }
    }
    ++itSource;
  }
}

// works only for ordered columns
template <class Column, class Entry_range, typename F1, typename F2, typename F3, typename F4, typename F5>
bool _generic_add_to_column(const Entry_range& source,
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
    _generic_merge_entry_to_column(targetColumn, itSource, itTarget,
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

template <class Column, class Entry_range>
bool _add_to_column(const Entry_range& source, Column& targetColumn)
{
  return _generic_add_to_column(
      source,
      targetColumn,
      [&]([[maybe_unused]] typename Column::Entry* entryTarget) {},
      [&](typename Entry_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        if constexpr (Column::Master::Option_list::is_z2) {
          targetColumn._insert_entry(itSource->get_row_index(), itTarget);
        } else {
          targetColumn._insert_entry(itSource->get_element(), itSource->get_row_index(), itTarget);
        }
      },
      [&](typename Column::Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        if constexpr (!Column::Master::Option_list::is_z2)
          targetColumn.operators_->add_inplace(targetElement, itSource->get_element());
      },
      [&]([[maybe_unused]] typename Column::Entry* entryTarget) {},
      [&]([[maybe_unused]] typename Column::Column_support::iterator& itTarget) {}
    );
}

template <class Column, class Entry_range>
bool _multiply_target_and_add_to_column(const typename Column::Field_element& val,
                                        const Entry_range& source,
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
      [&](typename Column::Entry* entryTarget) {
        targetColumn.operators_->multiply_inplace(entryTarget->get_element(), val);
        // targetColumn.RA_opt::update_entry(*itTarget) produces an internal compiler error
        // even though it works in _generic_add_to_column... Probably because of the lambda.
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_entry(*entryTarget);
      },
      [&](typename Entry_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        targetColumn._insert_entry(itSource->get_element(), itSource->get_row_index(), itTarget);
      },
      [&](typename Column::Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        targetColumn.operators_->multiply_and_add_inplace_front(targetElement, val, itSource->get_element());
      },
      [&]([[maybe_unused]] typename Column::Entry* entryTarget) {},
      [&](typename Column::Column_support::iterator& itTarget) {
        while (itTarget != targetColumn.column_.end()) {
          typename Column::Entry* targetEntry = _get_entry<typename Column::Entry>(itTarget);
          targetColumn.operators_->multiply_inplace(targetEntry->get_element(), val);
          if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_entry(*targetEntry);
          itTarget++;
        }
      }
    );
}

template <class Column, class Entry_range>
bool _multiply_source_and_add_to_column(const typename Column::Field_element& val,
                                        const Entry_range& source,
                                        Column& targetColumn)
{
  if (val == 0u) {
    return false;
  }

  return _generic_add_to_column(
      source,
      targetColumn,
      []([[maybe_unused]] typename Column::Entry* entryTarget) {},
      [&](typename Entry_range::const_iterator& itSource, const typename Column::Column_support::iterator& itTarget) {
        typename Column::Entry* entry =
            targetColumn._insert_entry(itSource->get_element(), itSource->get_row_index(), itTarget);
        targetColumn.operators_->multiply_inplace(entry->get_element(), val);
        if constexpr (Column::Master::Option_list::has_row_access) targetColumn.update_entry(*entry);
      },
      [&](typename Column::Field_element& targetElement, typename Entry_range::const_iterator& itSource) {
        targetColumn.operators_->multiply_and_add_inplace_back(itSource->get_element(), val, targetElement);
      },
      [&]([[maybe_unused]] typename Column::Entry* entryTarget) {},
      []([[maybe_unused]] typename Column::Column_support::iterator& itTarget) {});
}

// column has to be ordered (ie. not suited for unordered_map and heap) and contain the exact values
// (ie. not suited for vector and heap). A same column but ordered differently will have another hash value.
template <class Column>
std::size_t hash_column(const Column& column)
{
  std::size_t seed = 0;
  for (auto& entry : column) {
    seed ^= std::hash<unsigned int>()(entry.get_row_index() * static_cast<unsigned int>(entry.get_element())) +
            0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_COLUMN_UTILITIES_H
