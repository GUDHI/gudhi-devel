/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file matrix_utilities.h
 * @author Hannah Schreiber
 * @brief Contains @ref Gudhi::persistence_matrix::Persistence_interval class.
 */

#ifndef PM_MATRIX_UTILITIES_INCLUDED
#define PM_MATRIX_UTILITIES_INCLUDED

#include <string>
#include <cstdint>    // std::uint8_t
#include <fstream>    // std::ifstream
#include <algorithm>  // std::sort
#include <utility>

#include <gudhi/Matrix.h>

namespace Gudhi {
namespace persistence_matrix {

/**
 * @ingroup persistence_matrix
 *
 * @brief Default options for @ref compute_dynamic_persistence_from_file.
 */
struct Dynamic_default_options : Default_options<> {
  static const bool is_of_boundary_type = true;                /**< RU matrix. */
  static const bool has_removable_columns = true;              /**< Enables removal. */
  static const bool has_column_pairings = true;                /**< Enables barcode output. */
  static const bool has_vine_update = true;                    /**< Enables insertion and removal at any position. */
  static const bool can_retrieve_representative_cycles = true; /**< Enables rep cycle output. */
};

/**
 * @ingroup persistence_matrix
 *
 * @private
 */
enum class Operation_types : std::uint8_t { INSERTION, REMOVAL, SWAP, PAIRS, CYCLES, COMMENT };

/**
 * @ingroup persistence_matrix
 *
 * @private
 */
template <typename Index>
Operation_types read_operation_(std::string& line, Index& position, std::vector<Index>& boundary)
{
  Operation_types type;
  boundary.clear();
  Index idx;

  std::size_t current = line.find_first_not_of(' ', 0);
  if (current == std::string::npos) return Operation_types::COMMENT;

  switch (line[current]) {
    case 'i':
      type = Operation_types::INSERTION;
      break;
    case 'r':
      type = Operation_types::REMOVAL;
      break;
    case 's':
      type = Operation_types::SWAP;
      break;
    case 'p':
      return Operation_types::PAIRS;
    case 'c':
      return Operation_types::CYCLES;
    case '#':
      return Operation_types::COMMENT;
    default:
      std::clog << "(1) Syntax error in file." << '\n';
      exit(0);
  }

  // skips whole key word, e.g. allows both "i" and "insert" as key word
  current = line.find_first_of(' ', current);
  if (current == std::string::npos) {
    std::clog << "(2) Syntax error in file." << '\n';
    exit(0);
  }
  // searches for next symbol
  current = line.find_first_not_of(' ', current + 1);
  if (current == std::string::npos) {
    std::clog << "(3) Syntax error in file." << '\n';
    exit(0);
  }
  std::size_t next = line.find_first_of(' ', current);
  position = std::stoi(line.substr(current, next - current));

  if (type != Operation_types::INSERTION) return type;

  current = line.find_first_not_of(' ', next);
  while (current != std::string::npos) {
    next = line.find_first_of(' ', current);
    idx = std::stoi(line.substr(current, next - current));
    boundary.push_back(idx);
    current = line.find_first_not_of(' ', next);
  }
  std::sort(boundary.begin(), boundary.end());

  return type;
}

/**
 * @ingroup persistence_matrix
 *
 * @brief Takes a sequence of operation from a file and applies them:
 * - insertion: inserts a cell at given position in the filtration
 * - removal: removes the cell at given position from the filtration
 * - swap: swaps the position of two adjacent cells in the filtration
 * - barcode: outputs the persistence pairs of the filtration in its current state. Note that the filtration is
 * discretized, that is, the filtration values correspond to the relative positions of the cells in the filtration.
 * - cycles: outputs the representative cycles (alive or not) of the filtration in its current state. When using the 
 * default option (see description of template parameter), a cycle is represented by the positions of the cells
 * composing it. So the bar corresponding to a cycle is the one whose birth time is equal to the highest value in
 * the cycle.
 *
 * File format:
 * - Any operation has to be on a separate line.
 * - Any empty line or line starting with `#` is ignored.
 * - insertion: `i p b0 b1 b2 ...`
 * where `p` is the position to insert to and the sequence of `b*` are the **current** positions of the cells in the
 * boundary of the cell to insert.
 * - removal: `r p`
 * where `p` is the position of the cell to remove.
 * - swap: `s p`
 * where `p` is the position of the cell which will be swapped with the cell at position `p + 1`.
 * - barcode: `p`
 * - cycles: `c`
 *
 * Any operation letter can be replaced by a word starting with this letter. E.g., `i` can be replaced by `insert`.
 * 
 * @tparam PersistenceMatrixOptions Default value: @ref Dynamic_default_options. Structure respecting the
 * @ref PersistenceMatrixOptions concept, with the following additional restrictions:
 * - `has_removable_columns` has to be `true`
 * - `has_column_pairings` has to be `true`
 * - `has_vine_update` has to be `true`
 * - `can_retrieve_representative_cycles` has to be `true`
 * - `column_indexation_type` can **not** be set to @ref Column_indexation_types::IDENTIFIER
 * - If `is_of_boundary_type` is set to `false`, `has_map_column_container` has to be set to `true`
 * - If `is_of_boundary_type` is set to `false`, `column_indexation_type` has to be set to
 * @ref Column_indexation_types::POSITION
 * - Note that the cycle format depends on the value of`is_of_boundary_type`. If `true`, it is represented by the
 * positions in the current filtration of the cells composing it. If `false`, it is represented by the cell number
 * of the cells composing it: e.g. if a cell was the \f$ n^{th} \f$ cell to be inserted in the sequence of operations,
 * then its cell number is \f$ n \f$.
 * @param filePath Path to the file to read. If the file cannot be found, a message is output in `std::clog`, but the
 * method will not throw.
 * @param output_pairs Method determining how to output the persistence pairs. Takes arguments in this order:
 * @ref Matrix::Index (number of modifications so far), const @ref Matrix::Barcode & (container of persistence pairs).
 * No returning value expected.
 * @param output_cycles  Method determining how to output the representative cycles. Takes arguments in this order:
 * @ref Matrix::Index (number of modifications so far), const std::vector< @ref Matrix::Cycle >& (cycle container),
 * const @ref Matrix & (current state of the reduced boundary matrix).
 * No returning value expected.
 */
template <class PersistenceMatrixOptions = Dynamic_default_options, class F1, class F2>
inline void compute_dynamic_persistence_from_file(const std::string& filePath, F1&& output_pairs, F2&& output_cycles)
{
  using Index = typename PersistenceMatrixOptions::Index;
  std::string line;
  std::ifstream file(filePath);

  if (file.is_open()) {
    Index step = 0;
    std::vector<Index> data;
    Index position;
    Operation_types type;
    Matrix<PersistenceMatrixOptions> matrix;

    while (getline(file, line, '\n')) {
      type = read_operation_(line, position, data);
      switch (type) {
        case Operation_types::INSERTION:
          matrix.insert_maximal_cell(position, data);
          ++step;
          break;
        case Operation_types::REMOVAL:
          matrix.remove_maximal_cell(position);
          ++step;
          break;
        case Operation_types::SWAP:
          matrix.vine_swap(position);
          ++step;
          break;
        case Operation_types::PAIRS:
          std::forward<F1>(output_pairs)(step, matrix.get_current_barcode());
          break;
        case Operation_types::CYCLES:
          matrix.update_all_representative_cycles();
          std::forward<F2>(output_cycles)(step, matrix.get_all_representative_cycles(), matrix);
          break;
        default:
          // type == COMMENT
          break;
      }
    }

    file.close();
  } else {
    // Exception instead?
    std::clog << "Unable to open input file." << '\n';
    file.setstate(std::ios::failbit);
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_MATRIX_UTILITIES_INCLUDED
