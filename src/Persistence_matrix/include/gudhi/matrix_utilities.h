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
#include <cstdint>  // std::uint8_t
#include <fstream>  // std::ifstream
#include <utility>

#include <gudhi/Matrix.h>

namespace Gudhi {
namespace persistence_matrix {

struct Dynamic_default_options : Default_options<> {
  static const bool is_of_boundary_type = true;
  static const bool has_removable_columns = true;
  static const bool has_column_pairings = true;
  static const bool has_vine_update = true;
  static const bool can_retrieve_representative_cycles = true;
};

/**
 * @private
 */
enum class Operation_types : std::uint8_t { INSERTION, REMOVAL, SWAP, PAIRS, CYCLES, COMMENT };

/**
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
      std::clog << "(1) Syntax error in file." << std::endl;
      exit(0);
  }

  // skips whole key word, e.g. allows both "i" and "insert" as key word
  current = line.find_first_of(' ', current);
  if (current == std::string::npos) {
    std::clog << "(2) Syntax error in file." << std::endl;
    exit(0);
  }
  // searches for next symbol
  current = line.find_first_not_of(' ', current + 1);
  if (current == std::string::npos) {
    std::clog << "(3) Syntax error in file." << std::endl;
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

  return type;
}

template <class F1, class F2, class PersistenceMatrixOptions = Dynamic_default_options>
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
    std::clog << "Unable to open input file." << std::endl;
    file.setstate(std::ios::failbit);
  }
}

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // PM_MATRIX_UTILITIES_INCLUDED
