/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2026 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <fstream>
#include <string>
#include <vector>

#include <gudhi/Matrix.h>
#include <gudhi/matrix_utilities.h>

using Matrix = Gudhi::persistence_matrix::Matrix<Gudhi::persistence_matrix::Dynamic_default_options>;
using Index = Matrix::Index;
using Barcode = Matrix::Barcode;
using Cycles = std::vector<Matrix::Cycle>;

void output_in_stream(const std::string& inputFilePath, std::ostream& stream)
{
  auto output_pairs = [&stream](Index step, const Barcode& barcode) {
    stream << "Barcode step " << step << ":\n";
    for (const auto& bar : barcode) stream << bar << "\n";
    stream << "\n";
  };
  auto output_cycles = [&stream](Index step, const Cycles& cycles, const Matrix& m) {
    stream << "Cycles step " << step << ":\n";
    for (const auto& cycle : cycles) {
      stream << m.get_column_dimension(cycle[0]);
      stream << "-cycle: ";
      for (auto index : cycle) {
        stream << index << ", ";
      }
      stream << "\n";
    }
    stream << "\n";
  };

  Gudhi::persistence_matrix::compute_dynamic_persistence_from_file(inputFilePath, output_pairs, output_cycles);
}

void output_in_terminal(const std::string& inputFilePath) { output_in_stream(inputFilePath, std::cout); }

void output_in_file(const std::string& inputFilePath, const std::string& outputFilePath)
{
  std::ofstream outfile(outputFilePath);

  if (!outfile.is_open()) {
    std::clog << "Unable to open input file." << std::endl;
    outfile.setstate(std::ios::failbit);
    return;
  }

  output_in_stream(inputFilePath, outfile);
  outfile.close();  // not necessary at the end of the method, but a good reflex
}

// example of input file: example/dynamic_persistence_example.txt
int main(int argc, char* const argv[])
{
  if (argc != 3) {
    if (argc < 3)
      std::clog << "Missing argument: input file name and output file name are needed." << std::endl;
    else
      std::clog << "Too many arguments: only input file name and output file name are needed." << std::endl;
    return 0;
  }

  // two examples of outputs
  // both examples output in streams, but one could also decide to store everything or just parts in a
  // vector for example: `output_pairs` and `output_cycles` can do whatever as long as the parameters remain the same.
  output_in_terminal(argv[1]);
  output_in_file(argv[1], argv[2]);
  // note that both output in terminal and in file could have been done at the same time, the separation is just
  // for illustration purposes.

  return 0;
}