/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <limits>
#include <vector>
#include <utility>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  if ((argc != 3) && (argc != 2)) {
    std::cout << "This program creates a gnuplot script from a single persistence diagram file (*.pers).\n"
              << "To run this program, please provide the name of a file with persistence diagram.\n"
              << "The second optional parameter of a program is the dimension of the persistence that is to be used. "
              << "If your file contains only birth-death pairs, you can skip this parameter.\n";
    return 1;
  }
  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = -1;
  if (argc == 3) {
    dim = atoi(argv[2]);
  }
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }
  std::vector<std::pair<double, double> > intervals =
      Gudhi::Persistence_representations::read_persistence_intervals_in_one_dimension_from_file(argv[1], dimension);
  Persistence_intervals b(intervals);
  b.plot(argv[1]);
  return 0;
}
