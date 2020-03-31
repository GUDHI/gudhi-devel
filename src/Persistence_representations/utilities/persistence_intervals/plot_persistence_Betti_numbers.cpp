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
#include <vector>
#include <limits>
#include <utility>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  if ((argc != 3) && (argc != 2)) {
    std::clog << "This program creates a gnuplot script of Betti numbers from a single persistence diagram file"
              << "(*.pers).\n"
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

  Persistence_intervals p(argv[1], dimension);
  std::vector<std::pair<double, size_t> > pbns = p.compute_persistent_betti_numbers();

  // set up the ranges so that we see the image well.
  double xRangeBegin = pbns[0].first;
  double xRangeEnd = pbns[pbns.size() - 1].first;
  double yRangeBegin = 0;
  double yRangeEnd = 0;
  for (size_t i = 0; i != pbns.size(); ++i) {
    if (pbns[i].second > yRangeEnd) yRangeEnd = pbns[i].second;
  }
  xRangeBegin -= (xRangeEnd - xRangeBegin) / 100.0;
  xRangeEnd += (xRangeEnd - xRangeBegin) / 100.0;
  yRangeEnd += yRangeEnd / 100;

  std::stringstream gnuplot_script;
  gnuplot_script << argv[1] << "_GnuplotScript";
  std::ofstream out;
  out.open(gnuplot_script.str().c_str());

  out << "set xrange [" << xRangeBegin << " : " << xRangeEnd << "]" << std::endl;
  out << "set yrange [" << yRangeBegin << " : " << yRangeEnd << "]" << std::endl;
  out << "plot '-' using 1:2 notitle with lp " << std::endl;
  double previous_y = 0;
  for (size_t i = 0; i != pbns.size(); ++i) {
    out << pbns[i].first << " " << previous_y << std::endl;
    out << pbns[i].first << " " << pbns[i].second << std::endl;
    previous_y = pbns[i].second;
  }
  out << std::endl;
  out.close();

  std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
            << gnuplot_script.str().c_str() << "\'\"" << std::endl;

  return 0;
}
