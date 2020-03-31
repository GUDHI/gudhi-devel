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
  std::clog << "This program computes a histogram of barcode's length. A number of bins in the histogram is a "
            << "parameter of this program. \n";
  if ((argc != 3) && (argc != 4)) {
    std::clog << "To run this program, please provide the name of a file with persistence diagram and number of "
              << "dominant intervals you would like to get. Set a negative number dominant intervals value "
              << "If your file contains only birth-death pairs.\n"
              << "The third parameter is the dimension of the persistence that is to be used. If your "
              << "file contains only birth-death pairs, you can skip this parameter\n";
    return 1;
  }

  unsigned dominant_interval_number = std::numeric_limits<unsigned>::max();
  int nbr = atoi(argv[2]);
  if (nbr >= 0) {
    dominant_interval_number = static_cast<unsigned>(nbr);
  }

  int persistence_dimension = -1;
  if (argc == 4) {
    persistence_dimension = atoi(argv[3]);
  }

  Persistence_intervals p(argv[1], dominant_interval_number);
  std::vector<std::pair<double, double> > dominant_intervals = p.dominant_intervals(persistence_dimension);
  std::vector<size_t> histogram = p.histogram_of_lengths(10);

  std::stringstream gnuplot_script;
  gnuplot_script << argv[1] << "_GnuplotScript";
  std::ofstream out;
  out.open(gnuplot_script.str().c_str());

  out << "set style data histogram" << std::endl;
  out << "set style histogram cluster gap 1" << std::endl;
  out << "set style fill solid border -1" << std::endl;
  out << "plot '-' notitle" << std::endl;
  for (size_t i = 0; i != histogram.size(); ++i) {
    out << histogram[i] << std::endl;
  }
  out << std::endl;
  out.close();

  std::clog << "To visualize, install gnuplot and type the command: gnuplot -persist -e \"load \'"
            << gnuplot_script.str().c_str() << "\'\"" << std::endl;
  return 0;
}
