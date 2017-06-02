/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gudhi/reader_utils.h>
#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <vector>
#include <limits>
#include <utility>

using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;

int main(int argc, char** argv) {
  std::cout << "This program compute a histogram of barcode's length. A number of bins in the histogram is a parameter "
               "of this program. \n";
  if (argc != 3) {
    std::cout << "To run this program, please provide the name of a file with persistence diagram and number of "
                 "dominant intervals you would like to get \n";
    std::cout << "The third parameter of a program is the dimension of the persistence that is to be used. If your "
                 "file contains only birth-death pairs, you can skip this parameter\n";
    return 1;
  }
  unsigned dimension = std::numeric_limits<unsigned>::max();
  int dim = -1;
  if (argc > 2) {
    dim = atoi(argv[2]);
  }
  if (dim >= 0) {
    dimension = (unsigned)dim;
  }

  Persistence_intervals p(argv[1], dimension);
  std::vector<std::pair<double, double> > dominant_intervals = p.dominant_intervals(atoi(argv[2]));
  std::vector<size_t> histogram = p.histogram_of_lengths(10);

  std::stringstream gnuplot_script;
  gnuplot_script << argv[1] << "_Gnuplot_script";
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
  std::cout << "To visualize, open gnuplot and type: load \'" << gnuplot_script.str().c_str() << "\'" << std::endl;
  out.close();
  return 0;
}
