/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

//for persistence algorithm
#include "gudhi/reader_utils.h"
#include "gudhi/Bitmap_cubical_complex.h"
#include "gudhi/Persistent_cohomology.h"

#include <boost/program_options.hpp>

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

//standard stuff
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
  cout << "This program computes persistent homology, by using Bitmap_cubical_complex class, of cubical complexes provided in text files in Perseus style (the only numbed in \
the first line is a dimension D of a cubical complex. In the lines I between 2 and D+1 there are numbers of top dimensional cells in the direction I. Let N denote product \
of the numbers in the lines between 2 and D. In the lines D+2 to D+2+N there are filtrations of top dimensional cells. We assume that the cells are in the \
lexicographical order. See CubicalOneSphere.txt or CubicalTwoSphere.txt for example." << endl;

  int p = 2;
  double min_persistence = 0;

  if (argc != 2) {
    cout << "Wrong number of parameters. Please provide the name of a file with a Perseus style cubical complex at the input. The program will now terminate.\n";
    return 1;
  }

  Bitmap_cubical_complex<double> b(argv[1]);


  // Compute the persistence diagram of the complex
  persistent_cohomology::Persistent_cohomology< Bitmap_cubical_complex<double>, Field_Zp > pcoh(b);
  pcoh.init_coefficients(p); //initilizes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);


  stringstream ss;
  ss << argv[1] << "_persistence";
  std::ofstream out((char*) ss.str().c_str());
  pcoh.output_diagram(out);
  out.close();

  return 0;
}
