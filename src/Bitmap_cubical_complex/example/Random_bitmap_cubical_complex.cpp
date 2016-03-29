/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Saclay (France)
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


// for persistence algorithm
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>

int main(int argc, char** argv) {
  srand(time(0));

  std::cout << "This program computes persistent homology, by using bitmap_cubical_complex class, of cubical " <<
      "complexes. The first parameter of the program is the dimension D of the bitmap. The next D parameters are " <<
      "number of top dimensional cubes in each dimension of the bitmap. The program will create random cubical " <<
      "complex of that sizes and compute persistent homology of it." << std::endl;

  int p = 2;
  double min_persistence = 0;

  if (argc < 3) {
    std::cerr << "Wrong number of parameters, the program will now terminate\n";
    return 1;
  }

  size_t dimensionOfBitmap = (size_t) atoi(argv[1]);
  std::vector< unsigned > sizes;
  size_t multipliers = 1;
  for (size_t dim = 0; dim != dimensionOfBitmap; ++dim) {
    unsigned sizeInThisDimension = (unsigned) atoi(argv[2 + dim]);
    sizes.push_back(sizeInThisDimension);
    multipliers *= sizeInThisDimension;
  }

  std::vector< double > data;
  for (size_t i = 0; i != multipliers; ++i) {
    data.push_back(rand() / static_cast<double>(RAND_MAX));
  }

  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::Cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
  Bitmap_cubical_complex b(sizes, data);

  // Compute the persistence diagram of the complex
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  Persistent_cohomology pcoh(b);
  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);

  std::stringstream ss;
  ss << "randomComplex_persistence";
  std::ofstream out(ss.str().c_str());
  pcoh.output_diagram(out);
  out.close();

  return 0;
}

