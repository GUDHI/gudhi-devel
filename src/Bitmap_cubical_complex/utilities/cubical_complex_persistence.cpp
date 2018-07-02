/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015 Inria
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
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main(int argc, char** argv) {
  std::cout
      << "This program computes persistent homology, by using bitmap_cubical_complex class, of cubical "
      << "complexes provided in text files in Perseus style (the only numbered in the first line is a dimension D of a"
      << "bitmap. In the lines I between 2 and D+1 there are numbers of top dimensional cells in the direction I. Let "
      << "N denote product of the numbers in the lines between 2 and D. In the lines D+2 to D+2+N there are "
      << "filtrations of top dimensional cells. We assume that the cells are in the lexicographical order. See "
      << "CubicalOneSphere.txt or CubicalTwoSphere.txt for example.\n"
      << std::endl;

  if (argc != 2) {
    std::cerr << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap at "
              << "the input. The program will now terminate.\n";
    return 1;
  }

  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;

  Bitmap_cubical_complex b(argv[1]);

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(b);
  int p = 11;
  double min_persistence = 0;

  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);

  std::string output_file_name(argv[1]);
  output_file_name += "_persistence";

  std::size_t last_in_path = output_file_name.find_last_of("/\\");

  if (last_in_path != std::string::npos) {
    output_file_name = output_file_name.substr(last_in_path + 1);
  }

  std::ofstream out(output_file_name.c_str());
  pcoh.output_diagram(out);
  out.close();

  std::cout << "Result in file: " << output_file_name << "\n";

  return 0;
}
