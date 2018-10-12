/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
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

#include <gudhi/Persistence_landscape.h>

#include <iostream>
#include <sstream>

using Persistence_landscape = Gudhi::Persistence_representations::Persistence_landscape;

int main(int argc, char** argv) {
  std::cout << "This program creates a gnuplot script from a persistence landscape stored in a file (the file needs "
            << "to be created beforehand). Please call the code with the name of a single landscape file.\n";
  if (argc != 2) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  Persistence_landscape l;
  l.load_landscape_from_file(argv[1]);
  l.plot(argv[1]);

  return 0;
}
