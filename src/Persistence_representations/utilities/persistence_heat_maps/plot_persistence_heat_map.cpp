/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_heat_maps.h>

#include <iostream>
#include <sstream>

using constant_scaling_function = Gudhi::Persistence_representations::constant_scaling_function;
using Persistence_heat_maps = Gudhi::Persistence_representations::Persistence_heat_maps<constant_scaling_function>;

int main(int argc, char** argv) {
  std::cout << "This program creates a gnuplot script from a persistence heat maps stored in a file (the file needs "
            << "to be created beforehand). Please call the code with the name of a single heat maps file \n";
  if (argc != 2) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }
  Persistence_heat_maps l;
  l.load_from_file(argv[1]);
  l.plot(argv[1]);
  return 0;
}
