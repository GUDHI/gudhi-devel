/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University UK
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
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Hasse_diagram_cell.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main(int argc, char** argv) 
{ 
  std::cout  << "This program computes (persistent) homology of a chain complex given in a file. The format of the file is the following:" << std::endl;
  std::cout << "In the first line, the number of all cells is given." << std::endl;
  std::cout << "The cells are assumed to be enumerated from 0 to the number of cells." << std::endl;
  std::cout << "Each cell is given in the following format:" << std::endl;
  std::cout << "Id of a cell followed by its dimension and optionally a filtration (in the first line)" << std::endl;
  std::cout << "Sequence of ids of boundary elements followed by the incidence coeficient between given cell and the boundary element (all of them in the second line)" << std::endl << std::endl;
  std::cout << "The input parameters of the program are: " << std::endl;
  std::cout << "(1) Name of the file with complex, " << std::endl;
  std::cout << "(2) Optional -- name of the output file, " << std::endl;
  std::cout << "(3) Optional -- prime number p such that (persistent) homology over Zp will be computed. " << std::endl;
  
  if (argc < 2) 
  {
    std::cerr << "Wrong number of parameters. The program will now terminate.\n";
    return 1;        
  }
  
  const char* filename = argv[1];
  const char* output_file = "output";
  unsigned field_characteristic = 11;
  if ( argc > 2 )
  {
	  output_file = argv[2];
  }
  if ( argc > 3 )
  {
	  field_characteristic = (unsigned)(atoi(argv[3]));
  }
  
  typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;	
  typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> Hasse_diag;
    
  Hasse_diag hd( filename ); 
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diag, Field_Zp> Persistent_cohomology;

  Persistent_cohomology pcoh(hd,true);  
  double min_persistence = 0;
  
  pcoh.init_coefficients(field_characteristic);    
  pcoh.compute_persistent_cohomology(min_persistence);
  
  std::ofstream out( output_file );
  pcoh.output_diagram(out);
  out.close();

  std::cout << "Result in file: " << output_file << "\n";
  
  return 0;
}

