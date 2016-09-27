/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi, small modifications by Pawel Dlotko
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

#define CGAL_HAS_THREADS

#include <gudhi/Graph_matching.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

std::vector< std::pair<double, double> > read_diagram_from_file( const char* filename )
{
    std::ifstream in;
    in.open( filename );
    std::vector< std::pair<double, double> > result;
    if ( !in.is_open() )
    {
        std::cerr << "File : " << filename << " do not exist. The program will now terminate \n";
        throw "File do not exist \n";
    }

    std::string line;
    while (!in.eof())
    {
        getline(in,line);
        if ( line.length() != 0 )
        {
            std::stringstream lineSS;
            lineSS << line;
            double beginn, endd;
            lineSS >> beginn;
            lineSS >> endd;
            result.push_back( std::make_pair( beginn , endd ) );
        }
    }
    in.close();
    return result;
}  //read_diagram_from_file

int main( int argc , char** argv )
{
  if ( argc < 3 )
  {
      std::cout << "To run this program please provide as an input two files with persistence diagrams. Each file " <<
          "should contain a birth-death pair per line. Third, optional parameter is an error bound on a bottleneck" <<
          " distance (set by default to zero). The program will now terminate \n";
  }
  std::vector< std::pair< double , double > > diag1 = read_diagram_from_file( argv[1] );
  std::vector< std::pair< double , double > > diag2 = read_diagram_from_file( argv[2] );
  double tolerance = 0;
  if ( argc == 4 )
  {
      tolerance = atof( argv[3] );
  }
  double b = Gudhi::bottleneck_distance::compute(diag1, diag2, tolerance);
  std::cout << "The distance between the diagrams is : " << b << ". The tolerace is : " << tolerance << std::endl;
}
