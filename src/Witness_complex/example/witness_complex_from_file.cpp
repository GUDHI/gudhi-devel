/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#include <iostream>
#include <fstream>
#include <ctime>
//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Witness_complex.h"


using namespace Gudhi;

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::vector< std::vector <double> > Point_Vector;
//typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
//typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;


/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust ( std::string file_name , std::vector< std::vector< double > > & points)
{  
  std::ifstream in_file (file_name.c_str(),std::ios::in);
  if(!in_file.is_open())
    {
      std::cerr << "Unable to open file " << file_name << std::endl;
      return;
    }
  std::string line;
  double x;
  while( getline ( in_file , line ) )
    {
      std::vector< double > point;
      std::istringstream iss( line );
      while(iss >> x) { point.push_back(x); }
      if (point.size() != 1)
        points.push_back(point);
    }
  in_file.close();
}

int main (int argc, char * const argv[])
{
    if (argc != 3)
      {
          std::cerr << "Usage: " << argv[0]
    << " path_to_point_file nbL \n";
    return 0;
  }
  std::string file_name   = argv[1];
  int nbL       = atoi(argv[2]);
  
  clock_t start, end;
  //Construct the Simplex Tree
  Witness_complex<> witnessComplex;
 
  std::cout << "Let the carnage begin!\n";
  start = clock();
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  witnessComplex.setNbL(nbL);
  witnessComplex.witness_complex_from_points(point_vector);
  end = clock();
  std::cout << "Howdy world! The process took "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

}
