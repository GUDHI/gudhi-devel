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
#define BOOST_PARAMETER_MAX_ARITY 12


#include <iostream>
#include <fstream>
#include <ctime>
#include <utility> 

#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Landmark_choice_by_random_point.h>
#include <gudhi/reader_utils.h>

#include "generators.h"

using namespace Gudhi;

typedef std::vector< Vertex_handle > typeVectorVertex;
//typedef std::vector< std::vector <double> > Point_Vector;

typedef Witness_complex< Simplex_tree<> > WitnessComplex;

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

/** Write a gnuplot readable file.
 *  Data range is a random access range of pairs (arg, value)
 */
template < typename Data_range >
void write_data( Data_range & data, std::string filename )
{
  std::ofstream ofs(filename, std::ofstream::out);
  for (auto entry: data)
    ofs << entry.first << ", " << entry.second << "\n";
  ofs.close();
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

  // Construct the Simplex Tree
  Simplex_tree<> simplex_tree;

  std::vector< std::pair<int, double> > l_time;      
  
  // Read the point file
  for (int nbP = 500; nbP < 10000; nbP += 500)
    {
      Point_Vector point_vector;
      generate_points_sphere(point_vector, nbP, 4);
      std::cout << "Successfully generated " << point_vector.size() << " points.\n";
      std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";
            
      // Choose landmarks
      start = clock();
      std::vector<std::vector< int > > knn;
      Landmark_choice_by_random_point(point_vector, nbL, knn);
      
      // Compute witness complex
      WitnessComplex(knn, simplex_tree, nbL, point_vector[0].size());
      end = clock();
      double time = (double)(end-start)/CLOCKS_PER_SEC;
      std::cout << "Witness complex for " << nbL << " landmarks took "
       << time << " s. \n";
      l_time.push_back(std::make_pair(nbP,time));
    }
  write_data(l_time, "w_time.dat");
}
