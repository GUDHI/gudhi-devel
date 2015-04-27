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

#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>

//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Witness_complex.h"
#include "gudhi/reader_utils.h"
//#include <boost/filesystem.hpp>

using namespace Gudhi;
//using namespace boost::filesystem;

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

void write_wl( std::string file_name, std::vector< std::vector <int> > & WL)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : WL)
    {
      for (auto l: w)
        ofs << l << " ";
      ofs << "\n";
    }
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
  /*
  boost::filesystem::path p;

  for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  std::string file_name   = argv[1];
  int nbL       = atoi(argv[2]);
  
  clock_t start, end;
  //Construct the Simplex Tree
  Witness_complex<> witnessComplex;
 
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  //std::cout << "Successfully read the points\n";
  witnessComplex.setNbL(nbL);
  //  witnessComplex.witness_complex_from_points(point_vector);
  std::vector<std::vector< int > > WL;
  std::set<int> L;
  start = clock();
  //witnessComplex.landmark_choice_by_furthest_points(point_vector, point_vector.size(), WL);
  witnessComplex.landmark_choice_by_random_points(point_vector, point_vector.size(), L);
  witnessComplex.nearest_landmarks(point_vector,L,WL);
  end = clock();
  std::cout << "Landmark choice took "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  // Write the WL matrix in a file
  mkdir("output", S_IRWXU);
  const size_t last_slash_idx = file_name.find_last_of("/");
  if (std::string::npos != last_slash_idx)
    {
      file_name.erase(0, last_slash_idx + 1);
    }
  std::string out_file = "output/"+file_name+"_"+argv[2]+".wl";
  //write_wl(out_file,WL);
  start = clock();
  witnessComplex.witness_complex(WL);
  //
  end = clock();
  std::cout << "Howdy world! The process took "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  /*
  char buffer[100];
  int i = sprintf(buffer,"%s_%s_result.txt",argv[1],argv[2]);
  if (i >= 0)
    {
      std::string out_file = (std::string)buffer;
      std::ofstream ofs (out_file, std::ofstream::out);
      witnessComplex.st_to_file(ofs);
      ofs.close();
    }
  */

  out_file = "output/"+file_name+"_"+argv[2]+".stree";
  std::ofstream ofs (out_file, std::ofstream::out);
  witnessComplex.st_to_file(ofs);
  ofs.close();

  out_file = "output/"+file_name+"_"+argv[2]+".badlinks";
  std::ofstream ofs2(out_file, std::ofstream::out);
  witnessComplex.write_bad_links(ofs2);
  ofs2.close();
}
