 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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
#include <ctime>
#include "gudhi/reader_utils.h"
#include "gudhi/Simplex_tree.h"

using namespace Gudhi;

int main (int argc, char * const argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
      << " path_to_file_graph max_dim \n";
    return 0;
  }
  std::string filegraph   = argv[1];
  int max_dim       = atoi(argv[2]);

  clock_t start, end;
  //Construct the Simplex Tree
  Simplex_tree<> st;

  start = clock();
  auto g = read_graph(filegraph);
  st.insert_graph (g); //insert the graph in the simplex tree as 1-skeleton
  end = clock();
  std::cout << "Insert the 1-skeleton in the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  start = clock();
  st.expansion ( max_dim ); //expand the 1-skeleton until dimension max_dim
  end = clock();
  std::cout << "max_dim = " << max_dim << "\n";
  std::cout << "Expand the simplex tree in "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";

  std::cout << "Information of the Simplex Tree: " << std::endl;
  std::cout << "  Number of vertices = " << st.num_vertices() << " ";
  std::cout << "  Number of simplices = " << st.num_simplices() << std::endl;
  std::cout << std::endl << std::endl;

  std::cout << "Iterator on vertices: ";
    for( auto vertex : st.complex_vertex_range() ) { std::cout << vertex << " "; }

  std::cout << std::endl;

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on simplices: " << std::endl;
  for( auto simplex : st.complex_simplex_range() )
  {
    std::cout << "   ";
    for( auto vertex : st.simplex_vertex_range(simplex) ) { std::cout << vertex << " "; }
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  { std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for( auto vertex : st.simplex_vertex_range(f_simplex) )
      { std::cout << vertex << " "; } std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on Simplices in the filtration, and their boundary simplices:" << std::endl;
  for( auto f_simplex : st.filtration_simplex_range() )
  {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for( auto vertex : st.simplex_vertex_range(f_simplex) )
      { std::cout << vertex << " "; } std::cout << std::endl;

    for( auto b_simplex : st.boundary_simplex_range(f_simplex) )
    {
      std::cout << "      " << "[" << st.filtration(b_simplex) << "] ";
    for( auto vertex : st.simplex_vertex_range(b_simplex) )
      { std::cout << vertex << " "; } std::cout << std::endl;
    }
  }
  return 0;
}
