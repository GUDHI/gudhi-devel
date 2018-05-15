/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017  INRIA
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
#include "gudhi/Fake_simplex_tree.h"

#include <iostream>
#include <ctime>
#include <string>
#include <utility>  // for std::pair

using Toplex_map = Gudhi::Fake_simplex_tree;
using Vertex_handle = Toplex_map::Vertex_handle;
using Filtration_value = Toplex_map::Filtration_value;

typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
                                boost::property < vertex_filtration_t, Filtration_value >,
                                boost::property < edge_filtration_t, Filtration_value > > Graph_t;

int main(int argc, char * const argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_file_graph max_dim \n";
    return 0;
  }
  std::string filegraph = argv[1];
  int max_dim = atoi(argv[2]);

  clock_t start, end;
  // Construct the Toplex Map
  Toplex_map t_map;

  start = clock();
  auto g = Gudhi::read_graph<Graph_t, Filtration_value, Vertex_handle>(filegraph);
  // insert the graph in the toplex map as 1-skeleton
  t_map.insert_graph(g);
  end = clock();
  std::cout << "Insert the 1-skeleton in the toplex map in "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

  start = clock();
  // expand the 1-skeleton until dimension max_dim
  t_map.expansion(max_dim);
  end = clock();
  std::cout << "max_dim = " << max_dim << "\n";
  std::cout << "Expand the toplex map in "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

  std::cout << "Information of the toplex map: " << std::endl;
  std::cout << "  Number of vertices = " << t_map.num_vertices() << " ";
  std::cout << "  Number of simplices = " << t_map.num_simplices() << std::endl;
  std::cout << std::endl << std::endl;

  std::cout << "Iterator on Simplices in the filtration:" << std::endl;
  for (auto f_simplex : t_map.filtration_simplex_range()) {
    std::cout << "   " << "[" << t_map.filtration(f_simplex) << "] ";
    for (auto vertex : t_map.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  std::cout << "Iterator on skeleton:" << std::endl;
  for (auto f_simplex : t_map.skeleton_simplex_range()) {
    std::cout << "   " << "[" << t_map.filtration(f_simplex) << "] ";
    for (auto vertex : t_map.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
}
