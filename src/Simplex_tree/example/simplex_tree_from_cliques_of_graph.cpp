/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <ctime>
#include <string>
#include <utility>  // for std::pair

using namespace Gudhi;

typedef int Vertex_handle;
typedef double Filtration_value;
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
  // Construct the Simplex Tree
  Simplex_tree<> st;

  start = clock();
  auto g = read_graph<Graph_t, Filtration_value, Vertex_handle>(filegraph);
  // insert the graph in the simplex tree as 1-skeleton
  st.insert_graph(g);
  end = clock();
  std::clog << "Insert the 1-skeleton in the simplex tree in "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

  start = clock();
  // expand the 1-skeleton until dimension max_dim
  st.expansion(max_dim);
  end = clock();
  std::clog << "max_dim = " << max_dim << "\n";
  std::clog << "Expand the simplex tree in "
      << static_cast<double>(end - start) / CLOCKS_PER_SEC << " s. \n";

  std::clog << "Information of the Simplex Tree: " << std::endl;
  std::clog << "  Number of vertices = " << st.num_vertices() << " ";
  std::clog << "  Number of simplices = " << st.num_simplices() << std::endl;
  std::clog << std::endl << std::endl;

  std::clog << "Iterator on vertices: ";
  for (auto vertex : st.complex_vertex_range()) {
    std::clog << vertex << " ";
  }

  std::clog << std::endl;

  std::clog << std::endl << std::endl;

  std::clog << "Iterator on simplices: " << std::endl;
  for (auto simplex : st.complex_simplex_range()) {
    std::clog << "   ";
    for (auto vertex : st.simplex_vertex_range(simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;
  }

  std::clog << std::endl << std::endl;

  std::clog << "Iterator on Simplices in the filtration, with [filtration value]:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;
  }

  std::clog << std::endl << std::endl;

  std::clog << "Iterator on Simplices in the filtration, and their boundary simplices:" << std::endl;
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::clog << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << std::endl;

    for (auto b_simplex : st.boundary_simplex_range(f_simplex)) {
      std::clog << "      " << "[" << st.filtration(b_simplex) << "] ";
      for (auto vertex : st.simplex_vertex_range(b_simplex)) {
        std::clog << vertex << " ";
      }
      std::clog << std::endl;
    }
  }
  return 0;
}
