#include <gudhi/Flag_complex_edge_collapser.h>

#include <iostream>
#include <vector>

int main() {
  // Type definitions
  using Filtration_value = float;
  using Vertex_handle = short;
  using Flag_complex_edge_collapser = Gudhi::collapse::Flag_complex_edge_collapser<Vertex_handle, Filtration_value>;
  using Filtered_edge = Flag_complex_edge_collapser::Filtered_edge;
  using Filtered_edge_list = std::vector<Filtered_edge>;

  //  1   2
  //  o---o
  //  |\ /|
  //  | x |
  //  |/ \|
  //  o---o
  //  0   3
  Filtered_edge_list graph = {{0, 1, 1.},
                              {1, 2, 1.},
                              {2, 3, 1.},
                              {3, 0, 1.},
                              {0, 2, 2.},
                              {1, 3, 2.}};

  Flag_complex_edge_collapser edge_collapser(graph.begin(), graph.end());

  Filtered_edge_list remaining_edges;
  // Retrieve collapse edges from the output iterator
  edge_collapser.process_edges(
    [&remaining_edges](Vertex_handle u, Vertex_handle v, Filtration_value filtration) {
        remaining_edges.emplace_back(Filtered_edge(u, v, filtration));
      });

  for (Filtered_edge filtered_edge_from_collapse : remaining_edges) {
    std::cout << "fn[" << std::get<0>(filtered_edge_from_collapse) << ", " << std::get<1>(filtered_edge_from_collapse)
              << "] = " << std::get<2>(filtered_edge_from_collapse) << std::endl;
  }

  return 0;
}
