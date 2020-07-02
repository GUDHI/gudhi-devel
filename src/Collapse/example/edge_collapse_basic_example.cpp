#include <gudhi/Flag_complex_edge_collapser.h>

#include <iostream>
#include <vector>
#include <tuple>

int main() {
  // Type definitions
  using Filtration_value = float;
  using Vertex_handle = short;
  using Filtered_edge = std::tuple<Vertex_handle, Vertex_handle, Filtration_value>;
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

  auto remaining_edges = Gudhi::collapse::flag_complex_collapse_edges(graph);

  for (auto filtered_edge_from_collapse : remaining_edges) {
    std::cout << "fn[" << std::get<0>(filtered_edge_from_collapse) << ", " << std::get<1>(filtered_edge_from_collapse)
              << "] = " << std::get<2>(filtered_edge_from_collapse) << std::endl;
  }

  return 0;
}
