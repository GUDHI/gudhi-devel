#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

int main() {
  // Type definitions
  using Point = std::vector<double>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

  std::vector<Point> points;
  points.push_back({1.0, 1.0});
  points.push_back({7.0, 0.0});
  points.push_back({4.0, 6.0});
  points.push_back({9.0, 6.0});
  points.push_back({0.0, 14.0});
  points.push_back({2.0, 19.0});
  points.push_back({9.0, 17.0});

  // ----------------------------------------------------------------------------
  // Init of a rips complex from points
  // ----------------------------------------------------------------------------
  double threshold = 12.0;
  Rips_complex rips_complex_from_points(points, threshold, Euclidean_distance());

  Simplex_tree stree;
  rips_complex_from_points.create_complex(stree, 1);
  // ----------------------------------------------------------------------------
  // Display information about the one skeleton rips complex
  // ----------------------------------------------------------------------------
  std::cout << "Rips complex is of dimension " << stree.dimension() <<
               " - " << stree.num_simplices() << " simplices - " <<
               stree.num_vertices() << " vertices." << std::endl;

  std::cout << "Iterator on rips complex simplices in the filtration order, with [filtration value]:" <<
               std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << stree.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  return 0;
}
