#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <string>
#include <vector>
#include <array>

int main() {
  // Type definitions
  using Point_cloud = std::vector<std::array<double, 2>>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;

  Point_cloud points;
  points.push_back({1., 0.});                 // 0
  points.push_back({0., 1.});                 // 1
  points.push_back({2., 1.});                 // 2
  points.push_back({3., 2.});                 // 3
  points.push_back({0., 3.});                 // 4
  points.push_back({3. + std::sqrt(3.), 3.}); // 5
  points.push_back({1., 4.});                 // 6
  points.push_back({3., 4.});                 // 7
  points.push_back({2., 4. + std::sqrt(3.)}); // 8
  points.push_back({0., 4.});                 // 9
  points.push_back({-0.5, 2.});               // 10

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from points
  // ----------------------------------------------------------------------------
  Filtration_value max_radius = 1.;
  Cech_complex cech_complex_from_points(points, max_radius);

  Simplex_tree stree;
  cech_complex_from_points.create_complex(stree, 2);
  // ----------------------------------------------------------------------------
  // Display information about the one skeleton Cech complex
  // ----------------------------------------------------------------------------
  std::cout << "Cech complex is of dimension " << stree.dimension() <<
               " - " << stree.num_simplices() << " simplices - " <<
               stree.num_vertices() << " vertices." << std::endl;

  std::cout << "Iterator on Cech complex simplices in the filtration order, with [filtration value]:" <<
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
