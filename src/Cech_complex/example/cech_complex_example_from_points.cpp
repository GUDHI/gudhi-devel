#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version

#include <iostream>
#include <string>
#include <vector>
#include <array>

int main() {
  // Type definitions
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using FT = typename Kernel::FT;
  using Point = typename Kernel::Point_d;
  using Point_cloud = std::vector<Point>;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Kernel, Simplex_tree>;

  Point_cloud points;

  std::vector<FT> point0({1., 0.});
  points.emplace_back(point0.begin(), point0.end());
  std::vector<FT> point1({0., 1.});
  points.emplace_back(point1.begin(), point1.end());
  std::vector<FT> point2({2., 1.});
  points.emplace_back(point2.begin(), point2.end());
  std::vector<FT> point3({3., 2.});
  points.emplace_back(point3.begin(), point3.end());
  std::vector<FT> point4({0., 3.});
  points.emplace_back(point4.begin(), point4.end());
  std::vector<FT> point5({3. + std::sqrt(3.), 3.});
  points.emplace_back(point5.begin(), point5.end());
  std::vector<FT> point6({1., 4.});
  points.emplace_back(point6.begin(), point6.end());
  std::vector<FT> point7({3., 4.});
  points.emplace_back(point7.begin(), point7.end());
  std::vector<FT> point8({2., 4. + std::sqrt(3.)});
  points.emplace_back(point8.begin(), point8.end());
  std::vector<FT> point9({0., 4.});
  points.emplace_back(point9.begin(), point9.end());
  std::vector<FT> point10({-0.5, 2.});
  points.emplace_back(point10.begin(), point10.end());

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
  std::clog << "Cech complex is of dimension " << stree.dimension() << " - " << stree.num_simplices() << " simplices - "
            << stree.num_vertices() << " vertices." << std::endl;

  std::clog << "Iterator on Cech complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> "
              << "[" << stree.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }
  return 0;
}
