#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epeck_d.h>  // For EXACT or SAFE version

#include <iostream>
#include <string>
#include <vector>
#include <array>

int main() {
  // Type definitions
//   using Point_cloud = std::vector<std::array<double, 2>>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Kernel = CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>;
  using FT = typename Kernel::FT;
  using Point = typename Kernel::Point_d;
  using Point_cloud = std::vector<Point>;
  using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud, Kernel>;

  Point_cloud points;
//   points.push_back({1., 0.});                  // 0
//   points.push_back({0., 1.});                  // 1
//   points.push_back({2., 1.});                  // 2
//   points.push_back({3., 2.});                  // 3
//   points.push_back({0., 3.});                  // 4
//   points.push_back({3. + std::sqrt(3.), 3.});  // 5

//   std::vector<FT> point({0.0, 0.0, 0.0, 0.0});
//   points.emplace_back(point.begin(), point.end());

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

//   points.emplace_back(Point(std::vector<FT>({1., 0.})));
//   points.emplace_back(Point(std::vector<FT>({0., 1.})));
//   points.emplace_back(Point(std::vector<FT>({2., 1.})));
//   points.emplace_back(Point(std::vector<FT>({3., 2.})));
//   points.emplace_back(Point(std::vector<FT>({0., 3.})));
//   points.emplace_back(Point(std::vector<FT>({3. + std::sqrt(3.), 3.})));


//   points.push_back(Point(1.0, 0.0));
//   points.push_back(Point(0.0, 1.0));
//   points.push_back(Point(2.0, 1.0));
//   points.push_back(Point(3.0, 2.0));
//   points.push_back(Point(0.0, 3.0));
//   points.push_back(Point(3.0 + std::sqrt(3.0), 3.0));


//   points.push_back({1., 4.});                  // 6
//   points.push_back({3., 4.});                  // 7
//   points.push_back({2., 4. + std::sqrt(3.)});  // 8
//   points.push_back({0., 4.});                  // 9
//   points.push_back({-0.5, 2.});                // 10

  // ----------------------------------------------------------------------------
  // Init of a Cech complex from points
  // ----------------------------------------------------------------------------
  Filtration_value max_radius = 10.;
  std::clog << "Hind: Just before the Cech constructor" << std::endl;
  Cech_complex cech_complex_from_points(points, max_radius);
  std::clog << "Hind: Just after the Cech constructor" << std::endl;

  Simplex_tree stree;
  cech_complex_from_points.create_complex(stree, 3);
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
