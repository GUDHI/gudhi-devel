#include <gudhi/Alpha_complex_3d.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits

// Complexity = FAST, weighted = true, periodic = false
using Weighted_alpha_complex_3d =
    Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>;
using Bare_point = Weighted_alpha_complex_3d::Bare_point_3;
using Weighted_point = Weighted_alpha_complex_3d::Weighted_point_3;

int main(int argc, char **argv) {
  // ----------------------------------------------------------------------------
  // Init of a list of points and weights from a small molecule
  // ----------------------------------------------------------------------------
  std::vector<Weighted_point> weighted_points;
  weighted_points.push_back(Weighted_point(Bare_point(1, -1, -1), 4.));
  weighted_points.push_back(Weighted_point(Bare_point(-1, 1, -1), 4.));
  weighted_points.push_back(Weighted_point(Bare_point(-1, -1, 1), 4.));
  weighted_points.push_back(Weighted_point(Bare_point(1, 1, 1), 4.));
  weighted_points.push_back(Weighted_point(Bare_point(2, 2, 2), 1.));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Weighted_alpha_complex_3d alpha_complex_from_points(weighted_points);

  Gudhi::Simplex_tree<> simplex;
  if (alpha_complex_from_points.create_complex(simplex)) {
    // ----------------------------------------------------------------------------
    // Display information about the alpha complex
    // ----------------------------------------------------------------------------
    std::clog << "Alpha complex is of dimension " << simplex.dimension() << " - " << simplex.num_simplices()
              << " simplices - " << simplex.num_vertices() << " vertices." << std::endl;

    std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
    for (auto f_simplex : simplex.filtration_simplex_range()) {
      std::clog << "   ( ";
      for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << " ";
      }
      std::clog << ") -> "
                << "[" << simplex.filtration(f_simplex) << "] ";
      std::clog << std::endl;
    }
  }
  return 0;
}
