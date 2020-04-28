#include <gudhi/Alpha_complex.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epeck_d.h>

#include <iostream>
#include <vector>

// Explicit dimension 2 Epeck_d kernel
using Kernel = CGAL::Epeck_d< CGAL::Dimension_tag<2> >;
using Point = Kernel::Point_d;
using Vector_of_points = std::vector<Point>;

int main() {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;
  points.push_back(Point(1.0, 1.0));
  points.push_back(Point(7.0, 0.0));
  points.push_back(Point(4.0, 6.0));
  points.push_back(Point(9.0, 6.0));
  points.push_back(Point(0.0, 14.0));
  points.push_back(Point(2.0, 19.0));
  points.push_back(Point(9.0, 17.0));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_points(points);

  Gudhi::Simplex_tree<> simplex;
  if (alpha_complex_from_points.create_complex(simplex)) {
    // ----------------------------------------------------------------------------
    // Display information about the alpha complex
    // ----------------------------------------------------------------------------
    std::clog << "Alpha complex is of dimension " << simplex.dimension() <<
        " - " << simplex.num_simplices() << " simplices - " <<
        simplex.num_vertices() << " vertices." << std::endl;

    std::clog << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
    for (auto f_simplex : simplex.filtration_simplex_range()) {
      std::clog << "   ( ";
      for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
        std::clog << vertex << " ";
      }
      std::clog << ") -> " << "[" << simplex.filtration(f_simplex) << "] ";
      std::clog << std::endl;
    }
  }
  return 0;
}
