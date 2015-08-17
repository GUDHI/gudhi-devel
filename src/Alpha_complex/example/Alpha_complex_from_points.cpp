#include <stdio.h>
#include <stdlib.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <string>
#include <vector>

// to construct a Delaunay_triangulation from a OFF file
#include "gudhi/Delaunay_triangulation_off_io.h"
#include "gudhi/Alpha_complex.h"

typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
typedef Kernel::Point_d Point;
typedef std::vector<Point> Vector_of_points;

int main(int argc, char **argv) {
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;

  std::vector<double> coords = { 0.0, 0.0, 0.0, 1.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 0.0, 1.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 1.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 1.0, 0.0, 0.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  double max_alpha_square_value = 1e10;
  Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_points(3, points.size(), points.begin(), points.end(),
                                                                       max_alpha_square_value);

  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  std::cout << "Alpha complex is of dimension " << alpha_complex_from_points.dimension() <<
      " - " << alpha_complex_from_points.num_simplices() << " simplices - " <<
      alpha_complex_from_points.num_vertices() << " vertices." << std::endl;

  std::cout << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : alpha_complex_from_points.filtration_simplex_range()) {
    std::cout << "   ( ";
    for (auto vertex : alpha_complex_from_points.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> " << "[" << alpha_complex_from_points.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }
  return 0;
}
