#include <CGAL/Epick_d.h>
#include <gudhi/Alpha_complex.h>

#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Epick_d< CGAL::Dimension_tag<2> > Kernel;
typedef Kernel::Point_d Point;
typedef std::vector<Point> Vector_of_points;

int main(int argc, char **argv) {
  double alpha_square_max_value = 32.0;

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
  Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_points(points, alpha_square_max_value);

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
