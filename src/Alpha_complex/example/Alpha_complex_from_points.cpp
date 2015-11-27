#include <stdlib.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>
#include <vector>

#include <gudhi/Alpha_complex.h>

typedef CGAL::Epick_d< CGAL::Dimension_tag<2> > Kernel;
typedef Kernel::Point_d Point;
typedef std::vector<Point> Vector_of_points;

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " alpha_square_max_value" << std::endl;
  std::cerr << "       i.e.: " << progName << " 32.0" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }

  double alpha_square_max_value = atof(argv[1]);
  
  // ----------------------------------------------------------------------------
  // Init of a list of points
  // ----------------------------------------------------------------------------
  Vector_of_points points;

  std::vector<double> coords = { 1.0, 1.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 7.0, 0.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 4.0, 6.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 9.0, 6.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 0.0, 14.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 2.0, 19.0 };
  points.push_back(Point(coords.begin(), coords.end()));
  coords = { 9.0, 17.0 };
  points.push_back(Point(coords.begin(), coords.end()));

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
