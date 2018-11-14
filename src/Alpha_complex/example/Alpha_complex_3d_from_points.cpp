#include <gudhi/Alpha_complex_3d.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits

using Alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>;
using Point = Alpha_complex_3d::Point_3;
using Vector_of_points = std::vector<Point>;

int main(int argc, char **argv) {
  if (argc != 1) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct\n";
    std::cerr << "Usage: " << (argv[0] - 1) << " \n";
    exit(-1);  // ----- >>
  }

  // ----------------------------------------------------------------------------
  // Init of a list of points from a small molecule
  // ----------------------------------------------------------------------------
  Vector_of_points points;
  points.push_back(Point(1, -1, -1));
  points.push_back(Point(-1, 1, -1));
  points.push_back(Point(-1, -1, 1));
  points.push_back(Point(1, 1, 1));
  points.push_back(Point(2, 2, 2));

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from the list of points
  // ----------------------------------------------------------------------------
  Alpha_complex_3d alpha_complex_from_points(points);

  Gudhi::Simplex_tree<> simplex;
  if (alpha_complex_from_points.create_complex(simplex)) {
    // ----------------------------------------------------------------------------
    // Display information about the alpha complex
    // ----------------------------------------------------------------------------
    std::cout << "Alpha complex is of dimension " << simplex.dimension() << " - " << simplex.num_simplices()
              << " simplices - " << simplex.num_vertices() << " vertices." << std::endl;

    std::cout << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
    for (auto f_simplex : simplex.filtration_simplex_range()) {
      std::cout << "   ( ";
      for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
        std::cout << vertex << " ";
      }
      std::cout << ") -> "
                << "[" << simplex.filtration(f_simplex) << "] ";
      std::cout << std::endl;
    }
  }
  return 0;
}
