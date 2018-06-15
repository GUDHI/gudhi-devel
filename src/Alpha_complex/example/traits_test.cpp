#include <gudhi/Alpha_complex_3d.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " [alpha_square_max_value]\n";
  std::cerr << "       i.e.: " << progName << " 60.0\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  //if ((argc != 1) && (argc != 2)) usage(argc, (argv[0] - 1));

  using Alpha_shapes_3d = Gudhi::alpha_complex::Alpha_shapes_3d;
  std::vector<Alpha_shapes_3d::Point_3> points;
  points.push_back(Alpha_shapes_3d::Point_3(1., 2., 3.));
  points.push_back(Alpha_shapes_3d::Point_3(6., 5., 4.));

  Gudhi::alpha_complex::Alpha_complex_3d<Alpha_shapes_3d> alpha_complex(points);

  using Weighted_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_alpha_shapes_3d;
  std::vector<Weighted_alpha_shapes_3d::Point_3> w_points;
  w_points.push_back(Alpha_shapes_3d::Point_3(1., 2., 3.));
  w_points.push_back(Alpha_shapes_3d::Point_3(6., 5., 4.));

  std::vector<double> weights = {1., 2.};

  Gudhi::alpha_complex::Alpha_complex_3d<Weighted_alpha_shapes_3d> weighted_alpha_complex(points, weights);

  return 0;
}
