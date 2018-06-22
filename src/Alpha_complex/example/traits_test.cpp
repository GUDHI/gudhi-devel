#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits
#include <random>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " [alpha_square_max_value]\n";
  std::cerr << "       i.e.: " << progName << " 60.0\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  //if ((argc != 1) && (argc != 2)) usage(argc, (argv[0] - 1));

  std::cout << "Alpha complex 3d" << std::endl;
  using Alpha_shapes_3d = Gudhi::alpha_complex::Alpha_shapes_3d;
  std::vector<Alpha_shapes_3d::Point_3> points;
  points.push_back(Alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  points.push_back(Alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  points.push_back(Alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  points.push_back(Alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  Gudhi::alpha_complex::Alpha_complex_3d<Alpha_shapes_3d> alpha_complex(points);

  Gudhi::Simplex_tree<> stree;
  alpha_complex.create_complex(stree);

  std::cout << "Exact alpha complex 3d" << std::endl;
  using Exact_alpha_shapes_3d = Gudhi::alpha_complex::Exact_alpha_shapes_3d;
  std::vector<Exact_alpha_shapes_3d::Point_3> e_points;
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  e_points.push_back(Exact_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  Gudhi::alpha_complex::Alpha_complex_3d<Exact_alpha_shapes_3d> exact_alpha_complex(e_points);

  Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence> exact_stree;
  exact_alpha_complex.create_complex(exact_stree);

  std::cout << "Weighted alpha complex 3d" << std::endl;
  using Weighted_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_alpha_shapes_3d;
  std::vector<Weighted_alpha_shapes_3d::Point_3> w_points;
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  w_points.push_back(Weighted_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));

  std::vector<double> weights = {0.01, 0.005, 0.006, 0.01, 0.009, 0.001};

  Gudhi::alpha_complex::Alpha_complex_3d<Weighted_alpha_shapes_3d> weighted_alpha_complex(w_points, weights);

  Gudhi::Simplex_tree<> w_stree;
  weighted_alpha_complex.create_complex(w_stree);

  std::cout << "Periodic alpha complex 3d" << std::endl;
  using Periodic_alpha_shapes_3d = Gudhi::alpha_complex::Periodic_alpha_shapes_3d;
  std::vector<Periodic_alpha_shapes_3d::Point_3> p_points;

  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.1, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.6));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.8));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.0));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.4));
  p_points.push_back(Periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.6));

  Gudhi::alpha_complex::Alpha_complex_3d<Periodic_alpha_shapes_3d> periodic_alpha_complex(p_points,
                                                                                          0., 0., 0.,
                                                                                          1., 1., 1.);

  Gudhi::Simplex_tree<> p_stree;
  periodic_alpha_complex.create_complex(p_stree);

  std::cout << "Weighted periodic alpha complex 3d" << std::endl;
  using Weighted_periodic_alpha_shapes_3d = Gudhi::alpha_complex::Weighted_periodic_alpha_shapes_3d;
  std::vector<Weighted_periodic_alpha_shapes_3d::Point_3> wp_points;
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.0, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.2, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.4, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.1, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.6, 0.8, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.0, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.2, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.4, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.6));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.6, 0.8));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.0));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.2));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.4));
  wp_points.push_back(Weighted_periodic_alpha_shapes_3d::Point_3(0.8, 0.8, 0.6));

  std::vector<double> p_weights;

  std::random_device rd;
  std::mt19937 mt(rd());
  // Weights must be in range [0, <1/64]
  std::uniform_real_distribution<double> dist(0.0, 0.0156245);

  for (std::size_t i = 0; i < wp_points.size(); ++i) {
    double value = dist(mt);
    std::cout << value << std::endl;
    p_weights.push_back(value);
  }

  Gudhi::alpha_complex::Alpha_complex_3d<Weighted_periodic_alpha_shapes_3d>
      weighted_periodic_alpha_complex(wp_points, p_weights,
                                      0., 0., 0.,
                                      1., 1., 1.);
  Gudhi::Simplex_tree<> wp_stree;
  weighted_periodic_alpha_complex.create_complex(wp_stree);

  return 0;
}
