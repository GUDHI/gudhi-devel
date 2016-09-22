/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2014  INRIA Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

#include <gudhi/Alpha_complex.h>
#include <gudhi/Persistent_cohomology.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>  // for std::ofstream
#include <algorithm>  // for std::sort


using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<3> >;
using Point = Kernel::Point_d;
using Alpha_complex = Gudhi::alpha_complex::Alpha_complex<Kernel>;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology< Simplex_tree,
      Gudhi::persistent_cohomology::Field_Zp >;

std::vector<Point> random_points() {
  // Instanciate a random point generator
  CGAL::Random rng(0);

  // Generate "points_number" random points in a vector
  std::vector<Point> points;

  // Generates 1000 random 3D points on a sphere of radius 4.0
  CGAL::Random_points_on_sphere_d<Point> rand_outside(3, 4.0, rng);
  CGAL::cpp11::copy_n(rand_outside, 1000, std::back_inserter(points));
  // Generates 2000 random 3D points in a sphere of radius 3.0
  CGAL::Random_points_in_ball_d<Point> rand_inside(3, 3.0, rng);
  CGAL::cpp11::copy_n(rand_inside, 2000, std::back_inserter(points));

  return points;
}

/*
 * Compare two intervals by dimension, then by length.
 */
struct cmp_intervals_by_dim_then_length {
  explicit cmp_intervals_by_dim_then_length(Simplex_tree * sc)
      : sc_(sc) { }

  template<typename Persistent_interval>
  bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
    if (sc_->dimension(get < 0 > (p1)) == sc_->dimension(get < 0 > (p2)))
      return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
              > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
    else
      return (sc_->dimension(get < 0 > (p1)) > sc_->dimension(get < 0 > (p2)));
  }
  Simplex_tree* sc_;
};

int main(int argc, char **argv) {
  std::vector<Point> points = random_points();

  std::cout << "Points size=" << points.size() << std::endl;
  // Alpha complex persistence computation from generated points
  Alpha_complex alpha_complex_from_points(points);
  std::cout << "alpha_complex_from_points" << std::endl;

  Simplex_tree simplex;
  std::cout << "simplex" << std::endl;
  if (alpha_complex_from_points.create_complex(simplex, 0.6)) {
    std::cout << "simplex" << std::endl;
    // ----------------------------------------------------------------------------
    // Display information about the alpha complex
    // ----------------------------------------------------------------------------
    std::cout << "Simplicial complex is of dimension " << simplex.dimension() <<
        " - " << simplex.num_simplices() << " simplices - " <<
        simplex.num_vertices() << " vertices." << std::endl;
    
    // Sort the simplices in the order of the filtration
    simplex.initialize_filtration();
  
    std::cout << "Simplex_tree dim: " << simplex.dimension() << std::endl;
    
    Persistent_cohomology pcoh(simplex);
  
    // initializes the coefficient field for homology - Z/3Z
    pcoh.init_coefficients(3);
    pcoh.compute_persistent_cohomology(0.2);
  
    // Custom sort and output persistence
    cmp_intervals_by_dim_then_length cmp(&simplex);
    auto persistent_pairs = pcoh.get_persistent_pairs();
    std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
    for (auto pair : persistent_pairs) {
      std::cout << simplex.dimension(get<0>(pair)) << " "
            << simplex.filtration(get<0>(pair)) << " "
            << simplex.filtration(get<1>(pair)) << std::endl;
    }
  
    // Persistent Betti numbers
    std::cout << "The persistent Betti numbers in interval [0.40, 0.41] are : ";
    for (int dim = 0; dim < simplex.dimension(); dim++)
      std::cout << "b" << dim << " = " << pcoh.persistent_betti_number(dim, 0.40, 0.41) << " ; ";
    std::cout << std::endl;
  
    // Betti numbers
    std::vector<int> betti_numbers = pcoh.betti_numbers();
    std::cout << "The Betti numbers are : ";
    for (std::size_t i = 0; i < betti_numbers.size(); i++)
      std::cout << "b" << i << " = " << betti_numbers[i] << " ; ";
    std::cout << std::endl;
  }
  return 0;
}

