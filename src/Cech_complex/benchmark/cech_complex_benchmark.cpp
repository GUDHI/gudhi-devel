/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
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

#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Clock.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Miniball.hpp>

#include "boost/filesystem.hpp"  // includes all needed Boost.Filesystem declarations

#include <string>
#include <vector>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<Filtration_value>;
using Point_cloud = std::vector<Point>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Proximity_graph = Gudhi::Proximity_graph<Simplex_tree>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Cech_complex = Gudhi::cech_complex::Cech_complex<Simplex_tree, Point_cloud>;

class Minimal_enclosing_ball_radius {
 public:
  // boost::range_value is not SFINAE-friendly so we cannot use it in the return type
  template <typename Point>
  typename std::iterator_traits<typename boost::range_iterator<Point>::type>::value_type operator()(
      const Point& p1, const Point& p2) const {
    // Type def
    using Point_cloud = std::vector<Point>;
    using Point_iterator = typename Point_cloud::const_iterator;
    using Coordinate_iterator = typename Point::const_iterator;
    using Min_sphere =
        typename Gudhi::Miniball::Miniball<Gudhi::Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;

    Point_cloud point_cloud;
    point_cloud.push_back(p1);
    point_cloud.push_back(p2);

    GUDHI_CHECK((p1.end() - p1.begin()) != (p2.end() - p2.begin()), "inconsistent point dimensions");
    Min_sphere min_sphere(p1.end() - p1.begin(), point_cloud.begin(), point_cloud.end());

    return std::sqrt(min_sphere.squared_radius());
  }
};

int main(int argc, char* argv[]) {
  std::string off_file_points = "tore3D_1307.off";
  Filtration_value threshold = 1e20;

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  Gudhi::Clock euclidean_clock("Gudhi::Euclidean_distance");
  // Compute the proximity graph of the points
  Proximity_graph euclidean_prox_graph = Gudhi::compute_proximity_graph<Simplex_tree>(
      off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());

  std::cout << euclidean_clock << std::endl;

  Gudhi::Clock miniball_clock("Minimal_enclosing_ball_radius");
  // Compute the proximity graph of the points
  Proximity_graph miniball_prox_graph = Gudhi::compute_proximity_graph<Simplex_tree>(
      off_reader.get_point_cloud(), threshold, Minimal_enclosing_ball_radius());
  std::cout << miniball_clock << std::endl;

  Gudhi::Clock common_miniball_clock("Gudhi::Minimal_enclosing_ball_radius()");
  // Compute the proximity graph of the points
  Proximity_graph common_miniball_prox_graph = Gudhi::compute_proximity_graph<Simplex_tree>(
      off_reader.get_point_cloud(), threshold, Gudhi::Minimal_enclosing_ball_radius());
  std::cout << common_miniball_clock << std::endl;

  boost::filesystem::path full_path(boost::filesystem::current_path());
  std::cout << "Current path is : " << full_path << std::endl;

  std::cout << "File name;Radius;Rips time;Cech time; Ratio Rips/Cech time;Rips nb simplices;Cech nb simplices;"
            << std::endl;
  boost::filesystem::directory_iterator end_itr;  // default construction yields past-the-end
  for (boost::filesystem::directory_iterator itr(boost::filesystem::current_path()); itr != end_itr; ++itr) {
    if (!boost::filesystem::is_directory(itr->status())) {
      if (itr->path().extension() == ".off")  // see below
      {
        Points_off_reader off_reader(itr->path().string());
        Point p0 = off_reader.get_point_cloud()[0];

        for (Filtration_value radius = 0.1; radius < 0.4; radius += 0.1) {
          std::cout << itr->path().stem() << ";";
          std::cout << radius << ";";
          Gudhi::Clock rips_clock("Rips computation");
          Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), radius,
                                                Gudhi::Minimal_enclosing_ball_radius());
          Simplex_tree rips_stree;
          rips_complex_from_points.create_complex(rips_stree, p0.size() - 1);
          // ------------------------------------------
          // Display information about the Rips complex
          // ------------------------------------------
          double rips_sec = rips_clock.num_seconds();
          std::cout << rips_sec << ";";

          Gudhi::Clock cech_clock("Cech computation");
          Cech_complex cech_complex_from_points(off_reader.get_point_cloud(), radius);
          Simplex_tree cech_stree;
          cech_complex_from_points.create_complex(cech_stree, p0.size() - 1);
          // ------------------------------------------
          // Display information about the Cech complex
          // ------------------------------------------
          double cech_sec = cech_clock.num_seconds();
          std::cout << cech_sec << ";";
          std::cout << cech_sec / rips_sec << ";";

          assert(rips_stree.num_simplices() >= cech_stree.num_simplices());
          std::cout << rips_stree.num_simplices() << ";";
          std::cout << cech_stree.num_simplices() << ";" << std::endl;
        }
      }
    }
  }

  return 0;
}
