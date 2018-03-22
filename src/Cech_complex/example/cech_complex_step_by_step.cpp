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

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_off_io.h>

#include <gudhi/Miniball.hpp>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity
#include <utility>  // for pair
#include <map>

// ----------------------------------------------------------------------------
// rips_persistence_step_by_step is an example of each step that is required to
// build a Rips over a Simplex_tree. Please refer to rips_persistence to see
// how to do the same thing with the Rips_complex wrapper for less detailed
// steps.
// ----------------------------------------------------------------------------

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<>;
using Simplex_handle = Simplex_tree::Simplex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Proximity_graph = Gudhi::Proximity_graph<Simplex_tree>;

class Cech_blocker {
 private:
  using Point_cloud = std::vector<Point>;
  using Point_iterator = Point_cloud::const_iterator;
  using Coordinate_iterator = Point::const_iterator;
  using Min_sphere =  Gudhi::Miniball::Miniball <Gudhi::Miniball::CoordAccessor<Point_iterator, Coordinate_iterator>>;
 public:
  bool operator()(Simplex_handle sh) {
    std::vector<Point> points;
    for (auto vertex : simplex_tree_.simplex_vertex_range(sh)) {
      points.push_back(point_cloud_[vertex]);
#ifdef DEBUG_TRACES
      std::cout << "#(" << vertex << ")#";
#endif  // DEBUG_TRACES
    }
    Filtration_value radius = Gudhi::Minimal_enclosing_ball_radius()(points);
#ifdef DEBUG_TRACES
    std::cout << "radius = " << radius << " - " << (radius > max_radius_) << std::endl;
#endif  // DEBUG_TRACES
    simplex_tree_.assign_filtration(sh, radius);
    return (radius > max_radius_);
  }
  Cech_blocker(Simplex_tree& simplex_tree, Filtration_value max_radius, const std::vector<Point>& point_cloud)
    : simplex_tree_(simplex_tree),
      max_radius_(max_radius),
      point_cloud_(point_cloud) {
    dimension_ = point_cloud_[0].size();
  }
 private:
  Simplex_tree simplex_tree_;
  Filtration_value max_radius_;
  std::vector<Point> point_cloud_;
  int dimension_;
};


void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , Filtration_value & max_radius
                     , int & dim_max);


int main(int argc, char * argv[]) {
  std::string off_file_points;
  Filtration_value max_radius;
  int dim_max;

  program_options(argc, argv, off_file_points, max_radius, dim_max);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  Proximity_graph prox_graph = Gudhi::compute_proximity_graph<Simplex_tree>(off_reader.get_point_cloud(),
                                                                            max_radius,
                                                                            Gudhi::Minimal_enclosing_ball_radius());

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion_with_blockers(dim_max, Cech_blocker(st, max_radius, off_reader.get_point_cloud()));

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

#if DEBUG_TRACES
  std::cout << "********************************************************************\n";
  std::cout << "* The complex contains " << st.num_simplices() << " simplices - dimension=" << st.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << static_cast<int>(vertex) << " ";
    }
    std::cout << std::endl;
  }
#endif  // DEBUG_TRACES

  return 0;
}

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , Filtration_value & max_radius
                     , int & dim_max) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&off_file_points),
       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()
      ("help,h", "produce help message")
      ("max-radius,r",
       po::value<Filtration_value>(&max_radius)->default_value(std::numeric_limits<Filtration_value>::infinity()),
       "Maximal length of an edge for the Rips complex construction.")
      ("cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
       "Maximal dimension of the Rips complex we want to compute.");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Construct a Cech complex defined on a set of input points.\n \n";

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}
