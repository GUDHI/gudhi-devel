/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

// #include <CGAL/Epick_d.h>
// #include <CGAL/Euclidean_distance.h>
// #include <CGAL/Search_traits_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Optimisation_d_traits_d.h>
#include <CGAL/Min_sphere_d.h>

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
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Simplex_handle = Simplex_tree::Simplex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Siblings = Simplex_tree::Siblings;
using Graph_t = boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
, boost::property < vertex_filtration_t, Filtration_value >
, boost::property < edge_filtration_t, Filtration_value >
>;
using Edge_t = std::pair< Vertex_handle, Vertex_handle >;

// using Kernel = CGAL::Epick_d< CGAL::Dimension_tag<2> >;// CGAL::Dynamic_dimension_tag >;
typedef CGAL::Cartesian_d<double>              Kernel;
typedef CGAL::Optimisation_d_traits_d<Kernel>  Traits;
typedef CGAL::Min_sphere_d<Traits>             Min_sphere;

using Point = Kernel::Point_d;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
// using Min_sphere = CGAL::Min_sphere_d<Kernel>;

class Cech_blocker {
 public:
  std::pair<bool, Filtration_value> operator()(Simplex_handle origin_sh, Simplex_handle dict1_sh, Simplex_handle dict2_sh, Siblings* siblings) {
    //std::vector<Vertex_handle> path = {dict1_sh->first, origin_sh->first};
    Siblings* sib_path = siblings;
    std::vector<Point> sphere_points = {point_cloud_[dict1_sh->first], point_cloud_[origin_sh->first]};
    do {
      //path.push_back(sib_path->parent());
      sphere_points.push_back(point_cloud_[sib_path->parent()]);
      sib_path = sib_path->oncles();
    } while (sib_path->oncles() != nullptr);
    /*std::cout << square_threshold_ << "-";
    for (auto vh : path) {
      std::cout << vh << " ";
    }
    std::cout << std::endl;*/
    Min_sphere min_sphere(sphere_points.begin(), sphere_points.end());
    //std::cout << min_sphere.squared_radius() << std::endl;
    Filtration_value squared_diameter = min_sphere.squared_radius() * 4.;
    // Default blocker is always insert with the maximal filtration value between
    // origin, dict1 and dict2
    return std::make_pair(squared_diameter < square_threshold_,
                          squared_diameter);
  }
  Cech_blocker(Filtration_value threshold, const std::vector<Point>& point_cloud)
    : square_threshold_(threshold * threshold),
      point_cloud_(point_cloud) { }
 private:
  Filtration_value square_threshold_;
  std::vector<Point> point_cloud_;
};

template< typename InputPointRange>
Graph_t compute_proximity_graph(InputPointRange &points, Filtration_value threshold);

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , Filtration_value & threshold
                     , int & dim_max);

int main(int argc, char * argv[]) {
  std::string off_file_points;
  Filtration_value threshold;
  int dim_max;

  program_options(argc, argv, off_file_points, threshold, dim_max);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(off_reader.get_point_cloud(), threshold);

  //Min_sphere sph1(off_reader.get_point_cloud()[0], off_reader.get_point_cloud()[1], off_reader.get_point_cloud()[2]);
  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion_with_blockers(dim_max, Cech_blocker(threshold, off_reader.get_point_cloud()));

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  std::cout << "********************************************************************\n";
  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "* The complex contains " << st.num_simplices() << " simplices\n";
  std::cout << "   - dimension " << st.dimension() << "   - filtration " << st.filtration() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << static_cast<int>(vertex) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}

void program_options(int argc, char * argv[]
                     , std::string & off_file_points
                     , Filtration_value & threshold
                     , int & dim_max) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("input-file", po::value<std::string>(&off_file_points),
       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()
      ("help,h", "produce help message")
      ("max-edge-length,r",
       po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
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

/** Output the proximity graph of the points.
 *
 * If points contains n elements, the proximity graph is the graph
 * with n vertices, and an edge [u,v] iff the distance function between
 * points u and v is smaller than threshold.
 *
 * The type PointCloud furnishes .begin() and .end() methods, that return
 * iterators with value_type Point.
 */
template< typename InputPointRange>
Graph_t compute_proximity_graph(InputPointRange &points, Filtration_value threshold) {
  std::vector< Edge_t > edges;
  std::vector< Filtration_value > edges_fil;

  Kernel k;
  Filtration_value square_threshold = threshold * threshold;

  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = k.squared_distance_d_object()(*it_u, *it_v);
      if (fil <= square_threshold) {
        edges.emplace_back(idx_u, idx_v);
        edges_fil.push_back(fil);
      }
    }
    ++idx_u;
  }

  Graph_t skel_graph(edges.begin()
                     , edges.end()
                     , edges_fil.begin()
                     , idx_u);  // number of points labeled from 0 to idx_u-1

  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph);
       vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}
