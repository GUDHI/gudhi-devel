/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Points_off_io.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_d.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>   // infinity
#include <utility>  // for pair
#include <map>

// -------------------------------------------------------------------------------
// cech_complex_cgal_mini_sphere_3d is an example of each step that is required to
// build a Cech over a Simplex_tree. Please refer to cech_persistence to see
// how to do the same thing with the Cech_complex wrapper for less detailed
// steps.
// -------------------------------------------------------------------------------

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Simplex_handle = Simplex_tree::Simplex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Siblings = Simplex_tree::Siblings;
using Graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                                      boost::property<Gudhi::vertex_filtration_t, Filtration_value>,
                                      boost::property<Gudhi::edge_filtration_t, Filtration_value> >;
using Edge_t = std::pair<Vertex_handle, Vertex_handle>;

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<3> >;
using Point = Kernel::Point_d;
using Traits = CGAL::Min_sphere_of_points_d_traits_d<Kernel, Filtration_value, 3>;
using Min_sphere = CGAL::Min_sphere_of_spheres_d<Traits>;

using Points_off_reader = Gudhi::Points_off_reader<Point>;

class Cech_blocker {
 public:
  bool operator()(Simplex_handle sh) {
    std::vector<Point> points;
#if DEBUG_TRACES
    std::cout << "Cech_blocker on [";
#endif  // DEBUG_TRACES
    for (auto vertex : simplex_tree_.simplex_vertex_range(sh)) {
      points.push_back(point_cloud_[vertex]);
#if DEBUG_TRACES
      std::cout << vertex << ", ";
#endif  // DEBUG_TRACES
    }
    Min_sphere ms(points.begin(), points.end());
    Filtration_value radius = ms.radius();
#if DEBUG_TRACES
    std::cout << "] - radius = " << radius << " - returns " << (radius > threshold_) << std::endl;
#endif  // DEBUG_TRACES
    simplex_tree_.assign_filtration(sh, radius);
    return (radius > threshold_);
  }
  Cech_blocker(Simplex_tree& simplex_tree, Filtration_value threshold, const std::vector<Point>& point_cloud)
      : simplex_tree_(simplex_tree), threshold_(threshold), point_cloud_(point_cloud) {}

 private:
  Simplex_tree simplex_tree_;
  Filtration_value threshold_;
  std::vector<Point> point_cloud_;
};

template <typename InputPointRange>
Graph_t compute_proximity_graph(InputPointRange& points, Filtration_value threshold);

void program_options(int argc, char* argv[], std::string& off_file_points, Filtration_value& threshold, int& dim_max);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  Filtration_value threshold;
  int dim_max;

  program_options(argc, argv, off_file_points, threshold, dim_max);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph(off_reader.get_point_cloud(), threshold);

  // Min_sphere sph1(off_reader.get_point_cloud()[0], off_reader.get_point_cloud()[1], off_reader.get_point_cloud()[2]);
  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;
  // insert the proximity graph in the simplex tree
  st.insert_graph(prox_graph);
  // expand the graph until dimension dim_max
  st.expansion_with_blockers(dim_max, Cech_blocker(st, threshold, off_reader.get_point_cloud()));

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

#if DEBUG_TRACES
  std::cout << "********************************************************************\n";
  // Display the Simplex_tree - Can not be done in the middle of 2 inserts
  std::cout << "* The complex contains " << st.num_simplices() << " simplices - dimension=" << st.dimension() << "\n";
  std::cout << "* Iterator on Simplices in the filtration, with [filtration value]:\n";
  for (auto f_simplex : st.filtration_simplex_range()) {
    std::cout << "   "
              << "[" << st.filtration(f_simplex) << "] ";
    for (auto vertex : st.simplex_vertex_range(f_simplex)) {
      std::cout << static_cast<int>(vertex) << " ";
    }
    std::cout << std::endl;
  }
#endif  // DEBUG_TRACES
  return 0;
}

void program_options(int argc, char* argv[], std::string& off_file_points, Filtration_value& threshold, int& dim_max) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a 3d point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "max-edge-length,r",
      po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal length of an edge for the Cech complex construction.")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Cech complex we want to compute.");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Construct a Cech complex defined on a set of input points.\n \n";

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
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
template <typename InputPointRange>
Graph_t compute_proximity_graph(InputPointRange& points, Filtration_value threshold) {
  std::vector<Edge_t> edges;
  std::vector<Filtration_value> edges_fil;

  Kernel k;
  Vertex_handle idx_u, idx_v;
  Filtration_value fil;
  idx_u = 0;
  for (auto it_u = points.begin(); it_u != points.end(); ++it_u) {
    idx_v = idx_u + 1;
    for (auto it_v = it_u + 1; it_v != points.end(); ++it_v, ++idx_v) {
      fil = k.squared_distance_d_object()(*it_u, *it_v);
      // For Cech Complex, threshold is a radius (distance /2)
      fil = std::sqrt(fil) / 2.;
      if (fil <= threshold) {
        edges.emplace_back(idx_u, idx_v);
        edges_fil.push_back(fil);
      }
    }
    ++idx_u;
  }

  Graph_t skel_graph(edges.begin(), edges.end(), edges_fil.begin(),
                     idx_u);  // number of points labeled from 0 to idx_u-1

  auto vertex_prop = boost::get(Gudhi::vertex_filtration_t(), skel_graph);

  boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}
