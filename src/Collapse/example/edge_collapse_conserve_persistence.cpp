/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Flag_complex_edge_collapser.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/graph_simplicial_complex.h>

#include <boost/range/adaptor/transformed.hpp>

#include<utility>  // for std::pair
#include<vector>
#include<tuple>

// Types definition

using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Point = std::vector<Filtration_value>;
using Vector_of_points = std::vector<Point>;

using Flag_complex_edge_collapser = Gudhi::collapse::Flag_complex_edge_collapser<Vertex_handle, Filtration_value>;
using Proximity_graph = Gudhi::Proximity_graph<Flag_complex_edge_collapser>;

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

using Persistence_interval = std::tuple<int, Filtration_value, Filtration_value>;
/*
 * Compare two intervals by dimension, then by length.
 */
struct cmp_intervals_by_length {
  explicit cmp_intervals_by_length(Simplex_tree * sc)
      : sc_(sc) { }

  template<typename Persistent_interval>
  bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
    return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
            > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
  }
  Simplex_tree* sc_;
};

std::vector<Persistence_interval> get_persistence_intervals(Simplex_tree& st, int ambient_dim) {
  std::vector<Persistence_interval> persistence_intervals;
  st.expansion(ambient_dim);
  
  // Sort the simplices in the order of the filtration
  st.initialize_filtration();
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(st);
  // initializes the coefficient field for homology - must be a prime number
  int p = 11;
  pcoh.init_coefficients(p);

  // Default min_interval_length = 0.
  pcoh.compute_persistent_cohomology();
  // Custom sort and output persistence
  cmp_intervals_by_length cmp(&st);
  auto persistent_pairs = pcoh.get_persistent_pairs();
  std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
  for (auto pair : persistent_pairs) {
    persistence_intervals.emplace_back(st.dimension(get<0>(pair)),
                                       st.filtration(get<0>(pair)),
                                       st.filtration(get<1>(pair)));
  }
  return persistence_intervals;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "This program requires an OFF file and minimal threshold value as parameter\n";
    std::cerr << "For instance: ./Edge_collapse_conserve_persistence ../../data/points/tore3D_300.off 1.\n";
    exit(-1);  // ----- >>
  }

  std::string off_file_points {argv[1]};
  double threshold {atof(argv[2])};

  Gudhi::Points_off_reader<Point> off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file_points << "\n";
    exit(-1);  // ----- >>
  }

  Vector_of_points point_vector = off_reader.get_point_cloud();
  if (point_vector.size() <= 0) {
    std::cerr << "Empty point cloud." << std::endl;
    exit(-1);  // ----- >>
  }

  Proximity_graph proximity_graph = Gudhi::compute_proximity_graph<Simplex_tree>(off_reader.get_point_cloud(),
                                                                                 threshold,
                                                                                 Gudhi::Euclidean_distance());

  if (num_edges(proximity_graph) <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  int ambient_dim = point_vector[0].size();

  // ***** Simplex tree from a flag complex built after collapse *****
  Flag_complex_edge_collapser edge_collapser(
    boost::adaptors::transform(edges(proximity_graph), [&](auto&&edge){
      return std::make_tuple(source(edge, proximity_graph),
                             target(edge, proximity_graph),
                             get(Gudhi::edge_filtration_t(), proximity_graph, edge));
      })
  );

  Simplex_tree stree_from_collapse;
  for (Vertex_handle vertex = 0; static_cast<std::size_t>(vertex) < point_vector.size(); vertex++) {
    // insert the vertex with a 0. filtration value just like a Rips
    stree_from_collapse.insert_simplex({vertex}, 0.);
  }
  edge_collapser.process_edges(
    [&stree_from_collapse](Vertex_handle u, Vertex_handle v, Filtration_value filtration) {
        // insert the edge
        stree_from_collapse.insert_simplex({u, v}, filtration);
      });

  std::vector<Persistence_interval> persistence_intervals_from_collapse = get_persistence_intervals(stree_from_collapse, ambient_dim);

  // ***** Simplex tree from the complete flag complex *****
  Simplex_tree stree_wo_collapse;
  stree_wo_collapse.insert_graph(proximity_graph);

  std::vector<Persistence_interval> persistence_intervals_wo_collapse = get_persistence_intervals(stree_wo_collapse, ambient_dim);

  // ***** Comparison *****
  if (persistence_intervals_wo_collapse.size() != persistence_intervals_from_collapse.size()) {
    std::cerr << "Number of persistence pairs with    collapse is " << persistence_intervals_from_collapse.size() << std::endl;
    std::cerr << "Number of persistence pairs without collapse is " << persistence_intervals_wo_collapse.size()   << std::endl;
    exit(-1);
  }

  int return_value = 0;
  auto ppwoc_ptr = persistence_intervals_wo_collapse.begin();
  for (auto ppfc: persistence_intervals_from_collapse) {
    if (ppfc != *ppwoc_ptr) {
      return_value++;
      std::cerr << "Without collapse: "
                << std::get<0>(*ppwoc_ptr) << " " << std::get<1>(*ppwoc_ptr) << " " << std::get<2>(*ppwoc_ptr)
                << " - With collapse: "
                << std::get<0>(ppfc) << " " << std::get<1>(ppfc) << " " << std::get<2>(ppfc) << std::endl;
    }
    ppwoc_ptr++;
  }
  return return_value;
}
