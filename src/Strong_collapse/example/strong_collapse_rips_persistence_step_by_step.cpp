/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
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
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/Flag_complex_sparse_matrix.h>
#include <gudhi/Tower_assembler.h>
#include <gudhi/Rips_complex.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>   // infinity
#include <utility>  // for pair
#include <map>

// -------------------------------------------------------------------------------------------------------------------
// strong_collapse_rips_persistence_step_by_step is an example of each step that is required to collapse a Rips over a
// Simplex_tree in prder to compute its persistence
// -------------------------------------------------------------------------------------------------------------------

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using Proximity_graph = Gudhi::Proximity_graph<Simplex_tree>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string filediag;
  Filtration_value threshold;
  int dim_max;
  int p;
  Filtration_value min_persistence;

  program_options(argc, argv, off_file_points, filediag, threshold, dim_max, p, min_persistence);

  // Extract the points from the file filepoints
  Points_off_reader off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file " << off_file_points << "\n";
    exit(-1);  // ----- >>
  }

  std::size_t number_of_points = off_reader.get_point_cloud().size();
  // Compute the proximity graph of the points
  Gudhi::Filtered_edges_vector<Simplex_tree> edge_graph =
      Gudhi::compute_edge_graph<Gudhi::Filtered_edges_vector, Simplex_tree>(off_reader.get_point_cloud(), threshold,
                                                                            Gudhi::Euclidean_distance());

  Gudhi::strong_collapse::Tower_assembler twr_assembler(number_of_points);

  // The pipeline is:
  //
  //    alpha_1            <          alpha_2                      < ...               < alpha_10
  //                                                                                       |
  //                                                                                      \ /
  //                                                                           compute_edge_graph(alpha_10 = threshold)
  //                                                                                       |
  //       +-------------------------------------------------------------------------------+
  //       |                             |                         ...                     |
  //      \ /                           \ /                                               \ /
  // Edge_graph(alpha_1)   <       Edge_graph(alpha_2)             < ...           < Edge_graph(alpha_10 = threshold)
  //
  std::size_t edge_graph_size = edge_graph.size();

  // Creates a vector of index to sub filter on - just change the vector size to modify the number of loops
  std::vector<std::size_t> indices(4);
  std::size_t indices_size = indices.size();
  std::size_t index_iterator = 0;
  std::generate(indices.begin(), indices.end(),
                [=] () mutable { index_iterator++;
                                 return index_iterator * edge_graph_size / indices_size;
                               });
  // Force the last index value because of division rounding
  indices[indices_size] = edge_graph_size;

  std::cout << "min=" << edge_graph.get_filtration_min() << " - max=" << edge_graph.get_filtration_max()
            << " - index_step=" << indices[0] << " - size=" << edge_graph_size << std::endl;

  Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_prev_coll(number_of_points);

  for (std::size_t index : indices) {
    std::cout << "index=" << index << std::endl;
    Gudhi::strong_collapse::Flag_complex_sparse_matrix mat_coll(number_of_points,
                                                                edge_graph.sub_filter_edges_by_index(index));

    mat_coll.strong_collapse();
    Gudhi::strong_collapse::Map redmap = mat_coll.reduction_map();

    twr_assembler.build_tower_for_two_cmplxs(mat_prev_coll, mat_coll, redmap, threshold);
    mat_prev_coll = mat_coll;
  }
  Gudhi::strong_collapse::Distance_matrix sparse_distances = twr_assembler.distance_matrix();

  Rips_complex rips_complex_after_collapse(sparse_distances, threshold);

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;

  rips_complex_after_collapse.create_complex(st, dim_max);

  std::cout << "The complex contains " << st.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << st.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(st);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(p);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in filediag
  if (filediag.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(filediag);
    pcoh.output_diagram(out);
    out.close();
  }

  return 0;
}

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file_points),
                       "Name of an OFF file containing a point set.\n");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::cout")(
      "max-edge-length,r",
      po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal length of an edge for the Rips complex construction.")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Rips complex we want to compute.")(
      "field-charac,p", po::value<int>(&p)->default_value(11),
      "Characteristic p of the coefficient field Z/pZ for computing homology.")(
      "min-persistence,m", po::value<Filtration_value>(&min_persistence),
      "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
      "intervals");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::cout << "of a Rips complex defined on a set of input points.\n \n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
