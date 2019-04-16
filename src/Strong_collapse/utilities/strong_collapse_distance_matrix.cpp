/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
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

#include <gudhi/Points_off_io.h>
#include <gudhi/reader_utils.h>
#include <gudhi/distance_functions.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Flag_complex_strong_collapse.h>

#include <boost/program_options.hpp>

#include <string>
#include <vector>
#include <limits>  // infinity
#include <algorithm>  // for std::sort

// Types definition
struct Simplicial_complex {
  using Filtration_value = double;
  using Vertex_handle = int;
};
using Filtration_value = Simplicial_complex::Filtration_value;
using Point = std::vector<Filtration_value>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Flag_complex_strong_collapse = Gudhi::strong_collapse::Flag_complex_strong_collapse;
using Filtration_set = std::vector<Filtration_value>;

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& output_file,
                     bool& dry_run, Filtration_set& edge_length_set);

int main(int argc, char* argv[]) {
  std::string off_file_points;
  std::string output_file;
  bool dry_run = false;
  Filtration_set edge_length_set;

  program_options(argc, argv, off_file_points, output_file, dry_run, edge_length_set);
  std::sort (edge_length_set.begin(), edge_length_set.end());

  Points_off_reader off_reader(off_file_points);
  if (!off_reader.is_valid()) {
    std::cout << "OFF file " << off_file_points << "is not valid." << std::endl;
    exit(-1);
  }

  std::cout << "Number of points=" << off_reader.get_point_cloud().size() << std::endl;

  Gudhi::Filtered_edges_vector<Simplicial_complex> edge_graph =
      Gudhi::compute_edge_graph<Gudhi::Filtered_edges_vector, Simplicial_complex>(
          off_reader.get_point_cloud(),
          edge_length_set.back(),
          Gudhi::Euclidean_distance());

  std::cout << "Edge graph contains " << edge_graph.size() << " edges. ";
  std::cout << "Minimal edge length is '" << edge_graph.get_filtration_min() << "'. ";
  std::cout << "Maximal edge length is '" << edge_graph.get_filtration_max() << "'." << std::endl;

  if (dry_run) {
    exit(0);
  }

  Flag_complex_strong_collapse flag_complex(off_reader.get_point_cloud().size());
  if (edge_length_set.size() == 1)
    flag_complex.initialize_exact_version(edge_graph);
  else
    flag_complex.initialize_approximate_version(edge_graph, edge_length_set);

  auto distance_matrix = flag_complex.get_distance_matrix();

  std::ostream* output_stream = &std::cout;
  std::ofstream file_out;
  if (!output_file.empty()) {
    file_out.open(output_file);
    output_stream = &file_out;
  }

  for (auto line : distance_matrix) {
    for (auto value : line) {
      *output_stream << value << ";";
    }
    *output_stream << std::endl;
  }

  if (!output_file.empty()) {
    file_out.close();
  }


  return 0;
}

void program_options(int argc, char* argv[], std::string& off_file_points, std::string& output_file,
                     bool& dry_run, Filtration_set& edge_length_set) {
  namespace po = boost::program_options;
  po::options_description visible("Allowed options", 100);

  visible.add_options()("help,h", "produce help message")(
      "input-file,i", po::value<std::string>(&off_file_points),
      "Name of an OFF file containing a point set.")(
      "output-file,o", po::value<std::string>(&output_file)->default_value(std::string()),
      "Name of file in which the distance matrix is written. Default print in std::cout")(
      "dry-run",
      "If the dry-run is set, only a summary of the edge graph will be displayed. This can help to define edge-length-set")(
      "edge-length-set,s", po::value<Filtration_set>(&edge_length_set)->multitoken()->default_value(
          Filtration_set{std::numeric_limits<Filtration_value>::infinity()}, "{inf}"),
      "Edge length set. Default computes exact version. i.e. '-s 0.2 0.3 0.4'");

  po::options_description all;
  all.add(visible);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  dry_run = vm.count("dry-run");

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the strong collapsed distance matrix \n";
    std::cout << "of a set of input points." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] -i input-off-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
