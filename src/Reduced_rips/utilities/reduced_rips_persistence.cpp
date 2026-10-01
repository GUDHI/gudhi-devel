/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <gudhi/Reduced_rips.h>
#include <gudhi/Points_off_io.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Clock.h>

namespace po = boost::program_options;

using Reduced_rips = Gudhi::reduced_rips::Reduced_rips<>;

int main(int argc, char** argv) {
  std::string off_file;
  std::string matrix_file;
  std::string output_file;
  unsigned num_neighbors = 0;
  std::string search = "auto";
  double min_persistence = 0.0;

  po::options_description visible("Allowed options");
  visible.add_options()("help,h", "produce help message")(
      "distance-matrix,d", po::value<std::string>(&matrix_file)->default_value(""),
      "read a lower-triangular distance matrix (';'-separated CSV) instead of an OFF point cloud")(
      "output-file,o", po::value<std::string>(&output_file)->default_value(""),
      "name of an output file for the persistence diagram (default: standard output)")(
      "num-neighbors,k", po::value<unsigned>(&num_neighbors)->default_value(0),
      "initial neighbor budget per point (0 = sqrt(n))")(
      "search,s", po::value<std::string>(&search)->default_value("auto"),
      "spatial search strategy for the point-cloud input: auto | kd | brute (ignored for a distance matrix)")(
      "min-persistence,m", po::value<double>(&min_persistence)->default_value(0.0),
      "minimal lifetime (death - birth) of a bar to be recorded (the computation never emits zero-length bars)");

  po::options_description hidden("Hidden options");
  hidden.add_options()("input-file", po::value<std::string>(&off_file), "input OFF point-cloud file");
  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  bool have_matrix = !matrix_file.empty();
  if (have_matrix && vm.count("input-file") != 0U)
    std::cerr << "Warning: both a distance matrix (-d) and an OFF file were given; the OFF file is ignored.\n";
  if ((vm.count("help") != 0U) || (!have_matrix && (vm.count("input-file") == 0U))) {
    std::cout << "Usage: " << argv[0] << " [options] <input OFF file>\n";
    std::cout << "   or: " << argv[0] << " [options] -d <distance matrix CSV>\n\n";
    std::cout << "Computes the degree-1 Vietoris-Rips persistence diagram of a Euclidean point cloud or an\n";
    std::cout << "arbitrary symmetric distance matrix, using the Reduced Vietoris-Rips filtration.\n\n";
    std::cout << visible << "\n";
    return (vm.count("help") != 0U) ? 0 : 1;
  }

  Reduced_rips::Search method = Reduced_rips::Search::automatic;
  if (search == "kd")
    method = Reduced_rips::Search::kd_tree;
  else if (search == "brute")
    method = Reduced_rips::Search::brute_force;
  else if (search != "auto") {
    std::cerr << "Unknown --search value: " << search << " (use auto | kd | brute)\n";
    return 1;
  }

  Reduced_rips ph1;
  if (have_matrix) {
    auto distances = Gudhi::read_lower_triangular_matrix_from_csv_file<double>(matrix_file);
    if (distances.size() < 2) {
      std::cerr << "Could not read a distance matrix of at least two points from: " << matrix_file << '\n';
      return 1;
    }
    Gudhi::Clock clock("Reduced Vietoris-Rips degree-1 persistence");
    ph1 = Reduced_rips::from_distance_matrix(distances, num_neighbors, method);
    clock.end();
    std::clog << clock;
  } else {
    Gudhi::Points_off_reader<std::vector<double>> off_reader(off_file);
    if (!off_reader.is_valid()) {
      std::cerr << "Could not read OFF file: " << off_file << '\n';
      return 1;
    }
    Gudhi::Clock clock("Reduced Vietoris-Rips degree-1 persistence");
    ph1 = Reduced_rips::from_points(off_reader.get_point_cloud(), num_neighbors, method);
    clock.end();
    std::clog << clock;
  }

  auto barcode = ph1.persistence();
  // Return with longest bars first, rather than the engine's ascending-by-death order.
  std::sort(barcode.begin(), barcode.end(),
            [](const auto& x, const auto& y) { return (x[1] - x[0]) > (y[1] - y[0]); });

  std::ostream* out = &std::cout;
  std::ofstream ofs;
  if (!output_file.empty()) {
    ofs.open(output_file);
    out = &ofs;
  }
  // One bar per line in the GUDHI utilities convention `p dim birth death`. Coefficients are Z/2Z, so p = 2.
  for (const auto& bar : barcode)
    if (bar[1] - bar[0] > min_persistence) *out << "2 1 " << bar[0] << " " << bar[1] << "\n";

  std::clog << "1-simplices: " << ph1.num_one_simplices() << ", 2-simplices: " << ph1.num_two_simplices()
            << ", persistent pairs: " << ph1.num_persistence_pairs() << '\n';
  return 0;
}
