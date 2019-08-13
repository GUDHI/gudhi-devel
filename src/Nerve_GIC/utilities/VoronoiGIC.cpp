/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/GIC.h>

#include <string>
#include <vector>

void usage(int nbArgs, char *const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off N [-v] \n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/human.off 100 -v \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 3) && (argc != 4)) usage(argc, argv[0]);

  using Point = std::vector<float>;

  std::string off_file_name(argv[1]);
  int m = atoi(argv[2]);
  bool verb = 0;
  if (argc == 4) verb = 1;

  // ----------------------------------------------------------------------------
  // Init of a graph induced complex from an OFF file
  // ----------------------------------------------------------------------------

  Gudhi::cover_complex::Cover_complex<Point> GIC;
  GIC.set_verbose(verb);

  bool check = GIC.read_point_cloud(off_file_name);

  if (!check) {
    std::cout << "Incorrect OFF file." << std::endl;
  } else {
    GIC.set_type("GIC");

    GIC.set_color_from_coordinate();

    GIC.set_graph_from_OFF();
    GIC.set_cover_from_Voronoi(Gudhi::Euclidean_distance(), m);

    GIC.find_simplices();

    GIC.plot_OFF();

    Gudhi::Simplex_tree<> stree;
    GIC.create_complex(stree);

    // ----------------------------------------------------------------------------
    // Display information about the graph induced complex
    // ----------------------------------------------------------------------------

    if (verb) {
      std::cout << "Graph induced complex is of dimension " << stree.dimension() << " - " << stree.num_simplices()
                << " simplices - " << stree.num_vertices() << " vertices." << std::endl;

      std::cout << "Iterator on graph induced complex simplices" << std::endl;
      for (auto f_simplex : stree.filtration_simplex_range()) {
        for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
          std::cout << vertex << " ";
        }
        std::cout << std::endl;
      }
    }
  }

  return 0;
}
