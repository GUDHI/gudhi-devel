/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
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
  std::cerr << "Usage: " << progName << " filename.off coordinate resolution gain [-v] \n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/human.off 2 10 0.3 -v \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 5) && (argc != 6)) usage(argc, argv[0]);

  using Point = std::vector<float>;

  std::string off_file_name(argv[1]);
  int coord = atoi(argv[2]);
  int resolution = atoi(argv[3]);
  double gain = atof(argv[4]);
  bool verb = 0;
  if (argc == 6) verb = 1;

  // --------------------------------
  // Init of a Nerve from an OFF file
  // --------------------------------

  Gudhi::cover_complex::Cover_complex<Point> SC;
  SC.set_verbose(verb);

  bool check = SC.read_point_cloud(off_file_name);

  if (!check) {
    std::cout << "Incorrect OFF file." << std::endl;
  } else {
    SC.set_type("Nerve");

    SC.set_color_from_coordinate(coord);
    SC.set_function_from_coordinate(coord);

    SC.set_graph_from_OFF();
    SC.set_resolution_with_interval_number(resolution);
    SC.set_gain(gain);
    SC.set_cover_from_function();

    SC.find_simplices();

    SC.write_info();

    Gudhi::Simplex_tree<> stree;
    SC.create_complex(stree);
    SC.compute_PD();

    // ----------------------------------------------------------------------------
    // Display information about the graph induced complex
    // ----------------------------------------------------------------------------

    if (verb) {
      std::cout << "Nerve is of dimension " << stree.dimension() << " - " << stree.num_simplices() << " simplices - "
                << stree.num_vertices() << " vertices." << std::endl;

      std::cout << "Iterator on Nerve simplices" << std::endl;
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
