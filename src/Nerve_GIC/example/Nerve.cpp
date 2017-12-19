/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA
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

#include <gudhi/GIC.h>

#include <string>
#include <vector>

using namespace std;

void usage(int nbArgs, char *const progName) {
  cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  cerr << "Usage: " << progName << " filename.off coordinate resolution gain [--v] \n";
  cerr << "       i.e.: " << progName << " ../../data/points/human.off 2 10 0.3 --v \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 5) && (argc != 6)) usage(argc, argv[0]);

  using Point = vector<float>;

  string off_file_name(argv[1]);
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
    cout << "Incorrect OFF file." << endl;
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
    SC.compute_PD<Gudhi::Simplex_tree<> >();

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
