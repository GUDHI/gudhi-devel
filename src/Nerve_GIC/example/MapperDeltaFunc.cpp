/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA Saclay (France)
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

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off function [--v] \n";
  std::cerr << "       i.e.: " << progName << " ../../../../data/points/COIL_database/lucky_cat.off ../../../../data/points/COIL_database/lucky_cat_PCA1 --v \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 3) && (argc != 4)) usage(argc, (argv[0] - 1));

  std::string off_file_name(argv[1]);
  std::string func_file_name = argv[2];
  bool verb = 0; if(argc == 4)  verb = 1;

  // ---------------------------------------
  // Init of a Mapper Delta from an OFF file
  // ---------------------------------------

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;
  GIC.set_verbose(verb);

  bool check = GIC.read_point_cloud(off_file_name);

  if(!check)  std::cout << "Incorrect OFF file." << std::endl;
  else{

    GIC.set_color_from_file(func_file_name);
    GIC.set_function_from_file(func_file_name);

    GIC.set_graph_from_automatic_rips();

    GIC.set_automatic_resolution_for_GICMAP(); GIC.set_gain();
    GIC.set_cover_from_function(1);

    GIC.find_GICMAP_simplices_with_functional_minimal_cover();

    GIC.plot_pdf();

    Simplex_tree stree; GIC.create_complex(stree);

    std::streambuf* streambufffer = std::cout.rdbuf();
    std::ostream output_stream(streambufffer);

    // ------------------------------------------
    // Display information about the Mapper Delta
    // ------------------------------------------

    if(verb){
      output_stream << "Mapper Delta is of dimension " << stree.dimension() <<
                     " - " << stree.num_simplices() << " simplices - " <<
                     stree.num_vertices() << " vertices." << std::endl;

      output_stream << "Iterator on Mapper Delta simplices" << std::endl;
      for (auto f_simplex : stree.filtration_simplex_range()) {
        for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
          output_stream << vertex << " ";
        }
        output_stream << std::endl;
      }
    }
  }

  return 0;
}
