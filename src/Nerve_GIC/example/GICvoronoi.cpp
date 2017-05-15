#include <gudhi/GIC.h>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off N [--v] \n";
  std::cerr << "       i.e.: " << progName << " ../../../../data/points/human.off 100 --v \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 3) && (argc != 4)) usage(argc, (argv[0] - 1));

  std::string off_file_name(argv[1]);
  int m = atoi(argv[2]);
  bool verb = 0; if(argc == 4)  verb = 1;

  // ----------------------------------------------------------------------------
  // Init of a graph induced complex from an OFF file
  // ----------------------------------------------------------------------------

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;
  GIC.set_verbose(verb);

  GIC.read_point_cloud(off_file_name);

  GIC.set_color_from_coordinate();

  GIC.set_graph_from_OFF(off_file_name);

  GIC.set_cover_from_Voronoi(m);

  GIC.find_GIC_simplices();

  GIC.plot_with_Geomview();

  Simplex_tree stree; GIC.create_complex(stree);

  std::streambuf* streambufffer = std::cout.rdbuf();
  std::ostream output_stream(streambufffer);

  // ----------------------------------------------------------------------------
  // Display information about the graph induced complex
  // ----------------------------------------------------------------------------

  if(verb){
    output_stream << "Graph induced complex is of dimension " << stree.dimension() <<
                   " - " << stree.num_simplices() << " simplices - " <<
                   stree.num_vertices() << " vertices." << std::endl;

    output_stream << "Iterator on graph induced complex simplices" << std::endl;
    for (auto f_simplex : stree.filtration_simplex_range()) {
      for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
        output_stream << vertex << " ";
      }
      output_stream << std::endl;
    }
  }

  return 0;
}
