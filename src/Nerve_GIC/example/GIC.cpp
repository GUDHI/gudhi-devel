#include <gudhi/GIC.h>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off threshold coordinate resolution gain\n";
  std::cerr << "       i.e.: " << progName << " ../../../data/points/human.off 0.075 2 10 0.3 \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 6) && (argc != 7)) usage(argc, (argv[0] - 1));

  std::string off_file_name(argv[1]);
  double threshold = atof(argv[2]);
  int coord = atoi(argv[3]);
  double resolution = atof(argv[4]);
  double gain = atof(argv[5]);
  bool verb = 0; if(argc == 7)  verb = 1;

  // ----------------------------------------------------------------------------
  // Init of a graph induced complex from an OFF file
  // ----------------------------------------------------------------------------

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;
  GIC.set_verbose(verb);

  GIC.set_color_from_coordinate(off_file_name, coord);
  GIC.set_function_from_coordinate(coord, off_file_name);

  GIC.set_graph_from_rips(threshold, off_file_name);

  GIC.set_resolution_double(resolution); GIC.set_gain(gain);
  GIC.set_cover_from_function(1);

  GIC.find_GIC_simplices();

  GIC.plot_with_KeplerMapper();

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
