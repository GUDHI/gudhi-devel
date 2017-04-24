#include <gudhi/GIC.h>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off threshold function resolution gain [ouput_file.txt]\n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/test.off 1.5 test_cov \n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 6) && (argc != 7)) usage(argc, (argv[0] - 1));

  std::string off_file_name(argv[1]);
  double threshold = atof(argv[2]);
  int coord = atoi(argv[3]);
  double resolution = atof(argv[4]);
  double gain = atof(argv[5]);

  // Type definitions
  using Graph_t = boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,\
                                          boost::property < vertex_filtration_t, Filtration_value >,\
                                          boost::property < edge_filtration_t, Filtration_value > >;

  // ----------------------------------------------------------------------------
  // Init of a graph induced complex from an OFF file
  // ----------------------------------------------------------------------------

  Gudhi::graph_induced_complex::Graph_induced_complex GIC;

  GIC.set_graph_from_rips(threshold, off_file_name);
  //GIC.set_graph_from_OFF(off_file_name);
  GIC.set_function_from_coordinate(coord, off_file_name);
  GIC.set_cover_from_function(resolution,gain,0);
  GIC.find_GIC_simplices();
  //GIC.find_GIC_simplices_with_functional_minimal_cover();
  Simplex_tree stree; GIC.create_complex(stree);

  std::streambuf* streambufffer;
  std::ofstream ouput_file_stream;

  if (argc == 7) {
    ouput_file_stream.open(std::string(argv[4]));
    streambufffer = ouput_file_stream.rdbuf();
  } else {
    streambufffer = std::cout.rdbuf();
  }

  std::ostream output_stream(streambufffer);

  // ----------------------------------------------------------------------------
  // Display information about the graph induced complex
  // ----------------------------------------------------------------------------
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

  ouput_file_stream.close();

  return 0;
}
