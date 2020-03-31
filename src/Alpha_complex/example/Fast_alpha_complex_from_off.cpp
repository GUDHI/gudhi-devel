#include <gudhi/Alpha_complex.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.off alpha_square_max_value [ouput_file.txt]\n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/alphacomplexdoc.off 60.0\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 3) && (argc != 4)) usage(argc, (argv[0] - 1));

  std::string off_file_name {argv[1]};
  double alpha_square_max_value {atof(argv[2])};

  // WARNING : CGAL::Epick_d is fast but not safe (unlike CGAL::Epeck_d)
  // (i.e. when the points are on a grid)
  using Fast_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  // ----------------------------------------------------------------------------
  // Init of an alpha complex from an OFF file
  // ----------------------------------------------------------------------------
  Gudhi::alpha_complex::Alpha_complex<Fast_kernel> alpha_complex_from_file(off_file_name);

  std::streambuf* streambuffer;
  std::ofstream ouput_file_stream;

  if (argc == 4) {
    ouput_file_stream.open(std::string(argv[3]));
    streambuffer = ouput_file_stream.rdbuf();
  } else {
    streambuffer = std::clog.rdbuf();
  }

  Gudhi::Simplex_tree<> simplex;
  if (alpha_complex_from_file.create_complex(simplex, alpha_square_max_value)) {
    std::ostream output_stream(streambuffer);

    // ----------------------------------------------------------------------------
    // Display information about the alpha complex
    // ----------------------------------------------------------------------------
    output_stream << "Alpha complex is of dimension " << simplex.dimension() <<
        " - " << simplex.num_simplices() << " simplices - " <<
        simplex.num_vertices() << " vertices." << std::endl;

    output_stream << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" <<
        std::endl;
    for (auto f_simplex : simplex.filtration_simplex_range()) {
      output_stream << "   ( ";
      for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
        output_stream << vertex << " ";
      }
      output_stream << ") -> " << "[" << simplex.filtration(f_simplex) << "] ";
      output_stream << std::endl;
    }
  }
  ouput_file_stream.close();
  return 0;
}
