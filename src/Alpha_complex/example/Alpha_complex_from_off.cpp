#include <gudhi/Alpha_complex.h>
#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " filename.off alpha_square_max_value [ouput_file.txt]\n";
  std::cerr << "       i.e.: " << progName << " ../../data/points/alphacomplexdoc.off 60.0\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 3) && (argc != 4)) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct\n";
    usage(argv[0]);
  }

  std::string off_file_name(argv[1]);
  double alpha_square_max_value = atof(argv[2]);

  // ----------------------------------------------------------------------------
  // Init of an alpha complex from an OFF file
  // ----------------------------------------------------------------------------
  typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
  Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_file(off_file_name, alpha_square_max_value);

  std::streambuf* streambufffer;
  std::ofstream ouput_file_stream;
  
  if (argc == 4) {
    ouput_file_stream.open(std::string(argv[3]));
    streambufffer = ouput_file_stream.rdbuf();
  } else {
    streambufffer = std::cout.rdbuf();
  }

  std::ostream output_stream(streambufffer);

  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  output_stream << "Alpha complex is of dimension " << alpha_complex_from_file.dimension() <<
      " - " << alpha_complex_from_file.num_simplices() << " simplices - " <<
      alpha_complex_from_file.num_vertices() << " vertices." << std::endl;

  output_stream << "Iterator on alpha complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : alpha_complex_from_file.filtration_simplex_range()) {
    output_stream << "   ( ";
    for (auto vertex : alpha_complex_from_file.simplex_vertex_range(f_simplex)) {
      output_stream << vertex << " ";
    }
    output_stream << ") -> " << "[" << alpha_complex_from_file.filtration(f_simplex) << "] ";
    output_stream << std::endl;
  }
  ouput_file_stream.close();
  return 0;
}
