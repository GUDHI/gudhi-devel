#include <gudhi/Rips_complex.h>
// to construct Rips_complex from a csv file representing a distance matrix
#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <string>
#include <vector>

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " filename.csv threshold dim_max [ouput_file.txt]\n";
  std::cerr << "       i.e.: " << progName << " ../../data/distance_matrix/full_square_distance_matrix.csv 1.0 3\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if ((argc != 4) && (argc != 5)) usage(argc, (argv[0] - 1));

  std::string csv_file_name(argv[1]);
  double threshold = atof(argv[2]);
  int dim_max = atoi(argv[3]);

  // Type definitions
  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
  using Distance_matrix = std::vector<std::vector<Filtration_value>>;

  // ----------------------------------------------------------------------------
  // Init of a Rips complex from a distance matrix in a csv file
  // Default separator is ';'
  // ----------------------------------------------------------------------------
  Distance_matrix distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_file_name);
  Rips_complex rips_complex_from_file(distances, threshold);

  std::streambuf* streambufffer;
  std::ofstream ouput_file_stream;

  if (argc == 5) {
    ouput_file_stream.open(std::string(argv[4]));
    streambufffer = ouput_file_stream.rdbuf();
  } else {
    streambufffer = std::cout.rdbuf();
  }

  Simplex_tree stree;
  rips_complex_from_file.create_complex(stree, dim_max);
  std::ostream output_stream(streambufffer);

  // ----------------------------------------------------------------------------
  // Display information about the Rips complex
  // ----------------------------------------------------------------------------
  output_stream << "Rips complex is of dimension " << stree.dimension() <<
                   " - " << stree.num_simplices() << " simplices - " <<
                   stree.num_vertices() << " vertices." << std::endl;

  output_stream << "Iterator on Rips complex simplices in the filtration order, with [filtration value]:" <<
                   std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    output_stream << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      output_stream << vertex << " ";
    }
    output_stream << ") -> " << "[" << stree.filtration(f_simplex) << "] ";
    output_stream << std::endl;
  }

  ouput_file_stream.close();
  return 0;
}
