#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

void usage(int nbArgs, char * const progName) {
  std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
  std::cerr << "Usage: " << progName << " threshold\n";
  std::cerr << "       i.e.: " << progName << " 12.0\n";
  exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 2) usage(argc, argv[0]);

  double threshold = atof(argv[1]);

  // Type definitions
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
  using Distance_matrix = std::vector<std::vector<Filtration_value>>;

  // User defined distance matrix is:
  // | 0    0.94 0.77 0.99 0.11  |
  // | 0.94 0    0.26 0.99 0.39  |
  // | 0.77 0.26 0    0.28 0.97  |
  // | 0.99 0.99 0.28 0    0.30  |
  // | 0.11 0.39 0.97 0.30 0     |

  Distance_matrix distances;
  distances.push_back({});
  distances.push_back({0.94});
  distances.push_back({0.77, 0.26});
  distances.push_back({0.99, 0.99, 0.28});
  distances.push_back({0.11, 0.39, 0.97, 0.30});

  // ----------------------------------------------------------------------------
  // Init of a rips complex from points
  // ----------------------------------------------------------------------------
  Rips_complex rips_complex_from_points(distances, threshold);

  Simplex_tree simplex;
  if (rips_complex_from_points.create_complex(simplex, 1)) {
    // ----------------------------------------------------------------------------
    // Display information about the one skeleton rips complex
    // ----------------------------------------------------------------------------
    std::cout << "Rips complex is of dimension " << simplex.dimension() <<
        " - " << simplex.num_simplices() << " simplices - " <<
        simplex.num_vertices() << " vertices." << std::endl;

    std::cout << "Iterator on rips complex simplices in the filtration order, with [filtration value]:" <<
        std::endl;
    for (auto f_simplex : simplex.filtration_simplex_range()) {
      std::cout << "   ( ";
      for (auto vertex : simplex.simplex_vertex_range(f_simplex)) {
        std::cout << vertex << " ";
      }
      std::cout << ") -> " << "[" << simplex.filtration(f_simplex) << "] ";
      std::cout << std::endl;
    }
  }
  return 0;
}
