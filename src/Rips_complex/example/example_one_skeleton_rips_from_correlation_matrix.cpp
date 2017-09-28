#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

int main() {
  // Type definitions
  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
  using Distance_matrix = std::vector<std::vector<Filtration_value>>;

  // User defined correlation matrix is:
  // |1	    0.06	0.23	0.01	0.89|
  // |0.06	1	    0.74	0.01	0.61|
  // |0.23	0.74	1	    0.72	0.03|
  // |0.01	0.01	0.72	1	    0.7 |
  // |0.89	0.61	0.03	0.7	    1   |

  Distance_matrix correlations;
  correlations.push_back({});
  correlations.push_back({0.06});
  correlations.push_back({0.23, 0.74});
  correlations.push_back({0.01, 0.01, 0.72});
  correlations.push_back({0.89, 0.61, 0.03, 0.7});

  // ----------------------------------------------------------------------------
  // Convert correlation matrix to a distance matrix:
  // ----------------------------------------------------------------------------
  for (size_t i = 0; i != correlations.size(); ++i) {
    for (size_t j = 0; j != correlations[i].size(); ++j) {
      correlations[i][j] = 1 - correlations[i][j];
      if (correlations[i][j] < 0) {
        std::cerr << "The input matrix is not a correlation matrix. \n";
        throw "The input matrix is not a correlation matrix. \n";
      }
    }
  }

  //-----------------------------------------------------------------------------
  // Now the correlation matrix is really the distance matrix and can be processed further.
  //-----------------------------------------------------------------------------
  Distance_matrix distances = correlations;

  double threshold = 1.0;
  Rips_complex rips_complex_from_points(distances, threshold);

  Simplex_tree stree;
  rips_complex_from_points.create_complex(stree, 1);
  // ----------------------------------------------------------------------------
  // Display information about the one skeleton Rips complex
  // ----------------------------------------------------------------------------
  std::cout << "Rips complex is of dimension " << stree.dimension() << " - " << stree.num_simplices() << " simplices - "
            << stree.num_vertices() << " vertices." << std::endl;

  std::cout << "Iterator on Rips complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::cout << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::cout << vertex << " ";
    }
    std::cout << ") -> "
              << "[" << stree.filtration(f_simplex) << "] ";
    std::cout << std::endl;
  }

  return 0;
}
