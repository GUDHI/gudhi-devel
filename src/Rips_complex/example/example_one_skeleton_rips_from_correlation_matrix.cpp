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
  // |1     0.06  0.23  0.01  0.89|
  // |0.06  1     0.74  0.01  0.61|
  // |0.23  0.74  1     0.72  0.03|
  // |0.01  0.01  0.72  1     0.7 |
  // |0.89  0.61  0.03  0.7   1   |

  Distance_matrix correlations;
  correlations.push_back({});
  correlations.push_back({0.06});
  correlations.push_back({0.23, 0.74});
  correlations.push_back({0.01, 0.01, 0.72});
  correlations.push_back({0.89, 0.61, 0.03, 0.7});

  // ----------------------------------------------------------------------------
  // Convert correlation matrix to a distance matrix:
  // ----------------------------------------------------------------------------
  double threshold = 0;
  for (size_t i = 0; i != correlations.size(); ++i) {
    for (size_t j = 0; j != correlations[i].size(); ++j) {
      // Here we check if our data comes from corelation matrix.
      if ((correlations[i][j] < -1) || (correlations[i][j] > 1)) {
        std::cerr << "The input matrix is not a correlation matrix. The program will now terminate.\n";
        throw "The input matrix is not a correlation matrix. The program will now terminate.\n";
      }
      correlations[i][j] = 1 - correlations[i][j];
      // Here we make sure that we will get the treshold value equal to maximal
      // distance in the matrix.
      if (correlations[i][j] > threshold) threshold = correlations[i][j];
    }
  }

  //-----------------------------------------------------------------------------
  // Now the correlation matrix is a distance matrix and can be processed further.
  //-----------------------------------------------------------------------------
  Distance_matrix distances = correlations;

  Rips_complex rips_complex_from_points(distances, threshold);

  Simplex_tree stree;
  rips_complex_from_points.create_complex(stree, 1);
  // ----------------------------------------------------------------------------
  // Display information about the one skeleton Rips complex. Note that
  // the filtration displayed here comes from the distance matrix computed
  // above, which is 1 - initial correlation matrix. Only this way, we obtain
  // a complex with filtration. If a correlation matrix is used instead, we would
  // have a reverse filtration (i.e. filtration of boundary of each simplex S
  // is greater or equal to the filtration of S).
  // ----------------------------------------------------------------------------
  std::clog << "Rips complex is of dimension " << stree.dimension() << " - " << stree.num_simplices() << " simplices - "
            << stree.num_vertices() << " vertices." << std::endl;

  std::clog << "Iterator on Rips complex simplices in the filtration order, with [filtration value]:" << std::endl;
  for (auto f_simplex : stree.filtration_simplex_range()) {
    std::clog << "   ( ";
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
      std::clog << vertex << " ";
    }
    std::clog << ") -> "
              << "[" << stree.filtration(f_simplex) << "] ";
    std::clog << std::endl;
  }

  return 0;
}
