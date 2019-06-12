#define BOOST_PARAMETER_MAX_ARITY 12

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>

int main(int argc, char * const argv[]) {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
  using Witness_complex = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

  Simplex_tree simplex_tree;

  // Example contains 5 witnesses and 5 landmarks
  Nearest_landmark_table nlt = {
    {{0, 0.0}, {1, 0.1}, {2, 0.2}, {3, 0.3}, {4, 0.4}}, // witness 0
    {{1, 0.0}, {2, 0.1}, {3, 0.2}, {4, 0.3}, {0, 0.4}}, // witness 1
    {{2, 0.0}, {3, 0.1}, {4, 0.2}, {0, 0.3}, {1, 0.4}}, // witness 2
    {{3, 0.0}, {4, 0.1}, {0, 0.2}, {1, 0.3}, {2, 0.4}}, // witness 3
    {{4, 0.0}, {0, 0.1}, {1, 0.2}, {2, 0.3}, {3, 0.4}}  // witness 4
  };
  /* distance(witness3, landmark3) is 0, distance(witness3, landmark4) is 0.1, etc.  */

  Witness_complex witness_complex(nlt);
  witness_complex.create_complex(simplex_tree, .41);

  std::cout << "Number of simplices: " << simplex_tree.num_simplices() << std::endl;

  Persistent_cohomology pcoh(simplex_tree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(11);

  pcoh.compute_persistent_cohomology(-0.1);
  pcoh.output_diagram();
}
