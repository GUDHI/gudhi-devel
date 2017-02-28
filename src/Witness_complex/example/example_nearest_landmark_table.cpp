#define BOOST_PARAMETER_MAX_ARITY 12


#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistent_cohomology.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <string>
#include <vector>
#include <list>

int main(int argc, char * const argv[]) {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
  using Witness_complex = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
  using Simplex_tree = Gudhi::Simplex_tree<>;

  Simplex_tree simplex_tree;
  Nearest_landmark_table nlt;

  // Example contains 5 witnesses and 5 landmarks
  Nearest_landmark_range w0 = {std::make_pair(0, 0), std::make_pair(1, 1), std::make_pair(2, 2),
                               std::make_pair(3, 3), std::make_pair(4, 4)}; nlt.push_back(w0);
  Nearest_landmark_range w1 = {std::make_pair(1, 0), std::make_pair(2, 1), std::make_pair(3, 2),
                               std::make_pair(4, 3), std::make_pair(0, 4)}; nlt.push_back(w1);
  Nearest_landmark_range w2 = {std::make_pair(2, 0), std::make_pair(3, 1), std::make_pair(4, 2),
                               std::make_pair(0, 3), std::make_pair(1, 4)}; nlt.push_back(w2);
  Nearest_landmark_range w3 = {std::make_pair(3, 0), std::make_pair(4, 1), std::make_pair(0, 2),
                               std::make_pair(1, 3), std::make_pair(2, 4)}; nlt.push_back(w3);
  Nearest_landmark_range w4 = {std::make_pair(4, 0), std::make_pair(0, 1), std::make_pair(1, 2),
                               std::make_pair(2, 3), std::make_pair(3, 4)}; nlt.push_back(w4);

  Witness_complex witness_complex(nlt);
  witness_complex.create_complex(simplex_tree, 4.1);

  std::cout << "Number of simplices: " << simplex_tree.num_simplices() << std::endl;

  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp > pcoh(simplex_tree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(11);

  pcoh.compute_persistent_cohomology(-0.1);
  pcoh.output_diagram();
}
