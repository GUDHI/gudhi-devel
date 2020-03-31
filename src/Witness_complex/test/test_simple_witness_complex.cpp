#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simple_witness_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <gudhi/Simplex_tree.h>

#include <gudhi/Witness_complex.h>

#include <iostream>
#include <vector>
#include <utility>


BOOST_AUTO_TEST_CASE(simple_witness_complex) {
  using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
  using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
  using Witness_complex = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
  using Simplex_tree = Gudhi::Simplex_tree<>;

  Simplex_tree stree;
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
  BOOST_CHECK(witness_complex.create_complex(stree, 4.1));

  std::clog << "Number of simplices: " << stree.num_simplices() << std::endl;
  BOOST_CHECK(stree.num_simplices() == 31);

  // Check when complex not empty
  BOOST_CHECK(!witness_complex.create_complex(stree, 4.1));

  // Check when max_alpha_square negative
  Simplex_tree stree2;
  BOOST_CHECK(!witness_complex.create_complex(stree2, -0.02));

  witness_complex.create_complex(stree2, 4.1, 2);
  std::clog << "Number of simplices: " << stree2.num_simplices() << std::endl;
  BOOST_CHECK(stree2.num_simplices() == 25);

}
