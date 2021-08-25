#include <iostream>
#include <string>
#include <algorithm>
#include <utility> // std::pair, std::make_pair
#include <cmath> // float comparison
#include <limits>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistent_cohomology_multi_field"
#include <boost/test/unit_test.hpp>

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Multi_field.h>

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;
using namespace boost::unit_test;

typedef Simplex_tree<> typeST;

std::string test_persistence(int min_coefficient, int max_coefficient, double min_persistence) {
  // file is copied in CMakeLists.txt
  std::ifstream simplex_tree_stream;
  simplex_tree_stream.open("simplex_tree_file_for_multi_field_unit_test.txt");
  typeST st;
  simplex_tree_stream >> st;
  simplex_tree_stream.close();

  // Display the Simplex_tree
  std::clog << "The complex contains " << st.num_simplices() << " simplices" << " - dimension= " << st.dimension()
      << std::endl;

  // Check
  BOOST_CHECK(st.num_simplices() == 58);
  BOOST_CHECK(st.dimension() == 3);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology<Simplex_tree<>, Multi_field> pcoh(st);

  pcoh.init_coefficients(min_coefficient, max_coefficient); // initializes the coefficient field for homology
  // Compute the persistent homology of the complex
  pcoh.compute_persistent_cohomology(min_persistence); // Minimal lifetime of homology feature to be recorded.

  std::ostringstream ossPers;
  pcoh.output_diagram(ossPers);

  std::string strPers = ossPers.str();
  return strPers;
}

void test_persistence_with_coeff_field(int min_coefficient, int max_coefficient) {
  // there are 2 discontinued ensembles 
  std::string value0("  0 0.25 inf");
  std::string value1("  1 0.4 inf");
  // And a big hole - cut in 2 pieces after 0.3
  std::string value2("  0 0.2 0.3");

  // For dim <= 1 =>
  std::string value3("  1 0.25 inf");
  std::string value4("  2 0.25 inf");
  std::string value5("  1 0.3 inf");
  std::string value6("  2 0.3 inf");
  std::string value7("  2 0.4 inf");

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF PERSISTENT_COHOMOLOGY_MULTI_FIELD MIN_COEFF=" << min_coefficient << " MAX_COEFF=" << max_coefficient << " MIN_PERS=0" << std::endl;

  std::string str_persistence = test_persistence(min_coefficient, max_coefficient, 0.0);
  std::clog << "str_persistence=" << str_persistence << std::endl;

  BOOST_CHECK(str_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value2) != std::string::npos); // Check found

  if ((min_coefficient < 2) && (max_coefficient < 2)) {
    BOOST_CHECK(str_persistence.find(value3) != std::string::npos); // Check found
    BOOST_CHECK(str_persistence.find(value4) != std::string::npos); // Check found
    BOOST_CHECK(str_persistence.find(value5) != std::string::npos); // Check found
    BOOST_CHECK(str_persistence.find(value6) != std::string::npos); // Check found
    BOOST_CHECK(str_persistence.find(value7) != std::string::npos); // Check found
  } else {
    BOOST_CHECK(str_persistence.find(value3) == std::string::npos); // Check not found
    BOOST_CHECK(str_persistence.find(value4) == std::string::npos); // Check not found
    BOOST_CHECK(str_persistence.find(value5) == std::string::npos); // Check not found
    BOOST_CHECK(str_persistence.find(value6) == std::string::npos); // Check not found
    BOOST_CHECK(str_persistence.find(value7) == std::string::npos); // Check not found
  }

}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_0_0) {
  test_persistence_with_coeff_field(0, 0);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_0_1) {
  test_persistence_with_coeff_field(0, 1);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_0_6) {
  test_persistence_with_coeff_field(0, 6);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_1_2) {
  test_persistence_with_coeff_field(1, 2);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_1_3) {
  test_persistence_with_coeff_field(1, 3);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_1_5) {
  test_persistence_with_coeff_field(1, 5);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_2_3) {
  test_persistence_with_coeff_field(2, 3);
}

BOOST_AUTO_TEST_CASE(persistent_cohomology_multi_field_coeff_3_4) {
  test_persistence_with_coeff_field(3, 4);
}
