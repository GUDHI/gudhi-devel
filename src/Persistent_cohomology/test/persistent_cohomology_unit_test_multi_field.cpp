#define BOOST_TEST_MODULE const_string test
#include <boost/test/included/unit_test.hpp>
#include <boost/system/error_code.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <iostream>
#include <string>

#include <utility> // std::pair, std::make_pair

#include <cmath> // float comparison
#include <limits>

#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/reader_utils.h"
#include "gudhi/Simplex_tree.h"
#include "gudhi/Persistent_cohomology.h"
#include "gudhi/Persistent_cohomology/Multi_field.h"

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

typedef Simplex_tree<> typeST;

std::string test_rips_persistence(int min_coefficient, int max_coefficient, int min_persistence) {
  const std::string inputFile("simplex_tree_file_for_multi_field_unit_test.txt");
  std::ifstream simplex_tree_stream;
  simplex_tree_stream.open(inputFile.c_str());
  typeST st;
  simplex_tree_stream >> st;
  simplex_tree_stream.close();

  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << " - dimension= " << st.dimension()
      << " - filtration= " << st.filtration() << std::endl;

  // Check
  BOOST_CHECK(st.num_simplices() == 6142604);
  BOOST_CHECK(st.dimension() == 3);
  BOOST_CHECK(st.filtration() == 0.249999);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology<Simplex_tree<>, Multi_field> pcoh(st);

  pcoh.init_coefficients(min_coefficient, max_coefficient);  // initializes the coefficient field for homology
  // Check infinite rips
  pcoh.compute_persistent_cohomology(min_persistence);  // Minimal lifetime of homology feature to be recorded.

  std::ostringstream ossRips;
  pcoh.output_diagram(ossRips);

  std::string strRips = ossRips.str();
  return strRips;
}

void test_rips_persistence_in_dimension(int min_dimension, int max_dimension) {
  std::string value0("  0 0 inf");
  std::string value1("  1 0.0702103 inf");
  std::string value2("2  1 0.0702103 inf");
  std::string value3("2  2 0.159992 inf");

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_MULTI_FIELD MIN_DIM=" << min_dimension << " MAX_DIM=" << max_dimension << " MIN_PERS=0" << std::endl;

  std::string str_rips_persistence = test_rips_persistence(min_dimension, max_dimension, 1);

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_MULTI_FIELD DIM=" << min_dimension << " MAX_DIM=" << max_dimension << " MIN_PERS=2" << std::endl;

  str_rips_persistence = test_rips_persistence(min_dimension, max_dimension, 2);

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_MULTI_FIELD DIM=" << min_dimension << " MAX_DIM=" << max_dimension << " MIN_PERS=3" << std::endl;

  str_rips_persistence = test_rips_persistence(min_dimension, max_dimension, 3);

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_MULTI_FIELD DIM=" << min_dimension << " MAX_DIM=" << max_dimension << " MIN_PERS=Inf" << std::endl;

  str_rips_persistence = test_rips_persistence(min_dimension, max_dimension, std::numeric_limits<int>::max());

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_multi_field_dim_1_2 )
{
  test_rips_persistence_in_dimension(1, 2);
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_multi_field_dim_2_3 )
{
  test_rips_persistence_in_dimension(2, 3);
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_multi_field_dim_1_5 )
{
  test_rips_persistence_in_dimension(1, 5);
}

// TODO(VR): not working from 6
// std::string str_rips_persistence = test_rips_persistence(6, 0);
// TODO(VR): division by zero
// std::string str_rips_persistence = test_rips_persistence(0, 0);
// TODO(VR): is result OK of :
// test_rips_persistence_in_dimension(3, 4);

