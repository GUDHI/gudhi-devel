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

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

typedef Simplex_tree<> typeST;

std::string test_rips_persistence(int coefficient, int min_persistence) {
  const std::string inputFile("simplex_tree_file_for_unit_test.txt");
  std::ifstream simplex_tree_stream;
  simplex_tree_stream.open(inputFile.c_str());
  typeST st;
  simplex_tree_stream >> st;
  simplex_tree_stream.close();

  // Display the Simplex_tree
  std::cout << "The complex contains " << st.num_simplices() << " simplices" << " - dimension= " << st.dimension()
      << " - filtration= " << st.filtration() << std::endl;

  // Check
  BOOST_CHECK(st.num_simplices() == 98);
  BOOST_CHECK(st.dimension() == 3);
  BOOST_CHECK(st.filtration() == 1.89);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology<Simplex_tree<>, Field_Zp> pcoh(st);

  pcoh.init_coefficients( coefficient );  // initializes the coefficient field for homology
  // Check infinite rips
  pcoh.compute_persistent_cohomology( min_persistence );  // Minimal lifetime of homology feature to be recorded.
  std::ostringstream ossInfinite;

  pcoh.output_diagram(ossInfinite);
  std::string strInfinite = ossInfinite.str();
  return strInfinite;
}

void test_rips_persistence_in_dimension(int dimension) {
  std::string value0("  0 0.02 1.12");
  std::string value1("  0 0.03 1.13");
  std::string value2("  0 0.04 1.14");
  std::string value3("  0 0.05 1.15");
  std::string value4("  0 0.06 1.16");
  std::string value5("  0 0.07 1.17");
  std::string value6("  0 0.08 1.18");
  std::string value7("  0 0.09 1.19");
  std::string value8("  0 0 inf"    );
  std::string value9("  0 0.01 inf" );

  value0.insert(0,std::to_string(dimension));
  value1.insert(0,std::to_string(dimension));
  value2.insert(0,std::to_string(dimension));
  value3.insert(0,std::to_string(dimension));
  value4.insert(0,std::to_string(dimension));
  value5.insert(0,std::to_string(dimension));
  value6.insert(0,std::to_string(dimension));
  value7.insert(0,std::to_string(dimension));
  value8.insert(0,std::to_string(dimension));
  value9.insert(0,std::to_string(dimension));

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=0" << std::endl;

  std::string str_rips_persistence = test_rips_persistence(dimension, 0);

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value4) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value5) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value6) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value7) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value9) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=1" << std::endl;

  str_rips_persistence = test_rips_persistence(dimension, 1);

  BOOST_CHECK(str_rips_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value3) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value4) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value5) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value6) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value7) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value9) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=2" << std::endl;

  str_rips_persistence = test_rips_persistence(dimension, 2);

  BOOST_CHECK(str_rips_persistence.find(value0) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value1) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value2) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value3) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value4) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value5) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value6) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value7) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value9) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::cout << "********************************************************************" << std::endl;
  std::cout << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=Inf" << std::endl;

  str_rips_persistence = test_rips_persistence(dimension, std::numeric_limits<int>::max());

  BOOST_CHECK(str_rips_persistence.find(value0) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value1) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value2) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value3) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value4) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value5) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value6) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value7) == std::string::npos); // Check not found
  BOOST_CHECK(str_rips_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_rips_persistence.find(value9) != std::string::npos); // Check found
  std::cout << "str_rips_persistence=" << str_rips_persistence << std::endl;
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_single_field_dim_1 )
{
  test_rips_persistence_in_dimension(1);
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_single_field_dim_2 )
{
  test_rips_persistence_in_dimension(2);
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_single_field_dim_3 )
{
  test_rips_persistence_in_dimension(3);
}

BOOST_AUTO_TEST_CASE( rips_persistent_cohomology_single_field_dim_5 )
{
  test_rips_persistence_in_dimension(5);
}

// TODO(VR): not working from 6
// std::string str_rips_persistence = test_rips_persistence(6, 0);
// TODO(VR): division by zero
// std::string str_rips_persistence = test_rips_persistence(0, 0);
