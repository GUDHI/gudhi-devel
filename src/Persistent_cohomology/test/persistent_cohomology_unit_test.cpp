#include <iostream>
#include <string>
#include <algorithm>
#include <utility> // std::pair, std::make_pair
#include <cmath> // float comparison
#include <limits>
#include <cstdint>  // for std::uint8_t

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "persistent_cohomology"
#include <boost/test/unit_test.hpp>

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;
using namespace boost::unit_test;

typedef Simplex_tree<> typeST;

std::string test_persistence(int coefficient, int min_persistence) {
  // file is copied in CMakeLists.txt
  std::ifstream simplex_tree_stream;
  simplex_tree_stream.open("simplex_tree_file_for_unit_test.txt");
  typeST st;
  simplex_tree_stream >> st;
  simplex_tree_stream.close();

  // Display the Simplex_tree
  std::clog << "The complex contains " << st.num_simplices() << " simplices" << " - dimension= " << st.dimension()
      << std::endl;

  // Check
  BOOST_CHECK(st.num_simplices() == 98);
  BOOST_CHECK(st.dimension() == 3);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Compute the persistence diagram of the complex
  Persistent_cohomology<Simplex_tree<>, Field_Zp> pcoh(st);

  pcoh.init_coefficients( coefficient );  // initializes the coefficient field for homology
  // Compute the persistent homology of the complex
  pcoh.compute_persistent_cohomology( min_persistence );  // Minimal lifetime of homology feature to be recorded.
  std::ostringstream ossPers;

  pcoh.output_diagram(ossPers);
  std::string strPers = ossPers.str();
  return strPers;
}

void test_persistence_with_coeff_field(int coeff_field) {
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

  value0.insert(0,std::to_string(coeff_field));
  value1.insert(0,std::to_string(coeff_field));
  value2.insert(0,std::to_string(coeff_field));
  value3.insert(0,std::to_string(coeff_field));
  value4.insert(0,std::to_string(coeff_field));
  value5.insert(0,std::to_string(coeff_field));
  value6.insert(0,std::to_string(coeff_field));
  value7.insert(0,std::to_string(coeff_field));
  value8.insert(0,std::to_string(coeff_field));
  value9.insert(0,std::to_string(coeff_field));

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF PERSISTENT_COHOMOLOGY_SINGLE_FIELD COEFF_FIELD=" << coeff_field << " MIN_PERS=0" << std::endl;

  std::string str_persistence = test_persistence(coeff_field, 0);
  std::clog << str_persistence << std::endl;
  
  BOOST_CHECK(str_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value3) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value4) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value5) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value6) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value7) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value9) != std::string::npos); // Check found
  std::clog << "str_persistence=" << str_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF PERSISTENT_COHOMOLOGY_SINGLE_FIELD COEFF_FIELD=" << coeff_field << " MIN_PERS=1" << std::endl;

  str_persistence = test_persistence(coeff_field, 1);

  BOOST_CHECK(str_persistence.find(value0) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value1) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value2) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value3) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value4) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value5) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value6) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value7) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value9) != std::string::npos); // Check found
  std::clog << "str_persistence=" << str_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF PERSISTENT_COHOMOLOGY_SINGLE_FIELD COEFF_FIELD=" << coeff_field << " MIN_PERS=2" << std::endl;

  str_persistence = test_persistence(coeff_field, 2);

  BOOST_CHECK(str_persistence.find(value0) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value1) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value2) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value3) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value4) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value5) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value6) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value7) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value9) != std::string::npos); // Check found
  std::clog << "str_persistence=" << str_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF PERSISTENT_COHOMOLOGY_SINGLE_FIELD COEFF_FIELD=" << coeff_field << " MIN_PERS=Inf" << std::endl;

  str_persistence = test_persistence(coeff_field, (std::numeric_limits<int>::max)());

  BOOST_CHECK(str_persistence.find(value0) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value1) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value2) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value3) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value4) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value5) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value6) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value7) == std::string::npos); // Check not found
  BOOST_CHECK(str_persistence.find(value8) != std::string::npos); // Check found
  BOOST_CHECK(str_persistence.find(value9) != std::string::npos); // Check found
  std::clog << "str_persistence=" << str_persistence << std::endl;
}

BOOST_AUTO_TEST_CASE( persistent_cohomology_single_field_coeff_not_prime )
{
  for (auto non_prime : {0, 1, 4, 6})
    BOOST_CHECK_THROW(test_persistence_with_coeff_field(non_prime), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( persistent_cohomology_single_field_coeff_prime )
{
  for (auto prime : {2, 3, 5, 11, 13})
    test_persistence_with_coeff_field(prime);
}

BOOST_AUTO_TEST_CASE( persistent_cohomology_single_field_coeff_limit )
{
  BOOST_CHECK_THROW(test_persistence_with_coeff_field(46349), std::invalid_argument);
}

/** SimplexTree minimal options to test the limits.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint8_t>::max()<\CODE> = 256.*/
struct MiniSTOptions {
  typedef linear_indexing_tag Indexing_tag;
  typedef short Vertex_handle;
  typedef double Filtration_value;
  // Maximum number of simplices to compute persistence is 2^8 - 1 = 255. One is reserved for null_key
  typedef std::uint8_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = false;
  static const bool contiguous_vertices = false;
};

using Mini_simplex_tree = Gudhi::Simplex_tree<MiniSTOptions>;
using Mini_st_persistence =
    Gudhi::persistent_cohomology::Persistent_cohomology<Mini_simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;

BOOST_AUTO_TEST_CASE( persistence_constructor_exception )
{
  Mini_simplex_tree st;

  // To make number of simplices = 255
  const short simplex_0[] = {0, 1, 2, 3, 4, 5, 6, 7};
  st.insert_simplex_and_subfaces(simplex_0);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  BOOST_CHECK(st.num_simplices() <= std::numeric_limits<MiniSTOptions::Simplex_key>::max());
  // Class for homology computation
  BOOST_CHECK_NO_THROW(Mini_st_persistence pcoh(st));

  st.insert_simplex({8});
  BOOST_CHECK(st.num_simplices() > std::numeric_limits<MiniSTOptions::Simplex_key>::max());
  // Class for homology computation
  BOOST_CHECK_THROW(Mini_st_persistence pcoh2(st), std::out_of_range);

}
