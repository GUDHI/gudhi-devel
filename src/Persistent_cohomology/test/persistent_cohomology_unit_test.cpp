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

std::string test_rips_persistence(int coefficient, int min_persistence) {
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

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=0" << std::endl;

  std::string str_rips_persistence = test_rips_persistence(dimension, 0);
  std::clog << str_rips_persistence << std::endl;
  
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
  std::clog << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=1" << std::endl;

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
  std::clog << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=2" << std::endl;

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
  std::clog << "str_rips_persistence=" << str_rips_persistence << std::endl;

  std::clog << "********************************************************************" << std::endl;
  std::clog << "TEST OF RIPS_PERSISTENT_COHOMOLOGY_SINGLE_FIELD DIM=" << dimension << " MIN_PERS=Inf" << std::endl;

  str_rips_persistence = test_rips_persistence(dimension, (std::numeric_limits<int>::max)());

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
  std::clog << "str_rips_persistence=" << str_rips_persistence << std::endl;
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

/** SimplexTree minimal options to test the limits.
 * 
 * Maximum number of simplices to compute persistence is <CODE>std::numeric_limits<std::uint8_t>::max()<\CODE> = 256.*/
struct MiniSTOptions {
  typedef linear_indexing_tag Indexing_tag;
  static const bool is_zigzag = false;
  typedef short Vertex_handle;
  typedef double Filtration_value;
  // Maximum number of simplices to compute persistence is 2^8 - 1 = 255. One is reserved for null_key
  typedef std::uint8_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = false;
  static const bool contiguous_vertices = false;
  static const bool simplex_handle_strong_validity = false;
  static const bool link_nodes_by_label = false;
  static const bool store_morse_matching = false;
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
