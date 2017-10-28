#include <iostream>
#include <string>
#include <algorithm>
#include <utility> // std::pair, std::make_pair
#include <cmath> // float comparison
#include <limits>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "betti_numbers"
#include <boost/test/unit_test.hpp>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

struct MiniSTOptions : Gudhi::Simplex_tree_options_full_featured {
  // Implicitly use 0 as filtration value for all simplices
  static const bool store_filtration = false;
  // The persistence algorithm needs this
  static const bool store_key = true;
  // I have few vertices
  typedef short Vertex_handle;
};

using Mini_simplex_tree = Gudhi::Simplex_tree<MiniSTOptions>;
using Mini_st_persistence =
    Gudhi::persistent_cohomology::Persistent_cohomology<Mini_simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;

/*
 * Compare two intervals by dimension, then by length.
 */
template<class Simplicial_complex>
struct cmp_intervals_by_dim_then_length {
  explicit cmp_intervals_by_dim_then_length(Simplicial_complex * sc)
      : sc_(sc) { }

  template<typename Persistent_interval>
  bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
    if (sc_->dimension(get < 0 > (p1)) == sc_->dimension(get < 0 > (p2)))
      return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
              > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
    else
      return (sc_->dimension(get < 0 > (p1)) > sc_->dimension(get < 0 > (p2)));
  }
  Simplicial_complex* sc_;
};

BOOST_AUTO_TEST_CASE( plain_homology_betti_numbers )
{
  Mini_simplex_tree st;

  /* Complex to build. */
  /*    1   4          */
  /*    o---o          */
  /*   /3\ /           */
  /*  o---o   o        */
  /*  2   0   5        */
  const short tetra0123[] = {0, 1, 2, 3};
  const short edge04[] = {0, 4};
  const short edge14[] = {1, 4};
  const short vertex5[] = {5};
  st.insert_simplex_and_subfaces(tetra0123);
  st.insert_simplex_and_subfaces(edge04);
  st.insert_simplex(edge14);
  st.insert_simplex(vertex5);
  // FIXME: Remove this line
  st.set_dimension(3);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Class for homology computation
  Mini_st_persistence pcoh(st);

  // Initialize the coefficient field Z/3Z for homology
  pcoh.init_coefficients(3);

  // Compute the persistence diagram of the complex
  pcoh.compute_persistent_cohomology();

  // Print the result. The format is, on each line: 2 dim 0 inf
  // where 2 represents the field, dim the dimension of the feature.
  // 2  0 0 inf 
  // 2  0 0 inf 
  // 2  1 0 inf 
  // means that in Z/2Z-homology, the Betti numbers are b0=2 and b1=1.
  
  std::cout << "BETTI NUMBERS" << std::endl;

  BOOST_CHECK(pcoh.betti_number(0) == 2);
  BOOST_CHECK(pcoh.betti_number(1) == 1);
  BOOST_CHECK(pcoh.betti_number(2) == 0);

  std::vector<int> bns = pcoh.betti_numbers();
  BOOST_CHECK(bns.size() == 3);
  BOOST_CHECK(bns[0] == 2);
  BOOST_CHECK(bns[1] == 1);
  BOOST_CHECK(bns[2] == 0);

  std::cout << "GET PERSISTENT PAIRS" << std::endl;
  
  // Custom sort and output persistence
  cmp_intervals_by_dim_then_length<Mini_simplex_tree> cmp(&st);
  auto persistent_pairs = pcoh.get_persistent_pairs();
  
  std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
  
  BOOST_CHECK(persistent_pairs.size() == 3);
  // persistent_pairs[0] = 2  1 0 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[0])) == 1);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[0])) == 0);
  BOOST_CHECK(get<1>(persistent_pairs[0]) == st.null_simplex());

  // persistent_pairs[1] = 2  0 0 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[1])) == 0);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[1])) == 0);
  BOOST_CHECK(get<1>(persistent_pairs[1]) == st.null_simplex());

  // persistent_pairs[2] = 2  0 0 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[2])) == 0);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[2])) == 0);
  BOOST_CHECK(get<1>(persistent_pairs[2]) == st.null_simplex());

  std::cout << "INTERVALS IN DIMENSION" << std::endl;

  auto intervals_in_dimension_0 = pcoh.intervals_in_dimension(0);
  std::cout << "intervals_in_dimension_0.size() = " << intervals_in_dimension_0.size() << std::endl;
  for (std::size_t i = 0; i < intervals_in_dimension_0.size(); i++)
    std::cout << "intervals_in_dimension_0[" << i << "] = [" << intervals_in_dimension_0[i].first << "," <<
                 intervals_in_dimension_0[i].second << "]" << std::endl;
  BOOST_CHECK(intervals_in_dimension_0.size() == 2);
  BOOST_CHECK(intervals_in_dimension_0[0].first == 0);
  BOOST_CHECK(intervals_in_dimension_0[0].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());
  BOOST_CHECK(intervals_in_dimension_0[1].first == 0);
  BOOST_CHECK(intervals_in_dimension_0[1].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());


  auto intervals_in_dimension_1 = pcoh.intervals_in_dimension(1);
  std::cout << "intervals_in_dimension_1.size() = " << intervals_in_dimension_1.size() << std::endl;
  for (std::size_t i = 0; i < intervals_in_dimension_1.size(); i++)
    std::cout << "intervals_in_dimension_1[" << i << "] = [" << intervals_in_dimension_1[i].first << "," <<
                 intervals_in_dimension_1[i].second << "]" << std::endl;
  BOOST_CHECK(intervals_in_dimension_1.size() == 1);
  BOOST_CHECK(intervals_in_dimension_1[0].first == 0);
  BOOST_CHECK(intervals_in_dimension_1[0].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());

  auto intervals_in_dimension_2 = pcoh.intervals_in_dimension(2);
  std::cout << "intervals_in_dimension_2.size() = " << intervals_in_dimension_2.size() << std::endl;
  BOOST_CHECK(intervals_in_dimension_2.size() == 0);
}

using Simplex_tree = Gudhi::Simplex_tree<>;
using St_persistence =
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;

BOOST_AUTO_TEST_CASE( betti_numbers )
{
  Simplex_tree st;

  /* Complex to build. */
  /*    1   4          */
  /*    o---o          */
  /*   /3\ /           */
  /*  o---o   o        */
  /*  2   0   5        */
  const short tetra0123[] = {0, 1, 2, 3};
  const short edge04[] = {0, 4};
  const short edge14[] = {1, 4};
  const short vertex5[] = {5};
  st.insert_simplex_and_subfaces(tetra0123, 4.0);
  st.insert_simplex_and_subfaces(edge04, 2.0);
  st.insert_simplex(edge14, 2.0);
  st.insert_simplex(vertex5, 1.0);
  // FIXME: Remove this line
  st.set_dimension(3);

  // Sort the simplices in the order of the filtration
  st.initialize_filtration();

  // Class for homology computation
  St_persistence pcoh(st);

  // Initialize the coefficient field Z/3Z for homology
  pcoh.init_coefficients(3);

  // Compute the persistence diagram of the complex
  pcoh.compute_persistent_cohomology();

  // Check the Betti numbers are b0=2, b1=1 and b2=0.
  BOOST_CHECK(pcoh.betti_number(0) == 2);
  BOOST_CHECK(pcoh.betti_number(1) == 1);
  BOOST_CHECK(pcoh.betti_number(2) == 0);

  // Check the Betti numbers are b0=2, b1=1 and b2=0.
  std::vector<int> bns = pcoh.betti_numbers();
  BOOST_CHECK(bns.size() == 3);
  BOOST_CHECK(bns[0] == 2);
  BOOST_CHECK(bns[1] == 1);
  BOOST_CHECK(bns[2] == 0);
  
  // Check the persistent Betti numbers in [4., 10.] are b0=2, b1=1 and b2=0.
  BOOST_CHECK(pcoh.persistent_betti_number(0, 4., 10.) == 2);
  BOOST_CHECK(pcoh.persistent_betti_number(1, 4., 10.) == 1);
  BOOST_CHECK(pcoh.persistent_betti_number(2, 4., 10.) == 0);

  // Check the persistent Betti numbers in [2., 100.] are b0=2, b1=0 and b2=0.
  BOOST_CHECK(pcoh.persistent_betti_number(0, 2., 100.) == 2);
  BOOST_CHECK(pcoh.persistent_betti_number(1, 2., 100.) == 0);
  BOOST_CHECK(pcoh.persistent_betti_number(2, 2., 100.) == 0);

  // Check the persistent Betti numbers in [1., 1000.] are b0=1, b1=0 and b2=0.
  BOOST_CHECK(pcoh.persistent_betti_number(0, 1., 1000.) == 1);
  BOOST_CHECK(pcoh.persistent_betti_number(1, 1., 1000.) == 0);
  BOOST_CHECK(pcoh.persistent_betti_number(2, 1., 1000.) == 0);

  // Check the persistent Betti numbers in [.9, 1000.] are b0=0, b1=0 and b2=0.
  BOOST_CHECK(pcoh.persistent_betti_number(0, .9, 1000.) == 0);
  BOOST_CHECK(pcoh.persistent_betti_number(1, .9, 1000.) == 0);
  BOOST_CHECK(pcoh.persistent_betti_number(2, .9, 1000.) == 0);

  // Check the persistent Betti numbers in [4.1, 10000.] are b0=2, b1=1 and b2=0.
  bns = pcoh.persistent_betti_numbers(4.1, 10000.);
  BOOST_CHECK(bns[0] == 2);
  BOOST_CHECK(bns[1] == 1);
  BOOST_CHECK(bns[2] == 0);
  
  // Check the persistent Betti numbers in [2.1, 100000.] are b0=2, b1=0 and b2=0.
  bns = pcoh.persistent_betti_numbers(2.1, 100000.);
  BOOST_CHECK(bns[0] == 2);
  BOOST_CHECK(bns[1] == 0);
  BOOST_CHECK(bns[2] == 0);

  // Check the persistent Betti numbers in [1.1, 1000000.] are b0=1, b1=0 and b2=0.
  bns = pcoh.persistent_betti_numbers(1.1, 1000000.);
  BOOST_CHECK(bns[0] == 1);
  BOOST_CHECK(bns[1] == 0);
  BOOST_CHECK(bns[2] == 0);
  
  // Check the persistent Betti numbers in [.1, 10000000.] are b0=0, b1=0 and b2=0.
  bns = pcoh.persistent_betti_numbers(.1, 10000000.);
  BOOST_CHECK(bns[0] == 0);
  BOOST_CHECK(bns[1] == 0);
  BOOST_CHECK(bns[2] == 0);
  
  // Custom sort and output persistence
  cmp_intervals_by_dim_then_length<Simplex_tree> cmp(&st);
  auto persistent_pairs = pcoh.get_persistent_pairs();
  
  std::sort(std::begin(persistent_pairs), std::end(persistent_pairs), cmp);
  
  BOOST_CHECK(persistent_pairs.size() == 3);
  // persistent_pairs[0] = 2  1 4 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[0])) == 1);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[0])) == 4);
  BOOST_CHECK(get<1>(persistent_pairs[0]) == st.null_simplex());

  // persistent_pairs[1] = 2  0 2 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[1])) == 0);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[1])) == 2);
  BOOST_CHECK(get<1>(persistent_pairs[1]) == st.null_simplex());

  // persistent_pairs[2] = 2  0 1 inf
  BOOST_CHECK(st.dimension(get<0>(persistent_pairs[2])) == 0);
  BOOST_CHECK(st.filtration(get<0>(persistent_pairs[2])) == 1);
  BOOST_CHECK(get<1>(persistent_pairs[2]) == st.null_simplex());

  std::cout << "INTERVALS IN DIMENSION" << std::endl;

  auto intervals_in_dimension_0 = pcoh.intervals_in_dimension(0);
  std::cout << "intervals_in_dimension_0.size() = " << intervals_in_dimension_0.size() << std::endl;
  for (std::size_t i = 0; i < intervals_in_dimension_0.size(); i++)
    std::cout << "intervals_in_dimension_0[" << i << "] = [" << intervals_in_dimension_0[i].first << "," <<
                 intervals_in_dimension_0[i].second << "]" << std::endl;
  BOOST_CHECK(intervals_in_dimension_0.size() == 2);
  BOOST_CHECK(intervals_in_dimension_0[0].first == 2);
  BOOST_CHECK(intervals_in_dimension_0[0].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());
  BOOST_CHECK(intervals_in_dimension_0[1].first == 1);
  BOOST_CHECK(intervals_in_dimension_0[1].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());

  auto intervals_in_dimension_1 = pcoh.intervals_in_dimension(1);
  std::cout << "intervals_in_dimension_1.size() = " << intervals_in_dimension_1.size() << std::endl;
  for (std::size_t i = 0; i < intervals_in_dimension_1.size(); i++)
    std::cout << "intervals_in_dimension_1[" << i << "] = [" << intervals_in_dimension_1[i].first << "," <<
                 intervals_in_dimension_1[i].second << "]" << std::endl;
  BOOST_CHECK(intervals_in_dimension_1.size() == 1);
  BOOST_CHECK(intervals_in_dimension_1[0].first == 4);
  BOOST_CHECK(intervals_in_dimension_1[0].second == std::numeric_limits<Mini_simplex_tree::Filtration_value>::infinity());

  auto intervals_in_dimension_2 = pcoh.intervals_in_dimension(2);
  std::cout << "intervals_in_dimension_2.size() = " << intervals_in_dimension_2.size() << std::endl;
  BOOST_CHECK(intervals_in_dimension_2.size() == 0);
}
