#include <iostream>
#include <vector>
#include <gudhi/Lazy_Toplex_map.h>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "toplex map"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

using list_of_tested_variants = boost::mpl::list<Gudhi::Toplex_map, Gudhi::Lazy_Toplex_map>;

BOOST_AUTO_TEST_CASE_TEMPLATE(common_toplex_map_functionnalities, Toplex_map, list_of_tested_variants) {
  using Vertex = typename Toplex_map::Vertex;

  std::vector<Vertex> sigma1 = {1, 2, 3, 4};
  std::vector<Vertex> sigma2 = {5, 2, 3, 6};
  std::vector<Vertex> sigma3 = {5};
  std::vector<Vertex> sigma4 = {5, 2, 3};
  std::vector<Vertex> sigma5 = {5, 2, 7};
  std::vector<Vertex> sigma6 = {4, 5, 3};
  std::vector<Vertex> sigma7 = {4, 5, 9};
  std::vector<Vertex> sigma8 = {1, 2, 3, 6};

  Toplex_map K;
  K.insert_simplex(sigma1);
  K.insert_simplex(sigma2);
  K.insert_simplex(sigma3);
  K.insert_simplex(sigma6);
  K.insert_simplex(sigma7);

  std::cout << K.num_maximal_simplices();

  BOOST_CHECK(K.membership(sigma4));
  //BOOST_CHECK(!K.maximality(sigma3));
  BOOST_CHECK(!K.membership(sigma5));
  K.insert_simplex(sigma5);

  std::cout << K.num_maximal_simplices();

  BOOST_CHECK(K.membership(sigma5));
  std::vector<Vertex> sigma9 = {1, 2, 3};
  std::vector<Vertex> sigma10 = {2, 7};
  auto r = K.contraction(4,5);

  std::cout << K.num_maximal_simplices();

  sigma9.emplace_back(r);
  sigma10.emplace_back(r);
  BOOST_CHECK(!K.membership(sigma6));
  BOOST_CHECK(K.membership(sigma9));
  BOOST_CHECK(K.membership(sigma10));
  K.remove_simplex(sigma10);
  BOOST_CHECK(!K.membership(sigma10));

}

BOOST_AUTO_TEST_CASE(toplex_map_num_maximal_simplices) {
  using Vertex = Gudhi::Toplex_map::Vertex;

  Gudhi::Toplex_map K;
  K.insert_simplex({1, 2, 3, 4});
  K.insert_simplex({5, 2, 3, 6});
  K.insert_simplex({4, 5, 3});
  K.insert_simplex({4, 5, 9});

  std::cout << K.num_maximal_simplices();
  BOOST_CHECK(K.num_maximal_simplices() == 4);

  K.insert_simplex({5, 2, 7});

  std::cout << K.num_maximal_simplices();
  BOOST_CHECK(K.num_maximal_simplices() == 5);

  auto r = K.contraction(4,5);

  std::cout << K.num_maximal_simplices();
  BOOST_CHECK(K.num_maximal_simplices() == 4);

}

BOOST_AUTO_TEST_CASE(toplex_map_maximality) {
  using Vertex = Gudhi::Toplex_map::Vertex;

  std::vector<Vertex> sigma1 = {1, 2, 3, 4};
  std::vector<Vertex> sigma2 = {5, 2, 3, 6};
  std::vector<Vertex> sigma3 = {5};
  std::vector<Vertex> sigma4 = {4, 5, 3};
  std::vector<Vertex> sigma5 = {4, 5, 9};

  Gudhi::Toplex_map K;
  K.insert_simplex(sigma1);
  K.insert_simplex(sigma2);
  K.insert_simplex(sigma3);
  K.insert_simplex(sigma4);
  K.insert_simplex(sigma5);
  BOOST_CHECK(K.maximality(sigma1));
  BOOST_CHECK(K.maximality(sigma2));
  BOOST_CHECK(!K.maximality(sigma3));
  BOOST_CHECK(K.maximality(sigma4));
  BOOST_CHECK(K.maximality(sigma5));
}

