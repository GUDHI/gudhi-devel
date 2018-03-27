#include <iostream>
#include <gudhi/Filtered_toplex_map.h>
#include <gudhi/Fake_simplex_tree.h>
#include <gudhi/Sb_wrapper.h>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "toplex map"
#include <boost/test/unit_test.hpp>

using namespace Gudhi;

typedef Toplex_map::Vertex Vertex;

std::vector<Vertex> sigma1 = {1, 2, 3, 4};
std::vector<Vertex> sigma2 = {5, 2, 3, 6};
std::vector<Vertex> sigma3 = {5};
std::vector<Vertex> sigma4 = {5, 2, 3};
std::vector<Vertex> sigma5 = {5, 2, 7};
std::vector<Vertex> sigma6 = {4, 5, 3};
std::vector<Vertex> sigma7 = {4, 5, 9};
std::vector<Vertex> sigma8 = {1, 2, 3, 6};


BOOST_AUTO_TEST_CASE(toplexmap) {
    Toplex_map K;
    K.insert_simplex(sigma1);
    K.insert_simplex(sigma2);
    K.insert_simplex(sigma3);
    K.insert_simplex(sigma6);
    K.insert_simplex(sigma7);
    BOOST_CHECK(K.membership(sigma4));
    BOOST_CHECK(!K.maximality(sigma5));
    BOOST_CHECK(!K.membership(sigma5));
    K.contraction(4,5);
    BOOST_CHECK(!K.membership(sigma6));
}


BOOST_AUTO_TEST_CASE(ftoplexmap) {
    Filtered_toplex_map K;
    K.insert_simplex_and_subfaces(sigma1, 2.);
    K.insert_simplex_and_subfaces(sigma2, 2.);
    K.insert_simplex_and_subfaces(sigma6, 1.);
    K.insert_simplex_and_subfaces(sigma7, 1.);
    BOOST_CHECK(K.filtration(sigma4)==2.);
    BOOST_CHECK(K.filtration(sigma3)==1.);
}

