#define BOOST_TEST_MODULE bottleneck test

#include <boost/test/included/unit_test.hpp>
#include <random>
#include "gudhi/Graph_matching.h"

using namespace Gudhi::bottleneck;

Persistence_diagrams_graph* random_graph_generator(){
    int n1 = 100;
    int n2 = 120;
    double upper_bound = 80;
    // Random construction
    std::uniform_real_distribution<double> unif(0.,upper_bound);
    std::default_random_engine re;
    std::vector< std::pair<double, double> > v1, v2;
    for (int i = 0; i < n1; i++) {
        double a = unif(re);
        double b = unif(re);
        v1.emplace_back(std::min(a,b), std::max(a,b));
    }
    for (int i = 0; i < n2; i++) {
        double a = unif(re);
        double b = unif(re);
        v2.emplace_back(std::min(a,b), std::max(a,b));
    }
    return new Persistence_diagrams_graph(v1, v2, 0.);
}

BOOST_AUTO_TEST_CASE(global){
    int n = 100;
    // Random construction
    std::vector< std::pair<double, double> > v1, v2;
    for (int i = 0; i < n; i++) {
        int a = rand() % n;
        v1.emplace_back(a, a + rand() % (n - a));
        int b = rand() % n;
        v2.emplace_back(b, b + rand() % (n - b));
    }
    //
    BOOST_CHECK(bottleneck_distance(v1, v2, 1.) == 98);
}


BOOST_AUTO_TEST_CASE(persistence_diagrams_graph) {
    Persistence_diagrams_graph* g = random_graph_generator();
    std::vector<double>* d = g->sorted_distances();
    //
    BOOST_CHECK(!g->on_the_u_diagonal(99));
    BOOST_CHECK(!g->on_the_u_diagonal(100));
    BOOST_CHECK(!g->on_the_u_diagonal(119));
    BOOST_CHECK(g->on_the_u_diagonal(120));
    BOOST_CHECK(!g->on_the_v_diagonal(99));
    BOOST_CHECK(g->on_the_v_diagonal(100));
    BOOST_CHECK(g->on_the_v_diagonal(119));
    BOOST_CHECK(g->on_the_v_diagonal(120));
    //
    BOOST_CHECK(g->corresponding_point_in_u(0)==120);
    BOOST_CHECK(g->corresponding_point_in_u(100)==0);
    BOOST_CHECK(g->corresponding_point_in_v(0)==100);
    BOOST_CHECK(g->corresponding_point_in_v(120)==0);
    //
    BOOST_CHECK(g->size()==220);
    //
    BOOST_CHECK((int) d->size() <= 220*220 - 100*120 + 1); // could be more strict
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,99))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,100))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,119))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,120))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,219))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,99))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,100))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,119))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,120))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(100,219))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,99))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,100))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,119))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,120))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(219,219))==1);
    //
    delete g;
    delete d;
}

BOOST_AUTO_TEST_CASE(planar_nf) {
    Persistence_diagrams_graph* g = random_graph_generator();
    Planar_neighbors_finder pnf = Planar_neighbors_finder(*g,1.);
    for(int v_point_index=0; v_point_index<100; v_point_index+=2)
        pnf.add(v_point_index);
    pnf.remove(0);
    pnf.remove(1);
    //
    BOOST_CHECK(!pnf.contains(0));
    BOOST_CHECK(!pnf.contains(1));
    BOOST_CHECK(pnf.contains(2));
    BOOST_CHECK(!pnf.contains(3));
    //
    int v_point_index_1 = pnf.pull_near(120/2);
    BOOST_CHECK((v_point_index_1 == -1) || (g->distance(120/2,v_point_index_1)<1.));
    BOOST_CHECK(!pnf.contains(v_point_index_1));
    int v_point_index_2 = pnf.pull_near(120/2);
    BOOST_CHECK((v_point_index_2 == -1) || (g->distance(120/2,v_point_index_2)<1.));
    BOOST_CHECK(!pnf.contains(v_point_index_2));
    BOOST_CHECK((v_point_index_2 != -1) || (v_point_index_1 == -1));
    pnf.add(v_point_index_1);
    BOOST_CHECK(pnf.contains(v_point_index_1));
}

