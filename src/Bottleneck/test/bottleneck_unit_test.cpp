#define BOOST_TEST_MODULE bottleneck test

#include <boost/test/included/unit_test.hpp>
#include <random>
#include "../include/gudhi/Graph_matching.h"

using namespace Gudhi::bottleneck;

int n1 = 81; // a natural number >0
int n2 = 180; // a natural number >0
double upper_bound = 400.5; // any real >0

std::unique_ptr<Persistence_diagrams_graph> random_graph_generator(){
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
    return std::unique_ptr<Persistence_diagrams_graph>(new Persistence_diagrams_graph(v1, v2, 0.));
}


BOOST_AUTO_TEST_CASE(global){
    std::uniform_real_distribution<double> unif1(0.,upper_bound);
    std::uniform_real_distribution<double> unif2(upper_bound/1000.,upper_bound/100.);
    std::default_random_engine re;
    std::vector< std::pair<double, double> > v1, v2;
    for (int i = 0; i < n1; i++) {
        double a = unif1(re);
        double b = unif1(re);
        double x = unif2(re);
        double y = unif2(re);
        v1.emplace_back(std::min(a,b), std::max(a,b));
        v2.emplace_back(std::min(a,b)+std::min(x,y), std::max(a,b)+std::max(x,y));
        if(i%5==0)
            v1.emplace_back(std::min(a,b),std::min(a,b)+x);
        if(i%3==0)
            v2.emplace_back(std::max(a,b),std::max(a,b)+y);
    }
    BOOST_CHECK(bottleneck_distance(v1, v2) <= upper_bound/100.);
}


BOOST_AUTO_TEST_CASE(persistence_diagrams_graph) {
    std::unique_ptr<Persistence_diagrams_graph> g = std::move(random_graph_generator());
    std::unique_ptr< std::vector<double> > d = std::move(g->sorted_distances());
    //
    BOOST_CHECK(!g->on_the_u_diagonal(n1-1));
    BOOST_CHECK(!g->on_the_u_diagonal(n1));
    BOOST_CHECK(!g->on_the_u_diagonal(n2-1));
    BOOST_CHECK(g->on_the_u_diagonal(n2));
    BOOST_CHECK(!g->on_the_v_diagonal(n1-1));
    BOOST_CHECK(g->on_the_v_diagonal(n1));
    BOOST_CHECK(g->on_the_v_diagonal(n2-1));
    BOOST_CHECK(g->on_the_v_diagonal(n2));
    //
    BOOST_CHECK(g->corresponding_point_in_u(0)==n2);
    BOOST_CHECK(g->corresponding_point_in_u(n1)==0);
    BOOST_CHECK(g->corresponding_point_in_v(0)==n1);
    BOOST_CHECK(g->corresponding_point_in_v(n2)==0);
    //
    BOOST_CHECK(g->size()==(n1+n2));
    //
    BOOST_CHECK((int) d->size() <= (n1+n2)*(n1+n2) - n1*n2 + 1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,n1-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,n1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,n2-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,n2))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(0,(n1+n2)-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,n1-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,n1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,n2-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,n2))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance(n1,(n1+n2)-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,0))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,n1-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,n1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,n2-1))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,n2))==1);
    BOOST_CHECK(std::count(d->begin(), d->end(), g->distance((n1+n2)-1,(n1+n2)-1))==1);
}

BOOST_AUTO_TEST_CASE(planar_neighbors_finder) {
    std::unique_ptr<Persistence_diagrams_graph> g = std::move(random_graph_generator());
    Planar_neighbors_finder pnf = Planar_neighbors_finder(*g,1.);
    for(int v_point_index=0; v_point_index<n1; v_point_index+=2)
        pnf.add(v_point_index);
    //
    BOOST_CHECK(pnf.contains(0));
    BOOST_CHECK(!pnf.contains(1));
    BOOST_CHECK(pnf.contains(2));
    BOOST_CHECK(!pnf.contains(3));
    //
    pnf.remove(0);
    pnf.remove(1);
    //
    BOOST_CHECK(!pnf.contains(0));
    BOOST_CHECK(!pnf.contains(1));
    BOOST_CHECK(pnf.contains(2));
    BOOST_CHECK(!pnf.contains(3));
    //
    int v_point_index_1 = pnf.pull_near(n2/2);
    BOOST_CHECK((v_point_index_1 == -1) || (g->distance(n2/2,v_point_index_1)<=1.));
    BOOST_CHECK(!pnf.contains(v_point_index_1));
    std::list<int> l = *pnf.pull_all_near(n2/2);
    bool v = true;
    for(auto it = l.cbegin(); it != l.cend(); ++it)
        v = v && (g->distance(n2/2,*it)>1.);
    BOOST_CHECK(v);
    int v_point_index_2 = pnf.pull_near(n2/2);
    BOOST_CHECK(v_point_index_2 == -1);
    pnf.add(v_point_index_1);
    BOOST_CHECK(pnf.contains(v_point_index_1));
}


BOOST_AUTO_TEST_CASE(neighbors_finder) {
    std::unique_ptr<Persistence_diagrams_graph> g = std::move(random_graph_generator());
    Neighbors_finder nf = Neighbors_finder(*g,1.);
    for(int v_point_index=1; v_point_index<((n2+n1)*9/10); v_point_index+=2)
        nf.add(v_point_index);
    //
    int v_point_index_1 = nf.pull_near(n2/2);
    BOOST_CHECK((v_point_index_1 == -1) || (g->distance(n2/2,v_point_index_1)<=1.));
    std::list<int> l = *nf.pull_all_near(n2/2);
    bool v = true;
    for(auto it = l.cbegin(); it != l.cend(); ++it)
        v = v && (g->distance(n2/2,*it)>1.);
    BOOST_CHECK(v);
    int v_point_index_2 = nf.pull_near(n2/2);
    BOOST_CHECK(v_point_index_2 == -1);
}

BOOST_AUTO_TEST_CASE(layered_neighbors_finder) {
    std::unique_ptr<Persistence_diagrams_graph> g = std::move(random_graph_generator());
    Layered_neighbors_finder lnf = Layered_neighbors_finder(*g,1.);
    for(int v_point_index=1; v_point_index<((n2+n1)*9/10); v_point_index+=2)
        lnf.add(v_point_index, v_point_index % 7);
    //
    int v_point_index_1 = lnf.pull_near(n2/2,6);
    BOOST_CHECK((v_point_index_1 == -1) || (g->distance(n2/2,v_point_index_1)<=1.));
    int v_point_index_2 = lnf.pull_near(n2/2,6);
    BOOST_CHECK(v_point_index_2 == -1);
    v_point_index_1 = lnf.pull_near(n2/2,0);
    BOOST_CHECK((v_point_index_1 == -1) || (g->distance(n2/2,v_point_index_1)<=1.));
    v_point_index_2 = lnf.pull_near(n2/2,0);
    BOOST_CHECK(v_point_index_2 == -1);
}


BOOST_AUTO_TEST_CASE(graph_matching) {
    std::unique_ptr<Persistence_diagrams_graph> g = std::move(random_graph_generator());
    Graph_matching m1(*g);
    m1.set_r(0.);
    int e  = 0;
    while (m1.multi_augment())
        ++e;
    BOOST_CHECK(e <= 2*sqrt(2*(n1+n2)));
    Graph_matching m2 = m1;
    BOOST_CHECK(!m2.multi_augment());
    m2.set_r(upper_bound);
    e  = 0;
    while (m2.multi_augment())
        ++e;
    BOOST_CHECK(e <= 2*sqrt(2*(n1+n2)));
    BOOST_CHECK(m2.perfect());
    BOOST_CHECK(!m1.perfect());
}
