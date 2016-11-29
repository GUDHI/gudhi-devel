/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015  INRIA (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "bottleneck distance"
#include <boost/test/unit_test.hpp>

#include <random>
#include <gudhi/Bottleneck.h>

using namespace Gudhi::bottleneck_distance;

int n1 = 81; // a natural number >0
int n2 = 180; // a natural number >0
double upper_bound = 406.43; // any real >0


std::uniform_real_distribution<double> unif(0.,upper_bound);
std::default_random_engine re;
std::vector< std::pair<double, double> > v1, v2;

BOOST_AUTO_TEST_CASE(persistence_graph){
    // Random construction
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
    Persistence_graph g(v1, v2, 0.);
    std::vector<double> d(g.sorted_distances());
    //
    BOOST_CHECK(!g.on_the_u_diagonal(n1-1));
    BOOST_CHECK(!g.on_the_u_diagonal(n1));
    BOOST_CHECK(!g.on_the_u_diagonal(n2-1));
    BOOST_CHECK(g.on_the_u_diagonal(n2));
    BOOST_CHECK(!g.on_the_v_diagonal(n1-1));
    BOOST_CHECK(g.on_the_v_diagonal(n1));
    BOOST_CHECK(g.on_the_v_diagonal(n2-1));
    BOOST_CHECK(g.on_the_v_diagonal(n2));
    //
    BOOST_CHECK(g.corresponding_point_in_u(0)==n2);
    BOOST_CHECK(g.corresponding_point_in_u(n1)==0);
    BOOST_CHECK(g.corresponding_point_in_v(0)==n1);
    BOOST_CHECK(g.corresponding_point_in_v(n2)==0);
    //
    BOOST_CHECK(g.size()==(n1+n2));
    //
    BOOST_CHECK((int) d.size() == (n1+n2)*(n1+n2) + n1 + n2 + 1);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,0))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,n1-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,n1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,n2-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,n2))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(0,(n1+n2)-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,0))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,n1-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,n1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,n2-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,n2))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance(n1,(n1+n2)-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,0))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,n1-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,n1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,n2-1))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,n2))>0);
    BOOST_CHECK(std::count(d.begin(), d.end(), g.distance((n1+n2)-1,(n1+n2)-1))>0);
}

BOOST_AUTO_TEST_CASE(neighbors_finder) {
    Persistence_graph g(v1, v2, 0.);
    Neighbors_finder nf(g, 1.);
    for(int v_point_index=1; v_point_index<((n2+n1)*9/10); v_point_index+=2)
        nf.add(v_point_index);
    //
    int v_point_index_1 = nf.pull_near(n2/2);
    BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2/2,v_point_index_1)<=1.));
    std::vector<int> l = nf.pull_all_near(n2/2);
    bool v = true;
    for(auto it = l.cbegin(); it != l.cend(); ++it)
        v = v && (g.distance(n2/2,*it)>1.);
    BOOST_CHECK(v);
    int v_point_index_2 = nf.pull_near(n2/2);
    BOOST_CHECK(v_point_index_2 == -1);
}

BOOST_AUTO_TEST_CASE(layered_neighbors_finder) {
    Persistence_graph g(v1, v2, 0.);
    Layered_neighbors_finder lnf(g, 1.);
    for(int v_point_index=1; v_point_index<((n2+n1)*9/10); v_point_index+=2)
        lnf.add(v_point_index, v_point_index % 7);
    //
    int v_point_index_1 = lnf.pull_near(n2/2,6);
    BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2/2,v_point_index_1)<=1.));
    int v_point_index_2 = lnf.pull_near(n2/2,6);
    BOOST_CHECK(v_point_index_2 == -1);
    v_point_index_1 = lnf.pull_near(n2/2,0);
    BOOST_CHECK((v_point_index_1 == -1) || (g.distance(n2/2,v_point_index_1)<=1.));
    v_point_index_2 = lnf.pull_near(n2/2,0);
    BOOST_CHECK(v_point_index_2 == -1);
}

BOOST_AUTO_TEST_CASE(graph_matching) {
    Persistence_graph g(v1, v2, 0.);
    Graph_matching m1(g);
    m1.set_r(0.);
    int e  = 0;
    while (m1.multi_augment())
        ++e;
    BOOST_CHECK(e > 0);
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

BOOST_AUTO_TEST_CASE(global){
    std::uniform_real_distribution<double> unif1(0.,upper_bound);
    std::uniform_real_distribution<double> unif2(upper_bound/10000.,upper_bound/100.);
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
    BOOST_CHECK(compute(v1, v2) <= upper_bound/100.);
    BOOST_CHECK(compute(v1, v2, upper_bound/10000.) <= upper_bound/100. + upper_bound/10000.);
    BOOST_CHECK(std::abs(compute(v1, v2) - compute(v1, v2, upper_bound/10000.)) <=  upper_bound/10000.);
}
