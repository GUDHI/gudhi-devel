/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Sophia-Antipolis (France)
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

#include <gudhi/Graph_matching.h>
#include <iostream>

#include <chrono>
#include <fstream>

using namespace Gudhi::bipartite_graph_matching;


double upper_bound = 400.; // any real >0

int main(){
    std::ofstream objetfichier;
    objetfichier.open("results.csv", std::ios::out);

    for(int n =50; n<=1000; n+=100){
std::cout << n << "\n";
        std::uniform_real_distribution<double> unif1(0.,upper_bound);
        std::uniform_real_distribution<double> unif2(upper_bound/1000.,upper_bound/100.);
        std::default_random_engine re;
        std::vector< std::pair<double, double> > v1, v2;
        for (int i = 0; i < n; i++) {
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

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        double b = bottleneck_distance(v1,v2);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        typedef std::chrono::duration<int,std::milli> millisecs_t;
        millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end-start));
        objetfichier << n << ";" << duration.count() << ";" << b << std::endl;
    }
    objetfichier.close();
}


/*
int main() {
    std::vector< std::pair<double,double> > v1, v2;

    v1.push_back(std::pair<double,double>(2.7,3.7));
    v1.push_back(std::pair<double,double>(9.6,14));
    v1.push_back(std::pair<double,double>(34.2,34.974));

    v2.push_back(std::pair<double,double>(2.8,4.45));
    v2.push_back(std::pair<double,double>(9.5,14.1));


    double b = bottleneck_distance(v1, v2);

    std::cout << "Bottleneck distance = " << b << std::endl;

}*/
