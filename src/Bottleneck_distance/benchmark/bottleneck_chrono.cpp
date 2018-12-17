/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
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

#include <gudhi/Bottleneck.h>
#include <chrono>
#include <fstream>
#include <random>

using namespace Gudhi::persistence_diagram;


double upper_bound = 400.;  // any real > 0

int main() {
  std::ofstream result_file;
  result_file.open("results.csv", std::ios::out);

  for (int n = 1000; n <= 10000; n += 1000) {
    std::uniform_real_distribution<double> unif1(0., upper_bound);
    std::uniform_real_distribution<double> unif2(upper_bound / 1000., upper_bound / 100.);
    std::default_random_engine re;
    std::vector< std::pair<double, double> > v1, v2;
    for (int i = 0; i < n; i++) {
      double a = unif1(re);
      double b = unif1(re);
      double x = unif2(re);
      double y = unif2(re);
      v1.emplace_back(std::min(a, b), std::max(a, b));
      v2.emplace_back(std::min(a, b) + std::min(x, y), std::max(a, b) + std::max(x, y));
      if (i % 5 == 0)
        v1.emplace_back(std::min(a, b), std::min(a, b) + x);
      if (i % 3 == 0)
        v2.emplace_back(std::max(a, b), std::max(a, b) + y);
    }
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double b = bottleneck_distance(v1, v2);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    typedef std::chrono::duration<int, std::milli> millisecs_t;
    millisecs_t duration(std::chrono::duration_cast<millisecs_t>(end - start));
    result_file << n << ";" << duration.count() << ";" << b << std::endl;
  }
  result_file.close();
}
