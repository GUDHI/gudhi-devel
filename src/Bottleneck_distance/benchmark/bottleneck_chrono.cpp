/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - 2026/04 Vincent Rouvreau: Use Gudhi::random::get_default_random() in place of c++ custom use
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bottleneck.h>
#include <gudhi/Random.h>

#include <chrono>
#include <fstream>

using namespace Gudhi::persistence_diagram;


double upper_bound = 400.;  // any real > 0

int main() {
  std::ofstream result_file;
  result_file.open("results.csv", std::ios::out);

  auto rng = Gudhi::random::get_default_random();
  double delta_min = upper_bound / 1000.;
  double delta_max = upper_bound / 100.;
  for (int n = 1000; n <= 10000; n += 1000) {
    std::vector< std::pair<double, double> > v1, v2;
    for (int i = 0; i < n; i++) {
      double a = rng.get<double>(0., upper_bound);
      double b = rng.get<double>(0., upper_bound);
      double x = rng.get<double>(delta_min, delta_max);
      double y = rng.get<double>(delta_min, delta_max);
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
