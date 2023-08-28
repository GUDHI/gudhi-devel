/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Clock.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistence_on_rectangle.h>

#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <functional>
#include <limits>

// Set to true to test on a simple example with no finite interval.
const bool monotone = false;

int main() {
  std::vector<unsigned> sizes {1000, 999};
  std::vector<double> data(sizes[0] * sizes[1]);
  if (monotone) {
    std::iota(data.begin(), data.end(), std::size_t(0));
  } else {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0., 1.);
    std::generate(data.begin(), data.end(), std::bind(dist, gen));
  }

  Gudhi::Clock clock;
#ifndef ONLY_2D
  Gudhi::Clock clock_old;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Base> Cubical;
  Cubical complex_from_top_cells(sizes, data, true);
  std::clog << "Construction from top cells: " << clock;

  clock.begin();
  complex_from_top_cells.initialize_filtration();
  std::clog << "initialize_filtration: " << clock;

  clock.begin();
  using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
  Gudhi::persistent_cohomology::Persistent_cohomology<Cubical, Field_Zp> pers(complex_from_top_cells);
  pers.init_coefficients(2);
  pers.compute_persistent_cohomology();
  std::clog << "Compute persistent homology: " << clock;
  std::vector<std::pair<double, double>> res1;
  for (auto p: pers.get_persistent_pairs()){
    res1.emplace_back(complex_from_top_cells.filtration(std::get<0>(p)), complex_from_top_cells.filtration(std::get<1>(p)));
  }
  std::clog << "Total old code: " << clock_old << std::endl;
#endif

  clock.begin();
  std::vector<std::pair<double, double>> res2; res2.reserve(data.size() / 2);
  auto out = [&res2](double b, double d) { if (b < d) res2.emplace_back(b, d); };
  double global_min = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells(data.data(), sizes[1], sizes[0], out, out);
  res2.emplace_back(global_min, std::numeric_limits<double>::infinity());
  std::clog << "Total new code: " << clock << std::endl;

  clock.begin();
  std::vector<std::pair<double, double>> res3; res3.reserve(data.size() / 2);
  auto outi = [&res3, &data](double b, double d) { if (data[b] < data[d]) res3.emplace_back(data[b], data[d]); };
  std::size_t gm = Gudhi::cubical_complex::persistence_on_rectangle_from_top_cells<true>(data.data(), sizes[1], sizes[0], outi, outi);
  res3.emplace_back(data[gm], std::numeric_limits<double>::infinity());
  std::clog << "Total new code with index: " << clock << std::endl;

#ifndef ONLY_2D
  std::sort(res1.begin(), res1.end());
  std::sort(res2.begin(), res2.end());
  std::sort(res3.begin(), res3.end());
  if(res1 != res2) {
    std::cerr << "Bug 2!\n";
    std::exit(-2);
  }
  if(res1 != res3) {
    std::cerr << "Bug 3!\n";
    std::exit(-3);
  }
#endif

  return 0;
}
