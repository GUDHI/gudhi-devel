/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - 2026/04 Vincent Rouvreau: Use Gudhi::random in place of c++ custom use
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Clock.h>
#include <gudhi/Random.h>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>

int main() {

  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
  
  std::vector<unsigned> sizes_1d (1, 3000000);
  std::vector<double> data_1d = Gudhi::random::get_uniform_range<double>(sizes_1d[0], 0., 1.);

  std::vector<unsigned> sizes_5d_top_cells(5, 10);
  std::vector<unsigned> sizes_5d_vertices(5, 11);
  const unsigned multiplier_top_cells = round(pow(10., 5.));
  const unsigned multiplier_vertices = round(pow(11., 5.));
  std::vector<double> data_5d_top_cells = Gudhi::random::get_uniform_range<double>(multiplier_top_cells, 0., 1.);
  std::vector<double> data_5d_vertices = Gudhi::random::get_uniform_range<double>(multiplier_vertices, 0., 1.);

  Gudhi::Clock cub_1d_from_top_cells_creation_clock("Cubical complex creation from 3 000 000 top cells in 1D");
  cub_1d_from_top_cells_creation_clock.begin();
  Bitmap_cubical_complex complex_from_top_cells_1d(sizes_1d, data_1d, true);
  std::clog << cub_1d_from_top_cells_creation_clock << std::endl;

  Gudhi::Clock cub_1d_from_vertices_creation_clock("Cubical complex creation from 3 000 000 vertices in 1D");
  cub_1d_from_vertices_creation_clock.begin();
  Bitmap_cubical_complex complex_from_vertices_1d(sizes_1d, data_1d, false);
  std::clog << cub_1d_from_vertices_creation_clock << std::endl;

  Gudhi::Clock cub_5d_from_top_cells_creation_clock("Cubical complex creation from 100 000 top cells in 5D");
  cub_5d_from_top_cells_creation_clock.begin();
  Bitmap_cubical_complex complex_from_top_cells_5d(sizes_5d_top_cells, data_5d_top_cells, true);
  std::clog << cub_5d_from_top_cells_creation_clock << std::endl;

  Gudhi::Clock cub_5d_from_vertices1_creation_clock("Cubical complex creation from 100 000 vertices in 5D");
  cub_5d_from_vertices1_creation_clock.begin();
  Bitmap_cubical_complex complex_from_vertices1_5d(sizes_5d_top_cells, data_5d_top_cells, false);
  std::clog << cub_5d_from_vertices1_creation_clock << std::endl;

  Gudhi::Clock cub_5d_from_vertices2_creation_clock("Cubical complex creation from 161 051 vertices (equivalent to 100 000 top cells) in 5D");
  cub_5d_from_vertices2_creation_clock.begin();
  Bitmap_cubical_complex complex_from_vertices2_5d(sizes_5d_vertices, data_5d_vertices, false);
  std::clog << cub_5d_from_vertices2_creation_clock << std::endl;

  return 0;
}
