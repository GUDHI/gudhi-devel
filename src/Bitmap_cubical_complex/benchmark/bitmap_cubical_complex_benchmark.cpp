/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hind Montassif
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Clock.h>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());

double get_random()
{
    std::uniform_real_distribution<double> dist(0., 1.);
    return dist(gen);
}

int main() {

  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;

  std::vector<unsigned> sizes_1d (1, 100000);
  std::vector<double> data_1d;
  for (size_t i = 0; i != sizes_1d[0]; ++i) {
    data_1d.push_back(get_random());
  }

  std::vector<unsigned> sizes_5d_top_cells(5, 10);
  std::vector<unsigned> sizes_5d_vertices(5, 11);
  const unsigned multiplier_top_cells = round(pow(10., 5.));
  const unsigned multiplier_vertices = round(pow(11., 5.));
  std::vector<double> data_5d_top_cells, data_5d_vertices;
  for (size_t i = 0; i != multiplier_top_cells; ++i) {
    data_5d_top_cells.push_back(get_random());
  }
  for (size_t i = 0; i != multiplier_vertices; ++i) {
    data_5d_vertices.push_back(get_random());
  }

  Gudhi::Clock cub_1d_from_top_cells_creation_clock("Cubical complex creation from 100 000 top cells in 1D");
  cub_1d_from_top_cells_creation_clock.begin();
  Bitmap_cubical_complex complex_from_top_cells_1d(sizes_1d, data_1d, true);
  cub_1d_from_top_cells_creation_clock.end();
  std::clog << cub_1d_from_top_cells_creation_clock << std::endl;

  Gudhi::Clock cub_1d_from_vertices_creation_clock("Cubical complex creation from 100 000 vertices in 1D");
  cub_1d_from_vertices_creation_clock.begin();
  Bitmap_cubical_complex complex_from_vertices_1d(sizes_1d, data_1d, false);
  cub_1d_from_vertices_creation_clock.end();
  std::clog << cub_1d_from_vertices_creation_clock << std::endl;

  Gudhi::Clock cub_5d_from_top_cells_creation_clock("Cubical complex creation from 100 000 top cells in 5D");
  cub_5d_from_top_cells_creation_clock.begin();
  Bitmap_cubical_complex complex_from_top_cells_5d(sizes_5d_top_cells, data_5d_top_cells, true);
  cub_5d_from_top_cells_creation_clock.end();
  std::clog << cub_5d_from_top_cells_creation_clock << std::endl;

  Gudhi::Clock cub_5d_from_vertices_creation_clock("Cubical complex creation from 161 051 vertices (equivalent to 100 000 top cells) in 5D");
  cub_5d_from_vertices_creation_clock.begin();
  Bitmap_cubical_complex complex_from_vertices_5d(sizes_5d_vertices, data_5d_vertices, false);
  cub_5d_from_vertices_creation_clock.end();
  std::clog << cub_5d_from_vertices_creation_clock << std::endl;

  return 0;
}
