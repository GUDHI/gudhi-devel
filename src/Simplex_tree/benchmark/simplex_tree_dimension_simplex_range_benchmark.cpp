/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Clock.h>

#include <iostream>
#include <random>
#include <numeric>  // for std::iota
#include <algorithm>  // for std::shuffle and std::sample
#include <vector>
#include <cstdlib>

std::random_device rd;

template <class Vertex_handle>
std::vector<Vertex_handle> rand_int_range(int subset_min_size,
                                          int subset_max_size,
                                          Vertex_handle range_max_value)
{
#ifdef DEBUG_TRACES
  std::clog << "rand_int_range - subset_min_size = " << subset_min_size << " - subset_max_size = " << subset_max_size
            << " - range_max_value = " << range_max_value << std::endl;
#endif  // DEBUG_TRACES

  std::vector<Vertex_handle> range(range_max_value);
  std::iota(range.begin(), range.end(), 0); // range is { 0, 1, 2, ..., 99 } when range_max_value is 100

  std::shuffle(range.begin(), range.end(), std::mt19937 { rd() });

  std::uniform_int_distribution<int> dist(subset_min_size, subset_max_size);
  // Return a subset, which size is in between [subset_min_size; subset_max_size], of shuffled range
  // {1, 9, 13, 19, 86, 36}, for example
  return std::vector<Vertex_handle> {range.begin(), range.begin() + dist(rd)};
}

int main(int argc, char *argv[]) {
  int nb_vertices = 100000;
  if (argc > 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct\n";
    std::cerr << "Usage: " << (argv[0] - 1) << " [NB_VERTICES] \n";
    std::cerr << "    NB_VERTICES is 5.000 by default.\n";
    exit(EXIT_FAILURE);  // ----- >>
  }
  if (argc == 2)
    nb_vertices = atoi(argv[1]);

  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Vertex_handle  = typename Simplex_tree::Vertex_handle;

  std::clog << "... Simplex_tree construction\n";

  Simplex_tree st;
  // Insert 'nb_vertices' random simplices, of size in between [2; 5] and vertices in between [0; nb_vertices]
  for (Vertex_handle v=0; v < nb_vertices; v++)
    st.insert_simplex_and_subfaces(rand_int_range<Vertex_handle>(2, 5, nb_vertices));
  std::clog << "... " << st.num_vertices() << " vertices and " << st.num_simplices()
            << " simplices. Dimension is " << st.dimension() << "\n";

  for (int dim = 0; dim < 7; dim++) {
    std::clog << "... benchmark for dimension " << dim << "\n";
    int count = 0;
    Gudhi::Clock benchmark_skeleton_exact("Benchmark the skeleton range but with the exact dim test");
    for(auto sh : st.skeleton_simplex_range(dim)) {
      // Get only the simplices with the exact dimension
      if (st.dimension(sh) == dim)
        ++count;
    }
    std::clog << benchmark_skeleton_exact << " - count = " << count << "\n";

    count = 0;
    Gudhi::Clock benchmark_dimension("Benchmark the dimension range");
    for([[maybe_unused]] auto sh : st.dimension_simplex_range(dim)) {
      ++count;
    }
    std::clog << benchmark_dimension << " - count = " << count << "\n";
  }
  return EXIT_SUCCESS;
}
