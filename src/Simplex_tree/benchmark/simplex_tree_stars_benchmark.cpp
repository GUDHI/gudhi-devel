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



template <class Stree>
void benchmark_stars_computation() {
  using Vertex_handle  = typename Stree::Vertex_handle;
  using Simplex_handle = typename Stree::Simplex_handle;

  std::clog << "... Simplex_tree construction" << std::endl;

  Stree st;
  // Insert 5.000 random simplices, of size in between [2; 5] and vertices in between [0; 5.000]
  for (Vertex_handle v=0; v < 5000; v++)
    st.insert_simplex_and_subfaces(rand_int_range<Vertex_handle>(2, 5, 5000));
  std::cout << "... " << st.num_vertices() << " vertices and " << st.num_simplices() << " simplices." << std::endl;
  
  Gudhi::Clock benchmark_search("... Looking for random existing simplices");
  // 5.000 random simplices of size in between [2; 5] has been inserted - about 50.000 simplices in the Simplex_tree
  const int SH_SIZE = 20000;
  std::vector<Simplex_handle> sh_list;
  sh_list.reserve(SH_SIZE);

  std::sample(st.complex_simplex_range().begin(), st.complex_simplex_range().end(), std::back_inserter(sh_list),
              SH_SIZE, std::mt19937 {rd()});
  std::clog << benchmark_search << std::endl;

  Gudhi::Clock benchmark_stars("Benchmark the stars search of the random simplices");
  // Just browse the stars from the random simplices list
  for (auto& sh : sh_list)
    for ([[maybe_unused]] const auto& simplex : st.star_simplex_range(sh))
      [[maybe_unused]] volatile auto _ = simplex;

  std::clog << benchmark_stars << std::endl;
}

int main() {
  std::clog << "** Without stars computation optimization" << std::endl;
  benchmark_stars_computation<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_full_featured>>();

  std::clog << "** With stars computation optimization" << std::endl;
  benchmark_stars_computation<Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_cofaces>>();

  return 0;
}
