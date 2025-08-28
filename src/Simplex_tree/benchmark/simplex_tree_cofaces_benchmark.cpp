/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - 2025/08 Vincent Rouvreau: externalize rand_int_range in random_simplices.h for DRY purposes
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Clock.h>

#include "random_simplices.h"

#include <iostream>
#include <algorithm>  // for std::sample
#include <vector>
#include <cstdlib>

template <class Stree>
void benchmark_stars_computation(int nb_vertices) {
  using Vertex_handle  = typename Stree::Vertex_handle;
  using Simplex_handle = typename Stree::Simplex_handle;

  std::clog << "... Simplex_tree construction" << std::endl;

  Stree st;
  // Insert 'nb_vertices' random simplices, of size in between [2; 5] and vertices in between [0; nb_vertices]
  for (Vertex_handle v=0; v < nb_vertices; v++)
    st.insert_simplex_and_subfaces(random_simplex<Vertex_handle>(2, 5, nb_vertices));
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
    for (const auto& simplex : st.star_simplex_range(sh))
      [[maybe_unused]] volatile auto _ = simplex;

  std::clog << benchmark_stars << std::endl;

  Gudhi::Clock benchmark_cofaces("Benchmark the cofaces (codimension=1) search of the random simplices");
  // Just browse the stars from the random simplices list
  for (auto& sh : sh_list)
    for (const auto& simplex : st.cofaces_simplex_range(sh, 1))
      [[maybe_unused]] volatile auto _ = simplex;

  std::clog << benchmark_cofaces << std::endl;
}

struct Stree_basic_cofaces_options {
  typedef Gudhi::linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef double Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = true;
  static const bool contiguous_vertices = false;
  static const bool link_nodes_by_label = false;
  static const bool stable_simplex_handles = false;
};

struct Stree_fast_cofaces_options : Stree_basic_cofaces_options {
  static const bool link_nodes_by_label = true;
};

struct Stree_fast_cofaces_stable_simplex_handles_options : Stree_basic_cofaces_options {
  static const bool link_nodes_by_label = true;
  static const bool stable_simplex_handles = true;
};

int main(int argc, char *argv[]) {
  int nb_vertices = 5000;
  if (argc > 2) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct\n";
    std::cerr << "Usage: " << (argv[0] - 1) << " [NB_VERTICES] \n";
    std::cerr << "    NB_VERTICES is 5.000 by default.\n";
    exit(EXIT_FAILURE);  // ----- >>
  }
  if (argc == 2)
    nb_vertices = atoi(argv[1]);

  std::clog << "** Without cofaces computation optimization" << std::endl;
  benchmark_stars_computation<Gudhi::Simplex_tree<Stree_basic_cofaces_options>>(nb_vertices);

  std::clog << "** With cofaces computation optimization" << std::endl;
  benchmark_stars_computation<Gudhi::Simplex_tree<Stree_fast_cofaces_options>>(nb_vertices);

  std::clog << "** With cofaces computation optimization and stable simplex handles" << std::endl;
  benchmark_stars_computation<Gudhi::Simplex_tree<Stree_fast_cofaces_stable_simplex_handles_options> >(nb_vertices);

  return EXIT_SUCCESS;
}
