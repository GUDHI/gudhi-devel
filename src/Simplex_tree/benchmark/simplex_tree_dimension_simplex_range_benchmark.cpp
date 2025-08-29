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

#include "random_simplices.h"

#include <iostream>
#include <vector>
#include <cstdlib>

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
    st.insert_simplex_and_subfaces(random_simplex<Vertex_handle>(2, 5, nb_vertices));
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
