/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

#include <iostream>
#include <vector>
#include <cstdint>  // for std::uint8_t

/* We could perfectly well use the default Simplex_tree<> (which uses
 * Simplex_tree_options_full_featured), the following simply demonstrates
 * how to save on storage by not storing a filtration value.  */

struct MyOptions : Gudhi::Simplex_tree_options_full_featured {
  // Implicitly use 0 as filtration value for all simplices
  static const bool store_filtration = false;
  // The persistence algorithm needs this
  static const bool store_key = true;
  // I have few vertices
  typedef short Vertex_handle;
  // Maximum number of simplices to compute persistence is 2^8 - 1 = 255. One is reserved for null_key
  typedef std::uint8_t Simplex_key;
};

using ST = Gudhi::Simplex_tree<MyOptions>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<ST, Field_Zp>;

int main() {
  ST st;

  /* Complex to build. */
  /*    1   3   5      */
  /*    o---o---o      */
  /*   / \ /           */
  /*  o---o   o        */
  /*  2   0   4        */

  const short edge01[] = {0, 1};
  const short edge02[] = {0, 2};
  const short edge12[] = {1, 2};
  const short edge03[] = {0, 3};
  const short edge13[] = {1, 3};
  const short edge35[] = {3, 5};
  const short vertex4[] = {4};
  st.insert_simplex_and_subfaces(edge01);
  st.insert_simplex_and_subfaces(edge02);
  st.insert_simplex_and_subfaces(edge12);
  st.insert_simplex_and_subfaces(edge03);
  st.insert_simplex(edge13);
  st.insert_simplex_and_subfaces(edge35);
  st.insert_simplex(vertex4);

  // Class for homology computation
  // By default, since the complex has dimension 1, only 0-dimensional homology would be computed.
  // Here we also want persistent homology to be computed for the maximal dimension in the complex (persistence_dim_max = true)
  Persistent_cohomology pcoh(st, true);

  // Initialize the coefficient field Z/2Z for homology
  pcoh.init_coefficients(2);

  // Compute the persistence diagram of the complex
  pcoh.compute_persistent_cohomology();

  // Print the result. The format is, on each line: 2 dim 0 inf
  // where 2 represents the field, dim the dimension of the feature.
  // 2  0 0 inf
  // 2  0 0 inf
  // 2  1 0 inf
  // 2  1 0 inf
  // means that in Z/2Z-homology, the Betti numbers are b0=2 and b1=2.
  pcoh.output_diagram();

  // Print the Betti numbers are b0=2 and b1=2.
  std::clog << std::endl;
  std::clog << "The Betti numbers are : ";
  for (int i = 0; i < 3; i++)
    std::clog << "b" << i << " = " << pcoh.betti_number(i) << " ; ";
  std::clog << std::endl;
}
