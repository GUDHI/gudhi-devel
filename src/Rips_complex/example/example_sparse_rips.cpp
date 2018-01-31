#include <gudhi/Sparse_rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>

#include <iostream>
#include <vector>

int main() {
  using Point = std::vector<double>;
  using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
  using Filtration_value = Simplex_tree::Filtration_value;
  using Complex = Gudhi::rips_complex::Sparse_rips_complex<Filtration_value>;

  Point points[] = {
    {1.0, 1.0},
    {7.0, 0.0},
    {4.0, 6.0},
    {9.0, 6.0},
    {0.0, 14.0},
    {2.0, 19.0},
    {9.0, 17.0}};

  // ----------------------------------------------------------------------------
  // Init from Euclidean points
  // ----------------------------------------------------------------------------
  double epsilon = 2; // very rough, no guarantees
  Complex cpx(points, Gudhi::Euclidean_distance(), epsilon);

  Simplex_tree stree;
  cpx.create_complex(stree, 10);

  // ----------------------------------------------------------------------------
  // Display information about the complex
  // ----------------------------------------------------------------------------
  std::cout << "Sparse Rips complex is of dimension " << stree.dimension() <<
               " - " << stree.num_simplices() << " simplices - " <<
               stree.num_vertices() << " vertices." << std::endl;
}
