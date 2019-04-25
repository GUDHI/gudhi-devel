#include <iostream>

#include <gudhi/Freudenthal_triangulation.h>
#include <gudhi/Coxeter_triangulation.h>
#include "../example/cxx-prettyprint/prettyprint.hpp"

int main() {
  // typedef Gudhi::Freudenthal_triangulation Triangulation;
  // typedef typename Triangulation::Simplex_handle Simplex;
  // Triangulation triangulation(8);
  // Simplex s = triangulation.locate_point(std::vector<double>({5.5, 0.75, 4.5, 2.75, 3.15, -8, 6.25, 3.25}));
  // std::cout << s.vertex << " ";
  // for (auto p: s.partition)
  //   std::cout << p << " ";
  // std::cout << "\n";
  // std::cout << "Vertices:\n";
  // for (auto v: triangulation.vertex_range(s))
  //   std::cout << v << "\n" << triangulation.cartesian_coordinates(v) << "\n";
  // std::cout << triangulation.barycenter(s) << "\n";
  std::vector<double> x = {5.5, 0.75, 4.5};
  // std::vector<double> x = {5.5, 0.75, 4.5, 2.75, 3.15, -8, 6.25, 3.25};
  std::cout << Gudhi::Coxeter_triangulation(3).locate_point(x).vertex << "\n"; 
  
  // std::vector<int> v = {-6, 5, -1, 0};
  // std::cout << "vertex = " << v << "\n";
  // Simplex s = triangulation.locate_point(triangulation.cartesian_coordinates(v, 3), 3);
  // std::cout << s.vertex << " ";
  // for (auto p: s.partition)
  //   std::cout << p << " ";
  // std::cout << "\n";
  
}
