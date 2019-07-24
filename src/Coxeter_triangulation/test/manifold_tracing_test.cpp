#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>

using namespace Gudhi::coxeter_triangulation;

int main() {
  Function_Sm_in_Rd fun_sph(5.1111, 2);
  auto oracle = make_oracle(fun_sph, 0.3);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  Eigen::VectorXd seed;
  fun_sph.seed(seed);
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  typename MT::Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, out_simplex_map);

  for (auto si_pair: out_simplex_map)
    std::cout << "Simplex = \033[1;34m" << si_pair.first << "\033[0m:\n"
	      << "Intersection point:\n" << si_pair.second << "\n";
  
  std::cout << "\nSize of the output = " << out_simplex_map.size() << "\n";
  
  return 0;
}
