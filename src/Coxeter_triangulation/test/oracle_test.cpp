#include <iostream>
#include <string>

#include <gudhi/Implicit_manifold_intersection_oracle.h>

#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Cartesian_product.h>

#include <gudhi/Coxeter_triangulation.h>

#include <random>
#include <cstdlib>

using namespace Gudhi::coxeter_triangulation;

int main() {
  Function_Sm_in_Rd fun_sph(5.1111, 2);
  auto oracle = make_oracle(fun_sph, 0.3);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  
  std::cout << std::boolalpha;
  auto s = cox_tr.locate_point(fun_sph.seed());
  std::cout << "Full simplex = \033[1;31m" << s << "\033[0m:\n\n";
  
  for (auto f: s.face_range(oracle.cod_d())) {
    std::cout << "Query simplex = \033[1;34m" << f << "\033[0m:\n";
    auto qr = oracle.intersects(f, cox_tr);
    std::cout << "Face: " << qr.face << "\n";
    std::cout << "Intersection point:\n" << qr.intersection << "\n";
    std::cout << "Success: " << qr.success << "\n\n";
  }
  
  return 0;
}
