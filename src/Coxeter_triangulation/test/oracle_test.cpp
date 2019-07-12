#include <gudhi/Implicit_manifold_intersection_oracle.h>

#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Functions/Domain_from_function.h>

#include <iostream>
#include <string>

#include <random>
#include <cstdlib>

using namespace Gudhi::coxeter_triangulation;

int main() {
  Function_Sm_in_Rd fun_sph(1.0, 1);
  auto oracle = make_oracle(fun_sph);
  return 0;
}
