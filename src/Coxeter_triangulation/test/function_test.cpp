#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <iostream>

int main() {
  std::size_t m = 3, d = 5;
  Eigen::VectorXd center(d); center << 2, 1.5, -0.5, 4.5, -1;
  double radius = 5;
  Gudhi::coxeter_triangulation::Function_Sm_in_Rd fun_sphere(radius,
							     m,
							     d,
							     center);
  std::cout << "Value of fun_sphere at the origin:\n"
	    << fun_sphere(Eigen::VectorXd::Zero(d)) << "\n";
  std::cout << "Value of fun_sphere at the center:\n"
	    << fun_sphere(center) << "\n";
  std::cout << "Ambient dimension = " << fun_sphere.amb_d()
	    << ", codimension = " << fun_sphere.cod_d() << "\n";
  std::cout << "Seed point:\n" << fun_sphere.seed() << "\nValue at seed point:\n"
	    << fun_sphere(fun_sphere.seed()) << "\n";
  
  return 0;
}
