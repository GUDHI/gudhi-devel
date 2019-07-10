#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Constant_function.h>
#include <iostream>
#include <string>

template <class Function>
void print_test(const Function& fun) {
  std::cout << "Ambient dimension = " << fun.amb_d()
	    << ", codimension = " << fun.cod_d() << "\n";
  std::cout << "Value of fun at the origin:\n"
	    << fun(Eigen::VectorXd::Zero(fun.amb_d())) << "\n";
  try {
    std::cout << "Seed point:\n" << fun.seed() << "\nValue at seed point:\n"
	      << fun(fun.seed()) << "\n";
  }
  catch (...) {
    std::cout << "Caught an exception!\n";
  }
  std::cout << "\n";
}

int main() {
  {
    // the sphere testing part
    std::size_t m = 3, d = 5;
    Eigen::VectorXd center(d); center << 2, 1.5, -0.5, 4.5, -1;
    double radius = 5;
    typedef Gudhi::coxeter_triangulation::Function_Sm_in_Rd Function_sphere;
    Function_sphere fun_sphere(radius, m, d, center);
    print_test(fun_sphere);
  }
  {
    // the affine plane testing part
    std::size_t m = 0, d = 5;
    Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(d, d-m);
    for (std::size_t i = 0; i < d-m; ++i)
      normal_matrix(i,i) = 1;
    typedef Gudhi::coxeter_triangulation::Function_affine_plane_in_Rd Function_plane;
    Function_plane fun_plane(normal_matrix);
    print_test(fun_plane);
  }
  {
    // the constant function testing part
    std::size_t k = 2, d = 5;
    auto x = Eigen::VectorXd::Constant(2, 1);
    Gudhi::coxeter_triangulation::Constant_function fun_const(d, k, x);
    print_test(fun_const);
  }
  return 0;
}
