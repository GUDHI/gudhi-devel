
// workaround for the annoying boost message in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>
// end workaround 

#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/Function_chair_in_R3.h>
#include <gudhi/Functions/Function_torus_in_R3.h>
#include <gudhi/Functions/Function_whitney_umbrella_in_R3.h>
#include <gudhi/Functions/Function_lemniscate_revolution_in_R3.h>
#include <gudhi/Functions/Function_iron_in_R3.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <iostream>
#include <string>

#include <random>
#include <cstdlib>

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
    std::cout << "Function sphere:\n";
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
    std::cout << "Function plane:\n";
    print_test(fun_plane);
  }
  {
    // the constant function testing part
    std::size_t k = 2, d = 5;
    auto x = Eigen::VectorXd::Constant(2, 1);
    Gudhi::coxeter_triangulation::Constant_function fun_const(d, k, x);
    std::cout << "Constant function:\n";
    print_test(fun_const);
  }
  {
    // the chair function
    Gudhi::coxeter_triangulation::Function_chair_in_R3 fun_chair;
    std::cout << "Function chair:\n";
    print_test(fun_chair);
  }
  {
    // the torus function
    Gudhi::coxeter_triangulation::Function_torus_in_R3 fun_torus;
    std::cout << "Function torus:\n";
    print_test(fun_torus);
  }
  {
    // the whitney umbrella function
    Gudhi::coxeter_triangulation::Function_whitney_umbrella_in_R3 fun_umbrella;
    std::cout << "Function Whitney umbrella:\n";
    print_test(fun_umbrella);
  }
  {
    // the lemniscate revolution function
    Gudhi::coxeter_triangulation::Function_lemniscate_revolution_in_R3 fun_lemniscate;
    std::cout << "Function Revolution surface of Lemniscate of Bernoulli:\n";
    print_test(fun_lemniscate);
  }
  {
    // the iron function
    Gudhi::coxeter_triangulation::Function_iron_in_R3 fun_iron;
    std::cout << "Function iron:\n";
    print_test(fun_iron);
  }
  // {
  //   // random orthogonal matrix
  //   auto matrix = random_orthogonal_matrix(5);
  //   std::cout << "M:\n" << matrix << "\nM^t * M\n" << matrix.transpose() * matrix << "\n";
  // }
  return 0;
}
