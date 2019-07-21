#include <string>
#include <iostream>
#include <fstream>
#include <tuple>

#include <boost/filesystem.hpp>
#include <gudhi/Functions/Function.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/Function_chair_in_R3.h>
#include <gudhi/Functions/Function_torus_in_R3.h>
#include <gudhi/Functions/Function_whitney_umbrella_in_R3.h>
#include <gudhi/Functions/Function_lemniscate_revolution_in_R3.h>
#include <gudhi/Functions/Function_iron_in_R3.h>
#include <gudhi/Functions/Function_moment_curve_in_Rd.h>
#include <gudhi/Functions/Embed_in_Rd.h>
#include <gudhi/Functions/Translate.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Functions/Negation.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Functions/PL_approximation.h>

#include <gudhi/Cell_complex.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** const argv) {
  if (argc < 2) {
    std::cout << "Too few arguments. Usage: " << argv[0] << " input_file_name";
    return 1;
  }
  boost::filesystem::path path(argv[1]);
  std::string input_file_name = path.string();

  std::vector<Function*> functions;
  functions.emplace_back(new Function_Sm_in_Rd(5.0, 3, 5));
  functions.emplace_back(new Function_iron_in_R3());
  std::cout << "The result for amb_d = " << functions.back()->amb_d() << "\n";
  std::cout << "The result for cod_d = " << functions.back()->cod_d() << "\n";
  std::cout << "The result for f(0) =\n" << (*functions.back())(Eigen::VectorXd::Zero(5)) << "\n";
  
  std::ifstream stream(input_file_name);
  if (!stream.is_open()) {
    std::cerr << argv[0] << ": can't open file " << input_file_name << "\n";
    return 2;
  }

  for (auto fun: functions)
    delete fun;
  
  stream.close();
  return 0;
}
