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

using Function_range = std::vector<Function*>;
using It = typename Function_range::iterator;

bool read_symbol(std::ifstream& stream, char& symbol) {
  if (!(stream >> symbol)) {
    std::cerr << "Unexpected end of file.\n";
    return false;
  }
  std::cout << "Read " << symbol << "\n";
  return true;
}

double read_sphere_radius(std::ifstream& stream) {
  double radius;
  char symbol;
  if (!read_symbol(stream, symbol))
    return 0;
  if (symbol != '(') {
    std::cerr << "SR: Unrecognized symbol " << symbol << ". Expected '('\n";
    return 0;
  }
  if (!(stream >> radius))
    std::cerr << "SR: Invalid radius specification\n";
  std::cout << "SR: Read " << radius << "\n";
  if (!read_symbol(stream, symbol))
    return 0;
  if (symbol != ')') {
    std::cerr << "SR: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return 0;
  }
  return radius;
}

Eigen::VectorXd read_sphere_center(std::ifstream& stream) { 
  std::vector<double> coords;
  char symbol;
  if (!read_symbol(stream, symbol))
    return Eigen::VectorXd();
  if (symbol != '(') {
    std::cerr << "SC: Invalid center specification\n";
    return Eigen::VectorXd();
  }
  double x;
  while (stream >> x) {
    std::cout << "SC: Read " << x << "\n";
    coords.push_back(x);
    if (!read_symbol(stream, symbol))
      return Eigen::VectorXd();
    if (symbol != ')' && symbol != ',') {
      std::cerr << "SC: Unrecognized symbol " << symbol << ". Expected ')' or ','\n";
      return Eigen::VectorXd();
    }
    std::cout << "SC: Read " << symbol << "\n";
    if (symbol == ')') {
      Eigen::VectorXd v(coords.size());
      std::size_t i = 0;
      for (double& x: coords)
	v(i++) = x;
      return v;
    }
  }
  if (!read_symbol(stream, symbol))
    return Eigen::VectorXd();
  std::cout << "SC: Read " << symbol << "\n";
  if (symbol != ')') {
    std::cerr << "SC: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return Eigen::VectorXd();
  }
  Eigen::VectorXd v(coords.size());
  std::size_t i = 0;
  for (double& x: coords)
    v(i++) = x;
  return v;  
}

It parse_S(std::ifstream& stream, std::vector<Function*>& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center;
  double radius;
  std::size_t m;
  if (!(stream >> m)) {
    std::cerr << "S: Expected dimension\n";
    return fun;
  }
  std::cout << "S: Read " << m << "\n";
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "S: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "S: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      center = read_sphere_center(stream);
      break;
    case 'r':
      radius = read_sphere_radius(stream);
      break;
    case ')':
      fun_range.emplace_back(new Function_Sm_in_Rd(radius, m, center));
      return fun_range.end()-1;
    default:
      std::cerr << "S: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "S: Unexpected end of file.\n";
  return fun;
}

It parse_in_parenthesis(std::ifstream& stream, std::vector<Function*>& fun_range) {
  // fun_range.emplace_back(new Function());
  It fun = fun_range.begin();
  char symbol;
  while (stream >> symbol) {    
    std::cout << "P: Read " << symbol << "\n";
    switch (symbol) {
    case 'S':
      fun = parse_S(stream, fun_range);
      break;
    case ')':
      return fun;
    default:
      std::cerr << "P: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  return fun;
}

It parse_input(std::ifstream& stream, std::vector<Function*>& fun_range) {
  fun_range.emplace_back(new Constant_function());
  It fun = fun_range.begin();
  char symbol;
  while (stream >> symbol) {
    std::cout << "Read " << symbol << "\n";
    switch (symbol) {
    case 'S':
      fun = parse_S(stream, fun_range);
      break;
    case '(':
      fun = parse_in_parenthesis(stream, fun_range);
      break;
    default:
      std::cerr << "Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  return fun;
}

int main(int argc, char** const argv) {
  if (argc < 2) {
    std::cout << "Too few arguments. Usage: " << argv[0] << " input_file_name";
    return 1;
  }
  std::cout << "sizeof(Function*)" << sizeof(Function*) << "\n";
  std::cout << "sizeof(Function)" << sizeof(Function) << "\n";
  std::cout << "sizeof(Function_Sm_in_Rd*)" << sizeof(Function_Sm_in_Rd*) << "\n";
  std::cout << "sizeof(Function_Sm_in_Rd)" << sizeof(Function_Sm_in_Rd) << "\n";
  boost::filesystem::path path(argv[1]);
  std::string input_file_name = path.string();

  Function_range fun_range;
  
  std::ifstream stream(input_file_name);
  if (!stream.is_open()) {
    std::cerr << argv[0] << ": can't open file " << input_file_name << "\n";
    return 2;
  }
  It fun = parse_input(stream, fun_range);
  std::size_t amb_d = (*fun)->amb_d();
  // Eigen::VectorXd value = (**fun)(Eigen::VectorXd::Zero(amb_d));
  // std::cout << "amb_d = " << amb_d << "\n";
  // std::cout << "cod_d = " << (*fun)->cod_d() << "\n";
  // std::cout << "f(0) =\n" << value << "\n";
  std::cout << "f(0) =\n" << Function_Sm_in_Rd(1.0, 2)(Eigen::VectorXd::Zero(3)) << "\n";
  
  // fun_range.emplace_back(new Function_Sm_in_Rd(2.3, 2));
  // std::cout << "amb_d = " << (*fun_range.begin())->amb_d() << "\n";
  // Function* fun = (*fun_range.begin())->clone();
  // delete fun;
  // std::cout << "f(0) =\n" << (**fun)(Eigen::VectorXd::Zero((*fun)->amb_d())) << "\n";
  
  for (auto fun: fun_range)
    delete fun;
  
  stream.close();
  return 0;
}
