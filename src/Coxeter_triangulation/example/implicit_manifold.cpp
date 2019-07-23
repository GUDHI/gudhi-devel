#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <list>

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

using Function_range = std::list<Function*>;
using It = typename Function_range::iterator;

struct Cartesian_product_two : public Function {

  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    result.resize(cod_d_);
    std::size_t d1 = fun1_->amb_d(), d2 = fun2_->amb_d();
    std::size_t k1 = fun1_->cod_d(), k2 = fun2_->cod_d();
    Eigen::VectorXd p1(d1), p2(d2), res1, res2;
    for (std::size_t i = 0; i < d1; ++i)
      p1(i) = p(i);
    for (std::size_t i = 0; i < d2; ++i)
      p2(i) = p(i + d1);
    fun1_->evaluate(p1, res1);
    fun2_->evaluate(p2, res2);
    for (std::size_t i = 0; i < k1; ++i)
      result(i) = res1(i);
    for (std::size_t i = 0; i < k2; ++i)
      result(i+k1) = res2(i);
  }
  
  /* Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return amb_d_;}

  /* Returns the codomain dimension. */
  std::size_t cod_d() const {return cod_d_;}

  void seed(Eigen::VectorXd& result) {
    result.resize(amb_d_);
    std::size_t d1 = fun1_->amb_d(), d2 = fun2_->amb_d();
    Eigen::VectorXd res1(d1), res2(d2);
    fun1_->seed(res1);
    fun2_->seed(res2);
    for (std::size_t i = 0; i < d1; ++i)
      result(i) = res1(i);
    for (std::size_t i = 0; i < d2; ++i)
      result(i + d1) = res2(i);
  }

  Cartesian_product_two* clone() const {
    return new Cartesian_product_two(*this);
  }

  
  Cartesian_product_two(Function* fun1, Function* fun2)
    : fun1_(fun1),
      fun2_(fun2),
      amb_d_(fun1_->amb_d() + fun2_->amb_d()),
      cod_d_(fun1_->cod_d() + fun2_->cod_d())
  {}

private:
  Function *fun1_, *fun2_;
  std::size_t amb_d_, cod_d_;
};

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

It parse_S(std::ifstream& stream, Function_range& fun_range) {
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
      fun_range.emplace_front(new Function_Sm_in_Rd(radius, m, center));
      return fun_range.begin();
    default:
      std::cerr << "S: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "S: Unexpected end of file.\n";
  return fun;
}

It parse_in_parenthesis(std::ifstream& stream, Function_range& fun_range) {
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

It parse_input(std::ifstream& stream, Function_range& fun_range) {
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
    case 'x': {
      It fun2 = parse_input(stream, fun_range);
      fun_range.emplace_front(new Cartesian_product_two(*fun, *fun2));
      fun = fun_range.begin();
      break;
    }
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
  Eigen::VectorXd value;
  (*fun)->evaluate(Eigen::VectorXd::Zero(amb_d), value);
  
  std::cout << "amb_d = " << amb_d << "\n";
  std::cout << "cod_d = " << (*fun)->cod_d() << "\n";
  std::cout << "f(0) =\n" << value << "\n";
   
  for (auto fun: fun_range)
    delete fun;
  
  stream.close();
  return 0;
}
