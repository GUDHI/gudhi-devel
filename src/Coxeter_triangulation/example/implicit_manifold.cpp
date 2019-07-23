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

/* An alternative version of the Cartesian product that takes two
 * pointers as arguments. */
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

/* An alternative version of Translate that takes a
 * pointer as argument. */
struct Translate_function : public Function {

  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    fun_->evaluate(p - off_, result);
  }
  
  /* Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_->amb_d();}

  /* Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_->cod_d();}

  void seed(Eigen::VectorXd& result) {
    fun_->seed(result);
    result += off_;
  }

  Translate_function* clone() const {
    return new Translate_function(*this);
  }

  Translate_function(Function* fun, const Eigen::VectorXd& off) :
    fun_(fun), off_(off) {
  }

private:
  Function *fun_;
  Eigen::VectorXd off_;
};

/* An alternative version of Embed_in_Rd that takes a
 * pointer as argument. */
struct Embed_function : public Function {

  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    Eigen::VectorXd x = p;
    Eigen::VectorXd x_k(fun_->amb_d()), x_rest(d_ - fun_->amb_d());
    for (std::size_t i = 0; i < fun_->amb_d(); ++i)
      x_k(i) = x(i);
    for (std::size_t i = fun_->amb_d(); i < d_; ++i)
      x_rest(i - fun_->amb_d()) = x(i);
    fun_->evaluate(x_k, result);
    result.conservativeResize(this->cod_d());
    for (std::size_t i = fun_->cod_d(); i < this->cod_d(); ++i)
      result(i) = x_rest(i - fun_->cod_d());
  }
  
  /* Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return d_;}

  /* Returns the codomain dimension. */
  std::size_t cod_d() const {return d_ - (fun_->amb_d() - fun_->cod_d());}

  void seed(Eigen::VectorXd& result) {
    fun_->seed(result);
    result.conservativeResize(d_);
    for (std::size_t l = fun_->amb_d(); l < d_; ++l)
      result(l) = 0;
  }

  Embed_function* clone() const {
    return new Embed_function(*this);
  }

  Embed_function(Function* fun, std::size_t d) :
    fun_(fun), d_(d) {
  }

private:
  Function *fun_;
  std::size_t d_;
};

/* An alternative version of Negation that takes a
 * pointer as argument. */
struct Negation_function : public Function {

  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    fun_->evaluate(p, result);
    result = -result;
  }
  
  /* Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_->amb_d();}

  /* Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_->cod_d();}

  void seed(Eigen::VectorXd& result) {
    fun_->seed(result);
  }

  Negation_function* clone() const {
    return new Negation_function(*this);
  }

  Negation_function(Function* fun) :
    fun_(fun) {
  }

private:
  Function *fun_;
};


bool read_symbol(std::ifstream& stream, char& symbol) {
  if (!(stream >> symbol)) {
    std::cerr << "Unexpected end of file.\n";
    return false;
  }
  std::cout << "Read " << symbol << "\n";
  return true;
}

void read_vector(std::ifstream& stream, Eigen::VectorXd& result) {
  std::vector<double> coords;
  char symbol;
  if (!read_symbol(stream, symbol))
    return;
  if (symbol != '(') {
    std::cerr << "RV: Invalid vector specification\n";
    return;
  }
  double x;
  while (stream >> x) {
    std::cout << "RV: Read " << x << "\n";
    coords.push_back(x);
    if (!read_symbol(stream, symbol))
      return;
    if (symbol != ')' && symbol != ',') {
      std::cerr << "RV: Unrecognized symbol " << symbol << ". Expected ')' or ','\n";
      return;
    }
    std::cout << "RV: Read " << symbol << "\n";
    if (symbol == ')') {
      result.resize(coords.size());
      std::size_t i = 0;
      for (double& x: coords)
	result(i++) = x;
      return;
    }
  }
  if (!read_symbol(stream, symbol))
    return;
  std::cout << "RV: Read " << symbol << "\n";
  if (symbol != ')') {
    std::cerr << "RV: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return;
  }
  result.resize(coords.size());
  std::size_t i = 0;
  for (double& x: coords)
    result(i++) = x;
}

void read_scalar(std::ifstream& stream, double& result) {
  char symbol;
  if (!read_symbol(stream, symbol))
    return;
  if (symbol != '(') {
    std::cerr << "RS: Unrecognized symbol " << symbol << ". Expected '('\n";
    return;
  }
  if (!(stream >> result))
    std::cerr << "RS: Invalid radius specification\n";
  std::cout << "RS: Read " << result << "\n";
  if (!read_symbol(stream, symbol))
    return;
  if (symbol != ')') {
    std::cerr << "RS: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return;
  }
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
      read_vector(stream, center);
      break;
    case 'r':
      read_scalar(stream, radius);
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

It parse_T(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
  double Radius, radius;
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "T: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "T: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, center);
      break;
    case 'r':
      read_scalar(stream, radius);
      break;
    case 'R':
      read_scalar(stream, Radius);
      break;  
    case ')':
      fun_range.emplace_front(new Function_torus_in_R3(Radius, radius, center));
      return fun_range.begin();
    default:
      std::cerr << "T: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "T: Unexpected end of file.\n";
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
    case 'T':
      fun = parse_T(stream, fun_range);
      break;
    case '(':
      fun = parse_in_parenthesis(stream, fun_range);
      break;      
    case ')':
      return fun;
    case '-': {
      fun = parse_in_parenthesis(stream, fun_range);
      fun_range.emplace_front(new Negation_function(*fun));
      fun = fun_range.begin();
      return fun;
    }
    case 'x': {
      It fun2 = parse_in_parenthesis(stream, fun_range);
      fun_range.emplace_front(new Cartesian_product_two(*fun, *fun2));
      fun = fun_range.begin();
      return fun;
    }
    case 't': {
      Eigen::VectorXd v;
      read_vector(stream, v);
      fun_range.emplace_front(new Translate_function(*fun, v));
      fun = fun_range.begin();
      break;
    }
    case 'R': {
      std::size_t d;
      if (!(stream >> d)) {
	std::cerr << "R: Expected dimension\n";
	return fun;
      }
      std::cout << "R: Read " << d << "\n";
      fun_range.emplace_front(new Embed_function(*fun, d));
      fun = fun_range.begin();
      break;
    }
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
    case 'T':
      fun = parse_T(stream, fun_range);
      break;
    case '(':
      fun = parse_in_parenthesis(stream, fun_range);
      break;
    case '-': {
      fun = parse_in_parenthesis(stream, fun_range);
      fun_range.emplace_front(new Negation_function(*fun));
      fun = fun_range.begin();
      return fun;
    }
    case 'x': {
      It fun2 = parse_input(stream, fun_range);
      fun_range.emplace_front(new Cartesian_product_two(*fun, *fun2)); 
      fun = fun_range.begin();
      break;
    }
    case 't': {
      Eigen::VectorXd v;
      read_vector(stream, v);
      fun_range.emplace_front(new Translate_function(*fun, v));
      fun = fun_range.begin();
      break;
    }
    case 'R': {
      std::size_t d;
      if (!(stream >> d)) {
	std::cerr << "R: Expected dimension\n";
	return fun;
      }
      std::cout << "R: Read " << d << "\n";
      fun_range.emplace_front(new Embed_function(*fun, d));
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
