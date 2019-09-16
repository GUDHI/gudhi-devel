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

#include <gudhi/Query_result.h>
#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h> 

using namespace Gudhi::coxeter_triangulation;

using Function_range = std::list<Function*>;
using It = typename Function_range::iterator;

/******************************* DYNAMIC VERSIONS OF FUNCTION MODIFIERS ****************/

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

  void seed(Eigen::VectorXd& result) const {
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

  void seed(Eigen::VectorXd& result) const {
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

  void seed(Eigen::VectorXd& result) const {
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

  void seed(Eigen::VectorXd& result) const {
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

/* An alternative version of Linear_transformation that takes a
 * pointer as argument. */
struct Linear_transformation_function : public Function {

  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    fun_->evaluate(matrix_.householderQr().solve(p), result);
  }
  
  /* Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_->amb_d();}

  /* Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_->cod_d();}

  void seed(Eigen::VectorXd& result) const {
    fun_->seed(result);
    result = matrix_ * result;
  }

  Linear_transformation_function* clone() const {
    return new Linear_transformation_function(*this);
  }

  Linear_transformation_function(Function* fun, const Eigen::MatrixXd& matrix) :
    fun_(fun), matrix_(matrix) {
  }

private:
  Function *fun_;
  Eigen::MatrixXd  matrix_;
};

/******************************* PARSER *****************************************/

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
}

void read_matrix(std::ifstream& stream, Eigen::MatrixXd& result) {
  std::size_t cols = result.cols();
  std::vector<double> coords;
  char symbol;
  if (!read_symbol(stream, symbol))
    return;
  if (symbol != '(') {
    std::cerr << "RM: Invalid vector specification\n";
    return;
  }
  double x;
  while (stream >> x) {
    std::cout << "RM: Read " << x << "\n";
    coords.push_back(x);
    if (!read_symbol(stream, symbol))
      return;
    if (symbol != ')' && symbol != ',') {
      std::cerr << "RM: Unrecognized symbol " << symbol << ". Expected ')' or ','\n";
      return;
    }
    std::cout << "RM: Read " << symbol << "\n";
    if (symbol == ')') {
      std::size_t i = 0, j = 0;
      for (double& x: coords) {
	result(j, i++) = x;
	if (i == cols) {
	  j++;
	  i = 0;
	}
      }
      return;
    }
  }
  if (!read_symbol(stream, symbol))
    return;
  std::cout << "RM: Read " << symbol << "\n";
  if (symbol != ')') {
    std::cerr << "RM: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return;
  }
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

It parse_A(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  std::size_t m, d;
  if (!(stream >> m)) {
    std::cerr << "A: Expected intrinsic dimension\n";
    return fun;
  }
  std::cout << "A: Read " << m << "\n";
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != 'R') {
    std::cerr << "A: Unrecognized symbol " << symbol << ". Expected 'R'\n";
    return fun;
  }
  if (!(stream >> d)) {
    std::cerr << "A: Expected ambient dimension\n";
    return fun;
  }
  std::cout << "A: Read " << d << "\n";

  Eigen::VectorXd offset = Eigen::VectorXd::Zero(d);
  Eigen::MatrixXd matrix(d, d-m);
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "A: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "A: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, offset);
      break;
    case 'N':
      read_matrix(stream, matrix);
      break;
    case ')':
      fun_range.emplace_front(new Function_affine_plane_in_Rd(matrix, offset));
      return fun_range.begin();
    default:
      std::cerr << "A: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "A: Unexpected end of file.\n";
  return fun;
}

It parse_C(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
  double a = 0.8, b = 0.4, k = 1.0;
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "C: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "C: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, center);
      break;
    case 'a':
      read_scalar(stream, a);
      break;
    case 'b':
      read_scalar(stream, b);
      break;  
    case 'k':
      read_scalar(stream, k);
      break;  
    case ')':
      fun_range.emplace_front(new Function_chair_in_R3(a, b, k, center));
      return fun_range.begin();
    default:
      std::cerr << "C: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "C: Unexpected end of file.\n";
  return fun;
}

It parse_I(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "I: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "I: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, center);
      break;
    case ')':
      fun_range.emplace_front(new Function_iron_in_R3(center));
      return fun_range.begin();
    default:
      std::cerr << "I: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "I: Unexpected end of file.\n";
  return fun;
}

It parse_L(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
  double a = 1.0;
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "L: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "L: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, center);
      break;
    case 'a':
      read_scalar(stream, a);
      break;
    case ')':
      fun_range.emplace_front(new Function_lemniscate_revolution_in_R3(a, center));
      return fun_range.begin();
    default:
      std::cerr << "L: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "L: Unexpected end of file.\n";
  return fun;
}

It parse_M(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  // Eigen::VectorXd center = Eigen::Vector3d::Zero();
  double radius;
  std::size_t d;
  if (!(stream >> d)) {
    std::cerr << "M: Expected dimension\n";
    return fun;
  }
  std::cout << "M: Read " << d << "\n";
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "M: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "M: Read " << symbol << "\n";
    switch (symbol) {
    case 'r':
      read_scalar(stream, radius);
      break;
    case ')':
      fun_range.emplace_front(new Function_moment_curve_in_Rd(radius, d));
      return fun_range.begin();
    default:
      std::cerr << "M: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "M: Unexpected end of file.\n";
  return fun;
}

It parse_S(std::ifstream& stream, Function_range& fun_range) {
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
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
  double Radius = 1, radius = 0.5;
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
      fun_range.emplace_front(new Function_torus_in_R3(radius, Radius, center));
      return fun_range.begin();
    default:
      std::cerr << "T: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "T: Unexpected end of file.\n";
  return fun;
}

It parse_W(std::ifstream& stream, Function_range& fun_range) { 
  It fun = fun_range.begin();
  char symbol;
  Eigen::VectorXd center = Eigen::Vector3d::Zero();
  if (!read_symbol(stream, symbol))
    return fun;
  if (symbol != '(') {
    std::cerr << "W: Unrecognized symbol " << symbol << ". Expected ')'\n";
    return fun;
  }
  while (stream >> symbol) {
    std::cout << "W: Read " << symbol << "\n";
    switch (symbol) {
    case 'c':
      read_vector(stream, center);
      break;
    case ')':
      fun_range.emplace_front(new Function_whitney_umbrella_in_R3(center));
      return fun_range.begin();
    default:
      std::cerr << "W: Unrecognized symbol " << symbol << "\n";
      return fun;
    }
  }
  std::cerr << "W: Unexpected end of file.\n";
  return fun;
}

It parse_in_parenthesis(std::ifstream& stream, Function_range& fun_range) {
  // fun_range.emplace_back(new Function());
  It fun = fun_range.begin();
  char symbol;
  while (stream >> symbol) {    
    std::cout << "P: Read " << symbol << "\n";
    switch (symbol) {
    case 'A':
      fun = parse_A(stream, fun_range);
      break;
    case 'C':
      fun = parse_C(stream, fun_range);
      break;
    case 'I':
      fun = parse_I(stream, fun_range);
      break;
    case 'L':
      fun = parse_L(stream, fun_range);
      break;
    case 'M':
      fun = parse_M(stream, fun_range);
      break;
    case 'S':
      fun = parse_S(stream, fun_range);
      break;
    case 'T':
      fun = parse_T(stream, fun_range);
      break;
    case 'W':
      fun = parse_W(stream, fun_range);
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
    case 'm': {
      Eigen::MatrixXd matrix((*fun)->amb_d(), (*fun)->amb_d());
      read_matrix(stream, matrix);
      fun_range.emplace_front(new Linear_transformation_function(*fun, matrix));
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
    case 'A':
      fun = parse_A(stream, fun_range);
      break;
    case 'C':
      fun = parse_C(stream, fun_range);
      break;
    case 'I':
      fun = parse_I(stream, fun_range);
      break;
    case 'L':
      fun = parse_L(stream, fun_range);
      break;
    case 'M':
      fun = parse_M(stream, fun_range);
      break;
    case 'S':
      fun = parse_S(stream, fun_range);
      break;
    case 'T':
      fun = parse_T(stream, fun_range);
      break;
    case 'W':
      fun = parse_W(stream, fun_range);
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
    case 'm': {
      Eigen::MatrixXd matrix((*fun)->amb_d(), (*fun)->amb_d());
      read_matrix(stream, matrix);
      fun_range.emplace_front(new Linear_transformation_function(*fun, matrix));
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

/***************** DYNAMIC VERSION OF THE IMPLICIT MANIFOLD ORACLE ****************/

class Implicit_manifold_intersection_oracle {

  /* Computes the affine coordinates of the intersection point of the implicit manifold
   * and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_lambda(const Simplex_handle& simplex,
				 const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
    for (std::size_t i = 0; i < cod_d + 1; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords;
      fun_->evaluate(triangulation.cartesian_coordinates(v), v_coords);
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      j++;
    }
    Eigen::VectorXd z(cod_d + 1);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 1; ++i)
      z(i) = 0;
    return matrix.colPivHouseholderQr().solve(z);
  }

  /* Computes the affine coordinates of the intersection point of the boundary
   * of the implicit manifold and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_boundary_lambda(const Simplex_handle& simplex,
					  const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 2, cod_d + 2);
    for (std::size_t i = 0; i < cod_d + 2; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords;
      fun_->evaluate(triangulation.cartesian_coordinates(v), v_coords);
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      Eigen::VectorXd bv_coords;
      domain_fun_->evaluate(triangulation.cartesian_coordinates(v), bv_coords);
      matrix(cod_d + 1, j) = bv_coords(0);
      j++;
    }
    Eigen::VectorXd z(cod_d + 2);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 2; ++i)
      z(i) = 0;
    return matrix.colPivHouseholderQr().solve(z);
  }

  /* Computes the intersection result for a given simplex in a triangulation. */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersection_result(const Eigen::VectorXd& lambda,
						   const Simplex_handle& simplex,
						   const Triangulation& triangulation) const {
    using QR = Query_result<Simplex_handle>;
    std::size_t amb_d = triangulation.dimension();

    std::vector<std::size_t> snapping_indices;
    for (std::size_t i = 0; i < (std::size_t)lambda.size(); ++i) {
      if (lambda(i) < -threshold_ || lambda(i) > 1 + threshold_)
	return QR({Simplex_handle(), Eigen::VectorXd(), false});
      if (lambda(i) >= threshold_)
	snapping_indices.push_back(i);
    }

    std::size_t snap_d = snapping_indices.size();
    std::size_t i = 0;
    std::size_t num_line = 0;
    Eigen::MatrixXd vertex_matrix(snap_d, amb_d);
    Eigen::VectorXd reduced_lambda(snap_d);
    auto v_range = simplex.vertex_range();
    auto v_it = v_range.begin();
    for (; num_line < snap_d && v_it != v_range.end(); ++v_it, ++i) {
      if (i == snapping_indices[num_line]) {
	Eigen::VectorXd v_coords = triangulation.cartesian_coordinates(*v_it);
	for (std::size_t j = 0; j < amb_d; ++j)
	  vertex_matrix(num_line, j) = v_coords(j);
	reduced_lambda(num_line) = lambda(i);
	num_line++;
      }
    }
    reduced_lambda /= reduced_lambda.sum();
    Eigen::VectorXd intersection = reduced_lambda.transpose()*vertex_matrix;
    return QR({face_from_indices(simplex, snapping_indices), intersection, true});
  }
  
public:

  /** \brief Ambient dimension of the implicit manifold. */
  std::size_t amb_d() const {
    return fun_->amb_d();
  }
  
  /** \brief Codimension of the implicit manifold. */
  std::size_t cod_d() const {
    return fun_->cod_d();
  }

  /** \brief Intersection query with the relative interior of the manifold.
   *  
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the relative interior of the manifold exists, the intersection point
   *   and the face of the query simplex that contains 
   *   the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the same as the codimension of the manifold 
   *   (the codomain dimension of the function).
   *  @param[in] triangulation The ambient triangulation. The dimension of 
   *   the triangulation should be the same as the ambient dimension of the manifold 
   *   (the domain dimension of the function).
   */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects(const Simplex_handle& simplex,
					  const Triangulation& triangulation) const {
    Eigen::VectorXd lambda = compute_lambda(simplex, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }

  /** \brief Intersection query with the boundary of the manifold.
   *  
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the boundary of the manifold exists, the intersection point
   *   and the face of the query simplex that contains 
   *   the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the same as the codimension of the boundary of the manifold 
   *   (the codomain dimension of the function + 1).
   *  @param[in] triangulation The ambient triangulation. The dimension of 
   *   the triangulation should be the same as the ambient dimension of the manifold 
   *   (the domain dimension of the function).
   */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects_boundary(const Simplex_handle& simplex,
						   const Triangulation& triangulation) const {
    Eigen::VectorXd lambda = compute_boundary_lambda(simplex, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }

  
  /** \brief Returns true if the input point lies inside the piecewise-linear
   *   domain induced by the given ambient triangulation.
   *
   * @param p The input point. Needs to have the same dimension as the ambient
   *  dimension of the manifold (the domain dimension of the function).
   * @param triangulation The ambient triangulation. Needs to have the same
   *  dimension as the ambient dimension of the manifold 
   *  (the domain dimension of the function).
   */
  template <class Triangulation>
  bool lies_in_domain(const Eigen::VectorXd& p,
		      const Triangulation& triangulation) {
    return make_pl_approx(domain_fun_, triangulation)(p)(0) < 0;
  }
  
  Implicit_manifold_intersection_oracle(Function* function,
					Function* domain_function,
					double threshold = 0)
    : fun_(function), domain_fun_(domain_function->clone()), threshold_(threshold) {}

  Implicit_manifold_intersection_oracle(Function* function, double threshold = 0)
    : fun_(function),
      domain_fun_(new Constant_function(function->amb_d(), 1, Eigen::VectorXd::Constant(1,-1))),
      threshold_(threshold) {}

  ~Implicit_manifold_intersection_oracle() {
    delete domain_fun_;
  }
  
private:
  Function *fun_, *domain_fun_;
  double threshold_ = 0;
};


int main(int argc, char** const argv) {
  if (argc < 2) {
    std::cout << "Too few arguments. Usage: " << argv[0]
	      << " input_file_name [snapping_threshold] [triangulation_scale]\n";
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
  Function* function = *parse_input(stream, fun_range);
  stream.close();

  double threshold = 0;
  if (argc >= 3)
    threshold = atof(argv[2]);
  Implicit_manifold_intersection_oracle oracle(function, threshold);

  double lambda = 1;
  if (argc >= 4)
    lambda = atof(argv[3]);
  std::cout << "threshold = " << threshold << "\n";
  std::cout << "lambda = " << lambda << "\n";
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());

  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  Eigen::VectorXd seed;
  function->seed(seed);
  std::cout << "function->amb_d() = " << function->amb_d() << "\n";
  std::cout << "seed\n" << seed << "\n";
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, out_simplex_map);

  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cc(intr_d);
  cc.construct_complex(out_simplex_map);

  stream = std::ifstream(input_file_name);
  std::string output_file_name((std::istreambuf_iterator<char>(stream)),
			       std::istreambuf_iterator<char>());
  std::cout << "output_file_name = " << output_file_name << "\n";
  output_meshes_to_medit(3,
			 output_file_name,
			 build_mesh_from_cell_complex(cc,
						      Configuration({true, true, true, 1, 2, 3})));
  
  
  for (auto fun: fun_range)
    delete fun;

  return 0;
}
