#ifndef IMPLICIT_FUNCTION_INTERSECTION_ORACLE_H_
#define IMPLICIT_FUNCTION_INTERSECTION_ORACLE_H_

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

namespace Gudhi {

template<class Function>
class Implicit_function_intersection_oracle {

public:  
  Implicit_function_intersection_oracle(const Function& fun)
    : fun_(fun) {}

  template <class Triangulation>
  bool intersects(const typename Triangulation::Simplex_handle& c,
		  const Triangulation& tr) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
    for (std::size_t i = 0; i < cod_d + 1; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: tr.vertex_range(c)) {
      Eigen::VectorXd v_coords = fun_(tr.cartesian_coordinates(v));
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      j++;
    }
    Eigen::VectorXd z(cod_d + 1);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 1; ++i)
      z(i) = 0;
    Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
    for (std::size_t i = 0; i < cod_d + 1; ++i)
      if (lambda(i) < 0 || lambda(i) > 1)
	return false;
    return true;
  }

  std::size_t amb_d() const {
    return fun_.amb_d();
  }
  
  std::size_t cod_d() const {
    return fun_.cod_d();
  }
  
private:
  Function fun_;
};

}

#endif
