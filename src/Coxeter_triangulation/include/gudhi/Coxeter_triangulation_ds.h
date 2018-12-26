#ifndef COXETER_TRIANGULATION_DS_H_
#define COXETER_TRIANGULATION_DS_H_

#include <gudhi/Coxeter_triangulation/Cell_id.h>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

namespace Gudhi {
  
class Coxeter_triangulation_ds {

  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;
  
public:

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors
  //////////////////////////////////////////////////////////////////////////////////////////////////

  
  Coxeter_triangulation_ds() {
  }

  Coxeter_triangulation_ds(unsigned dimension)
    : root_t_(dimension, dimension), dimension_(dimension){
    unsigned short d = dimension;
    Matrix cartan(d,d);
    for (unsigned i = 0; i < d; i++) {
      cartan(i,i) = 1.0;
    }
    for (unsigned i = 1; i < d; i++) {
      cartan(i-1,i) = -0.5;
      cartan(i,i-1) = -0.5;
    }
    for (int i = 0; i < d; i++)
      for (int j = 0; j < d; j++)
	if (j < i-1 || j > i+1) 
	  cartan(i,j) = 0;    
    // std::cout << "cartan =" << std::endl << cartan << std::endl;
    Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
    Eigen::VectorXd sqrt_diag(d);
    for (int i = 0; i < d; ++i)
      sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
    root_t_ = saes.eigenvectors()*sqrt_diag.asDiagonal();
    colpivhouseholderqr_ = root_t_.colPivHouseholderQr();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Access functions
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  unsigned dimension() const {
    return dimension_;
  }

  unsigned pos_root_count() const {
    return dimension_*(dimension_ + 1)/2;
  }

  Eigen::MatrixXd const& simple_root_matrix() const {
    return root_t_;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Alcove dimension
  //////////////////////////////////////////////////////////////////////////////////////////////////
  unsigned alcove_dimension(const Cell_id& a_id) const {
    std::size_t i = 0, j = 1, k = 0;
    unsigned return_value = dimension_;
    while (k < a_id.size()) {
      if (a_id.mask(k)) {
	std::size_t l = i + 1;
	bool lin_independent = true;
	while (l < j && lin_independent) {
	  std::size_t k1 = (l*l+l-2)/2 - i;
	  std::size_t k2 = (j*j+j-2)/2 - l;
	  lin_independent = (!a_id.mask(k1)) || (!a_id.mask(k2));
	  l++;
	}
	if (lin_independent) {
	  return_value--;
	  if (!return_value)
	    return 0;
	}
      }
      k++;
      if (i == 0) {
	j++;
	i = j - 1;
      }
      else
	i--;
    }
    return return_value;
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Point location
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class Point_d>
  Cell_id locate_point(const Point_d& point, double level=1) {
    unsigned short d = point.size();
    assert(d == dimension_);
    double error = 1e-9;
    Cell_id c_id(level, d);
    Eigen::VectorXd p_vect(d);
    for (short i = 0; i < d; i++)
      p_vect(i) = point[i];
    Eigen::VectorXd scalprod_vect = root_t_ * p_vect;
    for (short i = 0; i < d; i++) {
      double root_scalprod = 0;
      for (short j = i; j >= 0; j--) {
	root_scalprod += scalprod_vect(j);
	double value = level * root_scalprod;
	if (std::abs(value - std::round(value)) >= error)
	  c_id.push_back(std::floor(level * root_scalprod), false);
	else
	  c_id.push_back(std::round(level * root_scalprod), true);
      }
    }
    c_id.set_dimension(alcove_dimension(c_id));
    return c_id;
  }
  
protected:
  Matrix root_t_;
  unsigned dimension_;
  Eigen::ColPivHouseholderQR<Matrix> colpivhouseholderqr_;
};

    
}

#endif
