#ifndef COXETER_TRIANGULATION_H_
#define COXETER_TRIANGULATION_H_

#include <stack>
#include <map>
#include <numeric> //iota

#include <boost/range/iterator_range.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include <gudhi/Coxeter_triangulation/Freudenthal_representation.h>
#include <gudhi/Coxeter_triangulation/Vertex_iterator.h>
#include <gudhi/Coxeter_triangulation/Face_iterator.h>
#include <gudhi/Coxeter_triangulation/Coface_iterator.h>

namespace Gudhi {
  
class Coxeter_triangulation : public Freudenthal_triangulation {

  typedef Eigen::MatrixXd Matrix;
  
public:
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors
  //////////////////////////////////////////////////////////////////////////////////////////////////

  Matrix root_matrix(unsigned d) {
    Matrix cartan(d,d);
    for (unsigned i = 0; i < d; i++) {
      cartan(i,i) = 1.0;
    }
    for (unsigned i = 1; i < d; i++) {
      cartan(i-1,i) = -0.5;
      cartan(i,i-1) = -0.5;
    }
    for (unsigned i = 0; i < d; i++)
      for (unsigned j = 0; j < d; j++)
	if (j+1 < i || j > i+1) 
	  cartan(i,j) = 0;    
    std::cout << "cartan =" << std::endl << cartan << std::endl;
    Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
    Eigen::VectorXd sqrt_diag(d);
    for (unsigned i = 0; i < d; ++i)
      sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
    std::cout << "sqrt_diag =" << std::endl << sqrt_diag << std::endl;

    Matrix lower(d,d);
    for (unsigned i = 0; i < d; i++)
      for (unsigned j = 0; j < d; j++)
	if (i < j)
	  lower(i,j) = 0;
	else
	  lower(i,j) = 1;
    Matrix result = lower*saes.eigenvectors()*sqrt_diag.asDiagonal();
    return result;
  }
  
  Coxeter_triangulation(unsigned dimension)
    : Freudenthal_triangulation(dimension, root_matrix(dimension)) {}
};

    
}

#endif
