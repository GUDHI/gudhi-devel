#ifndef FREUDENTHAL_TRIANGULATION_H_
#define FREUDENTHAL_TRIANGULATION_H_

#include <stack>
#include <map>
#include <numeric> //iota

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include <gudhi/Coxeter_triangulation/Freudenthal_representation.h>
#include <gudhi/Coxeter_triangulation/Vertex_iterator.h>
#include <gudhi/Coxeter_triangulation/Face_iterator.h>
#include <gudhi/Coxeter_triangulation/Coface_iterator.h>

namespace Gudhi {
  
class Freudenthal_triangulation {

  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;
  
public:

  typedef std::vector<int> Vertex_handle;
  typedef std::vector<std::vector<uint> > Ordered_partition_handle;
  typedef Freudenthal_representation<Vertex_handle, Ordered_partition_handle> Simplex_handle;
  typedef Gudhi::Vertex_iterator<Simplex_handle> Vertex_iterator;
  typedef Gudhi::Face_iterator<Simplex_handle> Face_iterator;
  typedef Gudhi::Coface_iterator<Simplex_handle> Coface_iterator;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors
  //////////////////////////////////////////////////////////////////////////////////////////////////

  
  Freudenthal_triangulation(unsigned dimension)
    : dimension_(dimension),
      matrix_(Matrix::Identity(dimension, dimension)),
      colpivhouseholderqr_(matrix_.colPivHouseholderQr()),
      is_freudenthal(true)  {
  }

  Freudenthal_triangulation(unsigned dimension, const Matrix& matrix)
    : dimension_(dimension),
      matrix_(matrix),
      colpivhouseholderqr_(matrix_.colPivHouseholderQr()),
      is_freudenthal(false) {
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Access functions
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  unsigned dimension() const {
    return dimension_;
  }

  Eigen::MatrixXd& matrix() {
    return matrix_;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Point location
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class Point_d>
  Simplex_handle locate_point(const Point_d& point, double scale=1) const {
    unsigned d = point.size();
    assert(d == dimension_);
    double error = 1e-9;
    Simplex_handle output;
    std::vector<double> z;
    if (is_freudenthal) {
      for (std::size_t i = 0; i < d; i++) {
	double x_i = scale * point[i];
	int y_i = std::floor(x_i);
	output.vertex.push_back(y_i);
	z.push_back(x_i - y_i);
      }
    }
    else {
      Eigen::VectorXd p_vect(d);
      for (std::size_t i = 0; i < d; i++)
	p_vect(i) = point[i];
      Eigen::VectorXd x_vect = colpivhouseholderqr_.solve(p_vect);
      for (std::size_t i = 0; i < d; i++) {
	double x_i = scale * x_vect(i);
	int y_i = std::floor(x_i);
	output.vertex.push_back(y_i);
	z.push_back(x_i - y_i);
      }
    }
    z.push_back(0);
    std::vector<std::size_t> indices(d+1);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&z](std::size_t i1, std::size_t i2) {return z[i1] > z[i2];});

    output.partition.push_back(std::vector<unsigned>(1, indices[0]));
    for (std::size_t i = 1; i <= d; ++i)
      if (z[indices[i-1]] > z[indices[i]] + error)
	output.partition.push_back(std::vector<unsigned>(1, indices[i]));
      else
	output.partition.back().push_back(indices[i]);
    return output;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Cartesian coordinates and barycenter
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  Eigen::VectorXd cartesian_coordinates(const Vertex_handle& v, double scale = 1) const {
    Eigen::VectorXd v_vect(dimension_);
    for (std::size_t j = 0; j < dimension_; j++)
      v_vect(j) = v[j] / scale;
    return matrix_ * v_vect;
  }

  Eigen::VectorXd barycenter(const Simplex_handle& simplex) const {
    Eigen::VectorXd res_vector(dimension_);
    for (size_t i = 0; i < dimension_; ++i)
      res_vector(i) = 0;
    for (auto v: vertex_range(simplex)) {
      res_vector += cartesian_coordinates(v);
    }
    return (1./(simplex.dimension()+1)) * res_vector;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Vertex computation
  //////////////////////////////////////////////////////////////////////////////////////////////////  
  typedef boost::iterator_range<Vertex_iterator> Vertex_range;
  Vertex_range vertex_range(const Simplex_handle& simplex) const {
    return Vertex_range(Vertex_iterator(simplex),
  			Vertex_iterator());
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Face computation
  //////////////////////////////////////////////////////////////////////////////////////////////////
  typedef boost::iterator_range<Face_iterator> Face_range;
  Face_range face_range(const Simplex_handle& simplex,
			std::size_t value_dim) const {
    return Face_range(Face_iterator(simplex, value_dim),
		      Face_iterator());
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Coface computation
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  typedef boost::iterator_range<Coface_iterator> Coface_range;
  Coface_range coface_range(const Simplex_handle& simplex,
			    std::size_t value_dim) const {
    return Coface_range(Coface_iterator(simplex, value_dim),
			Coface_iterator());
  }  

  
protected:
  unsigned dimension_;
  Matrix matrix_;
  Eigen::ColPivHouseholderQR<Matrix> colpivhouseholderqr_;
  bool is_freudenthal;
};

    
}

#endif
