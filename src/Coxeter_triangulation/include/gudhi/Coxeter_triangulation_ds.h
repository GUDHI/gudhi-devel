#ifndef COXETER_TRIANGULATION_DS_H_
#define COXETER_TRIANGULATION_DS_H_

#include <stack>

#include <gudhi/Coxeter_triangulation/Cell_id.h>

#include <boost/range/iterator_range.hpp>

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

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Cartesian coordinates and barycenter
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  Eigen::VectorXd cartesian_coordinates(const Cell_id& v_id) const {
    assert(v_id.dimension() == 0);
    Eigen::VectorXd val_vector(dimension_);
    std::size_t k = 0, j = 1;
    for (; j < (unsigned)dimension_ + 1; k += j, j++)
      val_vector(j-1) = v_id.value(k) / v_id.level();
    // return root_t_.colPivHouseholderQr().solve(val_vector);
    return colpivhouseholderqr_.solve(val_vector);
  }

  Eigen::VectorXd barycenter(const Cell_id& c_id) const {
    Eigen::VectorXd res_vector(dimension_);
    for (size_t i = 0; i < dimension_; ++i)
      res_vector(i) = 0;
    for (auto v: vertex_range(c_id)) {
      res_vector += cartesian_coordinates(v);
    }
    return (1/(c_id.dimension()+1)) * res_vector;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Vertex computation
  //////////////////////////////////////////////////////////////////////////////////////////////////

  class Vertex_iterator : public boost::iterator_facade< Vertex_iterator,
							 Cell_id const,
							 boost::forward_traversal_tag> {
  protected:
    friend class boost::iterator_core_access;

    bool equal(Vertex_iterator const& other) const {
      return (is_end_ && other.is_end_);
    }

    Cell_id const& dereference() const {
      return value_;
    }

    void update_value() {
      while (!is_end_) {
	std::size_t k = value_.size();
	std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
	if (state_stack_.top())
	  value_.push_back(c_id_.value(k)+1, true);
	else
	  value_.push_back(c_id_.value(k), true);
	std::size_t i = j-2;
	for (; i < j; --i) { // i is unsigned: if goes past 0, then becomes >j
	  std::size_t
	    k_s = i*(i+1)/2, // index of the simple root s_i
	    k_r = value_.size()-1; // index of the previous positive root r_{i+1,j}
	  value_.push_back(value_.value(k_r) + value_.value(k_s), true);
	  if (c_id_.mask(k_r+1)) {
	    if (value_.value(k_r+1) != c_id_.value(k_r+1))
	      break;
	  }
	  else // !c_id_.mask(k_r+1)
	    if (value_.value(k_r+1) < c_id_.value(k_r+1) ||
		value_.value(k_r+1) > c_id_.value(k_r+1) + 1)
	    break;
	}
	if (i >= j) { // success
	  if (value_.size() == c_id_.size()) // finished
	    return;
	  state_stack_.push(false);
	}
	if (i < j) { // fail
	  value_.resize(k);
	  elementary_increment();
	}
      }
    }

    void elementary_increment() {
      if (is_end_)
        return;
      while (!state_stack_.empty()) {
	bool p = state_stack_.top();
	state_stack_.pop();
	std::size_t j = state_stack_.size();
	std::size_t k_s = j*(j+1)/2;
	value_.resize(k_s);
	if (c_id_.mask(k_s))
	  continue;
	else if (p)
	  continue;
	else {
	  state_stack_.push(true);
	  return;
	}
      }
      if (state_stack_.empty())
	is_end_ = true;
    }
    
    void increment() {
      elementary_increment();
      update_value();
    }

  public:
    Vertex_iterator(const Cell_id& c_id,
		    const Coxeter_triangulation_ds& scs)
      :
      value_(c_id.level(), 0),
      is_end_(false),
      is_itself_(c_id.dimension() == 0)
    {
      if (is_itself_) {
	value_ = c_id_;
	return;
      }
      c_id_ = c_id;
      state_stack_.push(false);
      update_value();
    }

    Vertex_iterator() : is_end_(true) {}
  
  
  protected:
    Cell_id c_id_; 
    Cell_id value_;
    bool is_end_;
    bool is_itself_;
    std::stack<bool> state_stack_; // true means +1, false means +0
  };
  
  typedef boost::iterator_range<Vertex_iterator> Vertex_range;
  Vertex_range vertex_range(const Cell_id& a_id) const {
    return Vertex_range(Vertex_iterator(a_id, *this),
			Vertex_iterator());
  }  
  
protected:
  Matrix root_t_;
  unsigned dimension_;
  Eigen::ColPivHouseholderQR<Matrix> colpivhouseholderqr_;
};

    
}

#endif
