#ifndef SIMPLE_COXETER_SYSTEM_REMASTERED_H_
#define SIMPLE_COXETER_SYSTEM_REMASTERED_H_

#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <exception>

#include <gudhi/Coxeter_complex/Alcove_id.h>
#include <gudhi/Clock.h>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include "../../example/cxx-prettyprint/prettyprint.hpp"
// #include <Eigen/SPQRSupport>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpzf ET;
// typedef double ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

double prec = 1e-15;

class Simple_coxeter_system {

  typedef double FT;
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SparseMatrix<FT> SparseMatrix;
  typedef Eigen::Triplet<FT> Triplet;
    
  class wrong_family : public std::exception {  
  } wrong_family_exception_;
  
public:
  typedef Gudhi::Alcove_id Alcove_id;
  typedef Alcove_id Vertex_id;
  struct Filtered_alcove {
    Alcove_id id;
    double f;
    Filtered_alcove(const Alcove_id& id_in, double f_in) : id(id_in), f(f_in) {}
    Filtered_alcove(const Alcove_id& id_in) : id(id_in), f(0) {}
  };
  
public:
  
  Matrix root_t_;
  char family_;
  unsigned short dimension_;
  unsigned short vertex_level_ = 1;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  Simple_coxeter_system() {
  }
  
  Simple_coxeter_system(char family, unsigned short dimension)
    : root_t_(dimension, dimension), family_(family), dimension_(dimension){
    switch (family) {
    case 'A': {
      unsigned short d = dimension;
      vertex_level_ = 1;
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
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }
    case 'B': {
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 2);
      Matrix cartan(d,d);
      for (int i = 0; i < d-1; i++) {
        cartan(i,i) = 4.0;
      }
      cartan(d-1,d-1) = 2.0;
      for (int i = 1; i < d; i++) {
        cartan(i-1,i) = -2.0;
        cartan(i,i-1) = -2.0;
      }
      Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
      Eigen::VectorXd sqrt_diag(d);
      for (int i = 0; i < d; ++i)
        sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
      root_t_ = saes.eigenvectors()*sqrt_diag.asDiagonal();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'C': {
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 2);
      Matrix cartan(d,d);
      for (int i = 0; i < d-1; i++) {
        cartan(i,i) = 2.0;
      }
      cartan(d-1,d-1) = 4.0;
      for (int i = 1; i < d-1; i++) {
        cartan(i-1,i) = -1.0;
        cartan(i,i-1) = -1.0;
      }
      cartan(d-1,d-2) = -2.0;
      cartan(d-2,d-1) = -2.0;
      Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
      Eigen::VectorXd sqrt_diag(d);
      for (int i = 0; i < d; ++i)
        sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
      root_t_ = saes.eigenvectors()*sqrt_diag.asDiagonal();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'D': {   
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 3);
      Matrix cartan(d,d);
      for (int i = 0; i < d; i++) {
        cartan(i,i) = 2.0;
      }
      for (int i = 1; i < d-1; i++) {
        cartan(i-1,i) = -1.0;
        cartan(i,i-1) = -1.0;
      }
      cartan(d-1,d-3) = -1.0;
      cartan(d-3,d-1) = -1.0;
      Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
      Eigen::VectorXd sqrt_diag(d);
      for (int i = 0; i < d; ++i)
        sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
      root_t_ = saes.eigenvectors()*sqrt_diag.asDiagonal();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'E': {   
      // unsigned short d = dimension;
      // vertex_level_ = 2;
      // assert(d >= 6 && d <= 8);
      // std::vector<Triplet> cartan_triplets;
      // cartan_triplets.reserve(3*d-2);
      // for (unsigned i = 0; i < d; i++) {
      //   cartan_triplets.push_back(Triplet(i,i,2.0));
      // }
      // for (int i = 3; i < d; i++) {
      //   cartan_triplets.push_back(Triplet(i-1,i,-1.0));
      //   cartan_triplets.push_back(Triplet(i,i-1,-1.0));
      // }
      // cartan_triplets.push_back(Triplet(2,0,-1.0));
      // cartan_triplets.push_back(Triplet(3,1,-1.0));
      // cartan_triplets.push_back(Triplet(0,2,-1.0));
      // cartan_triplets.push_back(Triplet(1,3,-1.0));
      // Matrix cartan(d,d);
      // cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      // std::cout << "cartan =" << std::endl << cartan << std::endl;
      // Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      // root_t_ = chol.matrixL();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      // Eigen::MatrixXf base(8,8);      
      // std::vector<Triplet> base_triplets;
      // base(0,0) = 0.5;
      // base(7,0) = 0.5;
      // base(1,0) = -0.5;
      // base(2,0) = -0.5;
      // base(3,0) = -0.5;
      // base(4,0) = -0.5;
      // base(5,0) = -0.5;
      // base(6,0) = -0.5;
      // base(0,1) = 1.0;
      // base(1,1) = 1.0;
      // for (int i = 0; i < d-2; ++i) {
      //   base(i, i+2) =  -1.0;
      //   base(i+1, i+2) =  1.0;
      // }
      // if (d <= 6) {
      //   base(5,6) = 1.0;
      //   base(6,6) = -1.0;
      // }
      // if (d <= 7) {
      //   base(6,7) = 1.0;
      //   base(7,7) = 1.0;
      // }
      // std::cout << "base = " << std::endl << base << std::endl;
      // std::cout << "base^{-1} = " << std::endl << base.inverse() << std::endl;
      // for (int i = )
      break;
    }  
    default :
      std::cerr << "Simple_coxeter_system : The family " << family << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Access functions
  //////////////////////////////////////////////////////////////////////////////////////////////////

  
  unsigned short dimension() const {
    return dimension_;
  }

  char family() const {
    return family_;
  }

  unsigned pos_root_count() const {
    switch (family_) {
    case 'A': { return dimension_*(dimension_ + 1)/2; break; } 
    case 'B': { return dimension_*dimension_; break; } 
    case 'C': { return dimension_*dimension_; break; } 
    case 'D': { return dimension_*(dimension_ - 1); break; }
    case 'E': {
      switch (dimension_) {
      case 6: return 36;
      case 7: return 63;
      case 8: return 120;
      }
    }
    default :
      std::cerr << "Simple_coxeter_system::pos_root_count : The family " << family_ << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
  }

  Eigen::MatrixXd const& simple_root_matrix() const {
    return root_t_;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Alcove dimension
  //////////////////////////////////////////////////////////////////////////////////////////////////
  unsigned short alcove_dimension(const Alcove_id& a_id) const {
    switch (family_) {
    case 'A': {
      std::size_t i = 0, j = 1, k = 0;
      unsigned short return_value = dimension_;
      while (k < a_id.size()) {
        if (a_id.is_fixed(k)) {
          std::size_t l = i + 1;
          bool lin_independent = true;
          while (l < j && lin_independent) {
            std::size_t k1 = (l*l+l-2)/2 - i;
            std::size_t k2 = (j*j+j-2)/2 - l;
            lin_independent = (!a_id.is_fixed(k1)) || (!a_id.is_fixed(k2));
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
    default:
     std::cerr << "Simple_coxeter_system::dimension : The family " << family_ << " is not supported. "
                << "Please use A family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    } 
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Normal base iterator
  //////////////////////////////////////////////////////////////////////////////////////////////////

public:
  class Normal_basis_iterator : public boost::iterator_facade< Normal_basis_iterator,
                                                               std::size_t const,
                                                               boost::forward_traversal_tag > {
  protected:
    friend class boost::iterator_core_access;
    
    void update_value() {
      if (j_ != d_ + 1 && !*it_)
        elementary_increment();
    }
    
    bool equal(Normal_basis_iterator const& other) const {
      return value_ == other.value_;
    }

    void elementary_increment() {
      value_++;
      if (i_ == 0) {
        i_ = j_;
        j_++;
      }
      else
        i_--;
      it_++;
      update_value();
    }
    
    void increment() {
      while (i_ != 0) {
        value_++;
        i_--;
        it_++;
      }
      i_ = j_;
      j_++;
      value_++;
      it_++;
      update_value();
    }

    std::size_t const& dereference() const {
      return value_;
    }
    
  public:
    Normal_basis_iterator(const Alcove_id& a_id,
                          const Simple_coxeter_system& scs)
      : i_(0),
        j_(1),
        value_(0),
        d_(scs.dimension()),
        it_(a_id.fixed_begin())
    {
      update_value();
    }

    Normal_basis_iterator(const Simple_coxeter_system& scs)
      : value_(scs.pos_root_count()) { }
        
    std::size_t i_, j_, value_, d_;
    typename Alcove_id::fixed_const_iterator it_;
  };

  typedef boost::iterator_range<Normal_basis_iterator> Normal_basis_range;
  Normal_basis_range normal_basis_range(const Alcove_id& a_id) const {
    return Normal_basis_range(Normal_basis_iterator(a_id, *this),
                              Normal_basis_iterator(*this));
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Query point location
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point,
            class OutputIterator>
  void query_point_location(const Point& p, double level, OutputIterator output_it) const
  {
    unsigned short d = p.size();
    assert(d == dimension_);
    Eigen::VectorXd p_vect(d);
    for (short i = 0; i < d; i++)
      p_vect(i) = p[i];
    Eigen::VectorXd scalprod_vect = root_t_ * p_vect;
    switch (family_) {
    case 'A': {
      // e_i - e_j
      for (short i = 0; i < d; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      break;
    }
    case 'B': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      // e_i
      FT root_scalprod = 0;
      for (short i = d-1; i >= 0; i--) {
        root_scalprod += scalprod_vect(i);
        *output_it++ = std::floor(level * root_scalprod);
      }
      // e_i + e_j
      FT global_scalprod = 0;
      for (short i = d-1; i >= 0; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      break;
    }
    case 'C': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      // 2*e_i
      FT root_scalprod = -scalprod_vect(d-1);
      for (short i = d-1; i >= 0; i--) {
        root_scalprod += 2*scalprod_vect(i);
        *output_it++ = std::floor(level * root_scalprod);
      }
      // e_i + e_j
      FT global_scalprod = -scalprod_vect(d-1);
      for (short i = d-1; i >= 0; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      break;
    }
    case 'D': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      // e_i + e_j
      FT global_scalprod = -scalprod_vect(d-1) - scalprod_vect(d-2);
      for (short i = d-1; i >= 1; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      break;
    }
      
    case 'E': {
      // r1
      *output_it++ = std::floor(level * scalprod_vect(0));
      // // r2
      // *output_it++ = std::floor(level * scalprod_vect(1));
      for (short k = 1; k < d; k++) {
        // e_i - e_j OK
        FT root_scalprod = 0;
        for (short j = k; j >= 3; j--) {
          root_scalprod += scalprod_vect(j);
          *output_it++ = std::floor(level * root_scalprod);
        }
        // e_i + e_j OK
        FT global_scalprod = -scalprod_vect(1)-scalprod_vect(2);
        for (short i = 2; i <= k; i++) {
          global_scalprod += 2*scalprod_vect(i);
          FT root_scalprod = global_scalprod;
          for (short j = i+1; j <= k; j++) {
            root_scalprod += scalprod_vect(j);
          }
          *output_it++ = std::floor(level * root_scalprod);
        }
        // 1/2(prefix + sum of 1 - sum of (2d-8)) OK
        root_scalprod = scalprod_vect(0);
        for (short i = k; i >= 3; i--)
          root_scalprod += scalprod_vect(i);
        *output_it++ = std::floor(level * root_scalprod);
        // 1/2(prefix + sum of 3 - sum of (2d-10))
        global_scalprod = scalprod_vect(0) + scalprod_vect(3) + scalprod_vect(4) + scalprod_vect(5);
        for (short i = 1; i < k-1; i++) {
          global_scalprod += scalprod_vect(i) + scalprod_vect(i+1);
          *output_it++ = std::floor(level * global_scalprod);
          root_scalprod = global_scalprod;
          for (short j = i+2; j < k; j++) {
            root_scalprod += scalprod_vect(j);
            *output_it++ = std::floor(level * root_scalprod);
          }
        }
        // 1/2(prefix + sum of 5 - sum of (2d-12))
        for (short i = 1; i < k-4; i++) {
          root_scalprod = scalprod_vect(0) + 2*scalprod_vect(1) + 2*scalprod_vect(2) + 3*scalprod_vect(3) + 2*scalprod_vect(4) + scalprod_vect(5);
          *output_it++ = std::floor(level * root_scalprod);
        }
        // 1/2(prefix + sum of 7 - sum of (2d-14))
        if (d == 8) {
          root_scalprod = scalprod_vect(0) + 2*scalprod_vect(1) + 2*scalprod_vect(2) + 3*scalprod_vect(3) + 2*scalprod_vect(4) + scalprod_vect(5);
          *output_it++ = std::floor(level * root_scalprod);
        }
      }
      break;
    }
    default :
      std::cerr << "Simple_coxeter_system::alcove_coordinates : The family " << family_ << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
  }

  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point>
  Alcove_id query_point_location(const Point& p, double level) const {
    unsigned short d = p.size();
    assert(d == dimension_);
    Alcove_id a_id(level, d);
    Eigen::VectorXd p_vect(d);
    for (short i = 0; i < d; i++)
      p_vect(i) = p[i];
    Eigen::VectorXd scalprod_vect = root_t_ * p_vect;
    switch (family_) {
    case 'A': {
      // e_i - e_j
      for (short i = 0; i < d; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      break;
    }
    case 'B': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      // e_i
      FT root_scalprod = 0;
      for (short i = d-1; i >= 0; i--) {
        root_scalprod += scalprod_vect(i);
        a_id.push_back(std::floor(level * root_scalprod), false);
      }
      // e_i + e_j
      FT global_scalprod = 0;
      for (short i = d-1; i >= 0; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      break;
    }
    case 'C': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      // 2*e_i
      FT root_scalprod = -scalprod_vect(d-1);
      for (short i = d-1; i >= 0; i--) {
        root_scalprod += 2*scalprod_vect(i);
        a_id.push_back(std::floor(level * root_scalprod), false);
      }
      // e_i + e_j
      FT global_scalprod = -scalprod_vect(d-1);
      for (short i = d-1; i >= 0; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      break;
    }
    case 'D': {
      // e_i - e_j
      for (short i = 0; i < d-1; i++) {
        FT root_scalprod = 0;
        for (short j = i; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      // e_i + e_j
      FT global_scalprod = -scalprod_vect(d-1) - scalprod_vect(d-2);
      for (short i = d-1; i >= 1; i--) {
        global_scalprod += 2*scalprod_vect(i);
        FT root_scalprod = global_scalprod;
        for (short j = i-1; j >= 0; j--) {
          root_scalprod += scalprod_vect(j);
          a_id.push_back(std::floor(level * root_scalprod), false);
        }
      }
      break;
    }
    default :
      std::cerr << "Simple_coxeter_system::alcove_coordinates : The family " << family_ << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
    return a_id;
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Distance from a point to an alcove
  //////////////////////////////////////////////////////////////////////////////////////////////////
  template <class Point>
  double squared_distance_to_alcove(const Point& p,
                                    const Alcove_id& a_id) const {
    switch (family_) {
    case 'A': {
      double level = a_id.level();
      std::size_t d = dimension_;
      double sq_distance = 0;
      Eigen::VectorXd p_vect(d);
      for (unsigned short i = 0; i < d; i++)
        p_vect(i) = p[i];
      // std::cout << "a_id =" << a_id << "\n";
      // std::cout << "root_t_=\n" << root_t_ << "\n";
      std::list<std::pair<Eigen::VectorXd, double> > ah_hplanes, f_hplanes;
      // Find a basis for the affine hull of a_id
      std::size_t i = 0, j = 1, k = 0;
      while (k < a_id.size()) {
        if (a_id.is_fixed(k)) {
          std::size_t l = i + 1;
          bool lin_independent = true;
          while (l < j && lin_independent) {
            std::size_t k1 = (l*l+l-2)/2 - i;
            std::size_t k2 = (j*j+j-2)/2 - l;
            lin_independent = (!a_id.is_fixed(k1)) || (!a_id.is_fixed(k2));
            l++;
          }
          if (lin_independent) {
            ah_hplanes.push_back(std::make_pair(Eigen::VectorXd(d), a_id[k]/level));
            for (std::size_t m = 0; m < d; ++m)
              ah_hplanes.back().first(m) = 0;
            for (std::size_t l = i; l < j; ++l)
              for (std::size_t m = 0; m < d; ++m)
                ah_hplanes.back().first(m) += root_t_(l,m);
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
      // std::cout << "AH Hyperplanes:\n";
      // for (auto p: ah_hplanes) {
      //   std::cout << "n = " << p.first << "\na = " << p.second << "\n";
      // }
      // Add the hyperplanes that define the facets
      if (a_id.dimension() != 0)
        for (auto f_id: face_range(a_id, a_id.dimension()-1, a_id.dimension()-1))
          for (i = 0, j = 1, k = 0; k < a_id.size(); ++k) {
            if (f_id.is_fixed(k) && !a_id.is_fixed(k)) {
              // std::cout << "Facet : " << f_id << "\n";
              if (f_id[k] == a_id[k] + 1) {
                f_hplanes.push_back(std::make_pair(Eigen::VectorXd(d), (a_id[k] + 1)/level));
                for (std::size_t m = 0; m < d; ++m)
                  f_hplanes.back().first(m) = 0;
                for (std::size_t l = i; l < j; ++l)
                  for (std::size_t m = 0; m < d; ++m)
                    f_hplanes.back().first(m) += root_t_(l,m);
              }
              else {
                f_hplanes.push_back(std::make_pair(Eigen::VectorXd(d), -a_id[k]/level));
                for (std::size_t m = 0; m < d; ++m)
                  f_hplanes.back().first(m) = 0;
                for (std::size_t l = i; l < j; ++l)
                  for (std::size_t m = 0; m < d; ++m)
                    f_hplanes.back().first(m) -= root_t_(l,m);
              }
              break;
            }
            if (i == 0) {
              j++;
              i = j - 1;
            }
            else
              i--;
          }
      // std::cout << "F Hyperplanes:\n";
      // for (auto p: f_hplanes) {
      //   std::cout << "n = " << p.first << "\na = " << p.second << "\n";
      // }
      // Iteratively project the point on the affine space of the simplex
      for (auto ah_hplane_it = ah_hplanes.begin();
           ah_hplane_it != ah_hplanes.end();
           ) {
        Eigen::VectorXd& n1 = ah_hplane_it->first;
        double& a1 = ah_hplane_it->second;
        // std::cout << "Iteration for n =\n" << n1 << "\n a = " << a1 << "\n";
        double diff = a1 - p_vect.dot(n1);
        p_vect = p_vect + diff * n1;
        // std::cout << "diff = " << diff << "\n p =\n" << p_vect << "\n";
        sq_distance += diff*diff;
        // std::cout << "sq_dist = " << sq_distance << "\n";
        for (auto ah_hplane_after_it = ++ah_hplane_it;
             ah_hplane_after_it != ah_hplanes.end();
             ++ah_hplane_after_it) {
          Eigen::VectorXd& n2 = ah_hplane_after_it->first;
          double a2 = ah_hplane_after_it->second;
          double n1n2 = n1.dot(n2);
          n2 = n2 - n1n2 * n1;
          double inv_newnorm = 1/n2.norm();
          n2 = inv_newnorm * n2;
          a2 = inv_newnorm * (a2 - n1n2 * a1);
        }
        // std::cout << "AH Hyperplanes\n";
        // for (auto p: ah_hplanes) {
        //   std::cout << "n = " << p.first << "\na = " << p.second << "\n";
        // }
        for (auto f_hplane_it = f_hplanes.begin();
             f_hplane_it != f_hplanes.end();
             ++f_hplane_it) {
          Eigen::VectorXd& n2 = f_hplane_it->first;
          double a2 = f_hplane_it->second;
          double n1n2 = n1.dot(n2);
          n2 = n2 - n1n2 * n1;
          double inv_newnorm = 1/n2.norm();
          n2 = inv_newnorm * n2;
          a2 = inv_newnorm * (a2 - n1n2 * a1);
        }
        // std::cout << "F Hyperplanes\n";
        // for (auto p: f_hplanes) {
        //   std::cout << "n = " << p.first << "\na = " << p.second << "\n";
        // }
      }
      // Iteratively project the point on the hyperplanes that define the facets
      auto f_hplane_it = f_hplanes.begin();
      // std::cout << "p_vect:\n" << p_vect << "\n";
      while (f_hplane_it != f_hplanes.end()) {
        Eigen::VectorXd& n1 = f_hplane_it->first;
        double& a1 = f_hplane_it->second;
        // std::cout << "n1:\n" << n1 << "\n";
        // std::cout << "a1 = " << a1 << "\n";
        // std::cout << "p_vect.dot(n1) = " << p_vect.dot(n1) << "\n";
        // std::cout << "diff = " << a1 - p_vect.dot(n1) << "\n";
        double diff = a1 - p_vect.dot(n1);
        if (diff < 0) {
          p_vect = p_vect + diff * n1;
          sq_distance += diff*diff;
          // std::cout << "p_vect:\n" << p_vect << "\n";
          // std::cout << "sq_distance = " << sq_distance << "\n";
          for (auto f_hplane_other_it = f_hplanes.begin();
               f_hplane_other_it != f_hplanes.end();
               ++f_hplane_other_it) {
            if (f_hplane_other_it != f_hplane_it) {
              Eigen::VectorXd& n2 = f_hplane_other_it->first;
              double& a2 = f_hplane_other_it->second;
              double n1n2 = n1.dot(n2);
              n2 = n2 - n1n2 * n1;
              double inv_newnorm = 1./n2.norm();
              n2 = inv_newnorm * n2;
              // std::cout << "inv_newnorm = " << inv_newnorm << "\n";
              // std::cout << "a2 = " << a2 << "\n";
              // std::cout << "n1n2 = " << n1n2 << "\n";
              // std::cout << "a1 = " << a1 << "\n";
              a2 = inv_newnorm * (a2 - n1n2 * a1);
              // std::cout << "result a2 = " << a2 << "\n\n";
            }
          }
          f_hplanes.erase(f_hplane_it);
          f_hplane_it = f_hplanes.begin();
        }
        else
          f_hplane_it++;
      }
      // std::cout << "sq_dist = " << sq_distance << "\n";
      return sq_distance;
    }
    default :
      std::cerr << "Simple_coxeter_system::alcove_coordinates : The family " << family_ << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Query ball intersection
  //////////////////////////////////////////////////////////////////////////////////////////////////

  /** A function that outputs the alcoves that intersect the ball centered 
   *  at p of radius eps.
   *  The output pairs consist of the Alcove_id and the squared distance.
   */
  template <class Point>
  void alcoves_of_ball(const Point& p,
                       double level,
                       double eps,
                       std::vector<Filtered_alcove>& alcoves,
                       std::vector<std::vector<Vertex_id> >& vertices_per_alcove,
                       bool root_coords = false) const {
    unsigned d = dimension_;
    Alcove_id a_id(level, d);
    a_id.reserve(pos_root_count());
    Eigen::VectorXd p_vect(d);
    for (unsigned short i = 0; i < d; i++)
      p_vect(i) = p[i];
    Eigen::VectorXd scalprod_vect;
    if (root_coords)
      scalprod_vect = p_vect;
    else
      scalprod_vect = root_t_ * p_vect;
#ifdef CC_STAR_COMPLETION
    rec_alcoves_of_ball_A(a_id, scalprod_vect, eps, alcoves, vertices_per_alcove, 1, 0, 0, 0, p_vect);
#else
    rec_alcoves_of_ball_A(p, a_id, scalprod_vect, eps, alcoves, vertices_per_alcove, 1, 0, 0, p_vect);
#endif
  }
  
private:

  /** Construct the simplices that intersect a given ball.
   */
#ifdef CC_STAR_COMPLETION
  void rec_alcoves_of_ball_A(Alcove_id& a_id,
                             Eigen::VectorXd& scalprod_vect,
                             double eps,
                             std::vector<Filtered_alcove>& alcoves,
                             std::vector<std::vector<Vertex_id> >& vertices_per_alcove,
                             int j,
                             int i,
                             double root_scalprod,
                             double filtration,
                             Eigen::VectorXd& p_vect) const {
#else
  template <class Point>
  void rec_alcoves_of_ball_A(const Point& p,
                             Alcove_id& a_id,
                             Eigen::VectorXd& scalprod_vect,
                             double eps,
                             std::vector<Filtered_alcove>& alcoves,
                             std::vector<std::vector<Vertex_id> >& vertices_per_alcove,
                             int j,
                             int i,
                             double root_scalprod,
                             Eigen::VectorXd& p_vect) const {
#endif
    unsigned short d = dimension_;
    double level = a_id.level();
    double root_norm = level;
    if (j == d+1) {
#ifdef CC_STAR_COMPLETION
      alcoves.emplace_back(Filtered_alcove(a_id, filtration));
#else
      double sq_norm = squared_distance_to_alcove(p, a_id);
      sq_norm = std::round(sq_norm*10e10)/10e10;
      if (sq_norm <= eps*eps) {
        alcoves.emplace_back(Filtered_alcove(a_id, sq_norm));
#ifdef CC_A_V_VISITORS
        vertices_per_alcove.push_back(vertices);
#endif
      }
#endif

      return;
    }
    if (i == -1) {
#ifdef CC_STAR_COMPLETION
      rec_alcoves_of_ball_A(a_id,
                            scalprod_vect,
                            eps,
                            alcoves,
                            vertices_per_alcove,
                            j+1,
                            j,
                            0,
                            filtration,
                            p_vect);
#else
      rec_alcoves_of_ball_A(p,
                            a_id,
                            scalprod_vect,
                            eps,
                            alcoves,
                            vertices_per_alcove,
                            j+1,
                            j,
                            0,
                            p_vect);
#endif
      return;
    }
    root_scalprod += scalprod_vect(i);
    int min_lim = std::floor(level*root_scalprod - root_norm*eps);
    int max_lim = std::floor(level*root_scalprod + root_norm*eps);
    for (int val = min_lim; val <= max_lim; ++val) {
      bool valid = true;
      for (int l = i+1; l < j; ++l) {
        int pr_i_l = *(a_id.end() - ((j+l-1)*(j-l)/2+(j-l))),
          pr_l_j = *(a_id.end() - (l-i));
        // check if the floor(<x,r>) + floor(<x,s>) - 1 <= floor(<x,r+s>) <= floor(<x,r>) + floor(<x,s>)
        // for a point in the ball B(p,eps)
        if (val < pr_i_l + pr_l_j || val > pr_i_l + pr_l_j + 1) {
          valid = false;
          break;
        }
      }
      if (valid) {
        a_id.push_back(val);
#ifdef CC_STAR_COMPLETION
        double new_filtration = 0;
        int true_value = std::floor(level*root_scalprod);
        if (val > true_value) {
          new_filtration = (val - level*root_scalprod)/root_norm;        
        }
        else if (val < true_value)
          new_filtration = (level*root_scalprod  - val - 1)/root_norm;
        // std::cout << "val=" << val << ", true_value=" << true_value << ", filtration=" << new_filtration << "\n";
        rec_alcoves_of_ball_A(a_id,
                              scalprod_vect,
                              eps,
                              alcoves,
                              vertices_per_alcove,
                              j,
                              i-1,
                              root_scalprod,
                              std::max(filtration, new_filtration),
                              p_vect);
#else
        rec_alcoves_of_ball_A(p,
                              a_id,
                              scalprod_vect,
                              eps,
                              alcoves,
                              vertices_per_alcove,
                              j,
                              i-1,
                              root_scalprod,
                              p_vect);
#endif
        a_id.pop_back();
      }
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Face range
  //////////////////////////////////////////////////////////////////////////////////////////////////

  public:
  class Face_iterator : public boost::iterator_facade< Face_iterator,
                                                       Alcove_id const,
                                                       boost::forward_traversal_tag> {
  protected:
    // using State_pair = std::pair<std::size_t, bool>;
    friend class boost::iterator_core_access;

    bool triplet_check(std::size_t i, std::size_t l, std::size_t j) const {
      int k  = (j*j+j-2)/2 - i;
      int k1 = (l*l+l-2)/2 - i;
      int k2 = (j*j+j-2)/2 - l;
      if (value_.is_fixed(k))
        return (value_.is_fixed(k1) && value_.is_fixed(k2) && value_[k1] + value_[k2] == value_[k]) ||
               (!value_.is_fixed(k1) && !value_.is_fixed(k2) && value_[k1] + value_[k2] + 1 == value_[k]);
      else
        return ((value_.is_fixed(k1) xor value_.is_fixed(k2)) && value_[k] == value_[k1] + value_[k2]) ||
               (!value_.is_fixed(k1) && !value_.is_fixed(k2) && (value_[k] == value_[k1] + value_[k2] ||
                                                                 value_[k] == value_[k1] + value_[k2] + 1));
    }
    
    void update_value() {
      if (is_end_)
        return;
      switch (family_) {
      case 'A': {
        // current position variables
        while (!is_end_) {
          std::size_t k = value_.size();
          std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
          std::size_t i = (j*j + j - 2)/2 - k;
          if (stack_.top().f)
            if (stack_.top().p)
              value_.push_back(a_id_[k]+1, true);
            else
              value_.push_back(a_id_[k], true);
          else
            if (stack_.top().p)
              value_.push_back(a_id_[k]-1, false);
            else
              value_.push_back(a_id_[k], false);
          bool curr_state_is_valid = true;
          if (basis_k_.size() == j && i != basis_k_.back())
            for (std::size_t l = basis_k_.back(); l < j && curr_state_is_valid; ++l)
              curr_state_is_valid = triplet_check(i,l,j);
          else
            for (std::size_t l = i + 1; l < j && curr_state_is_valid; ++l)
              curr_state_is_valid = triplet_check(i,l,j);
          if (curr_state_is_valid &&
              value_.size() == a_id_.size()) {
            value_.set_dimension(curr_dim_lower_);
            return;
          }
          if (!curr_state_is_valid) {
            elementary_increment();
            continue;
          }
          if (a_id_.is_fixed(k+1)) {
            if (!stack_push(true, false)) {
              stack_pop();
              elementary_increment();
            }
          }
          else 
            if (!stack_push(false, false)) {
              stack_pop();
              if (!stack_push(true, false)) {
                stack_pop();
                elementary_increment();
              }
            }
        }
        break;
      }
      default:
        std::cerr << "Simple_coxeter_system::facet_iterator : The family " << family_ << " is not supported. "
                  << "Please use A family for the constructor (in capital).\n";
      }
    }
    
    bool equal(Face_iterator const& other) const {
      return (is_end_ && other.is_end_) || (stack_ == other.stack_);
    }
    Alcove_id const& dereference() const {
      return value_;
    }
    void elementary_increment() {
      if (is_end_)
        return;
      while (!stack_.empty()) {
        bool f = stack_.top().f;
        bool p = stack_.top().p;
        stack_pop();
        value_.pop_back();
        if (a_id_.is_fixed(stack_.size()))
          continue;
        else
          if (f)
            if (p)
              continue;
            else
              if (stack_push(true, true))
                break;
              else {
                stack_pop();
                continue;
              }
          else
            if (stack_push(true, false)) {
              break;
            }
            else {
              stack_pop();
              continue;
            }
      }
      if (stack_.empty())
        is_end_ = true;
    }
    
    void increment() {
      elementary_increment();
      update_value();
    }

    bool stack_push(bool fixed, bool plus_one) {
      std::size_t k = stack_.size();
      std::size_t j = std::floor(0.5*(1+std::sqrt(1+8*k)));
      std::size_t i = (j*j+j-2)/2 - k;
      if (basis_k_.size() == j - 1) {
        bool is_basis_vect = true;
        for (std::size_t l = i + 1; l < j && is_basis_vect; ++l) {
          std::size_t k1 = (l*l+l-2)/2 - i;
          is_basis_vect = !value_.is_fixed(k1);
        }
        if (!fixed)
          for (std::size_t l = 0; l < i && is_basis_vect; ++l) {
           is_basis_vect = false;
            for (std::size_t m = l+1; m < j && !is_basis_vect; ++m) {
              std::size_t k1 = (m*m+m-2)/2 - l;
              is_basis_vect = value_.is_fixed(k1);
            }
          }
        if (is_basis_vect) {
          basis_k_.push_back(i);
          if (fixed)
            curr_dim_upper_--;
          else
            curr_dim_lower_++;
        }
      }
      stack_.push({fixed, plus_one});        
      return curr_dim_lower_ <= curr_dim_upper_ &&
             curr_dim_lower_ <= value_dim_upper_ &&
             curr_dim_upper_ >= value_dim_lower_;
    }

    void stack_pop() {
      bool f = stack_.top().f;
      stack_.pop();
      std::size_t k = stack_.size();
      std::size_t j = std::floor(0.5*(1+std::sqrt(1+8*k)));
      std::size_t i = (j*j+j-2)/2 - k;
      if (basis_k_.size() == j && basis_k_.back() == i) {
        basis_k_.pop_back();
        if (f)
          curr_dim_upper_++;
        else
          curr_dim_lower_--;
      }
    }
    
  public:
    Face_iterator(const Alcove_id& a_id,
                  const Simple_coxeter_system& scs,
                  std::size_t value_dim_lower,
                  std::size_t value_dim_upper)
      : a_id_(a_id),
        value_(a_id.level()),
        family_(scs.family()),
        is_end_(false),
        ambient_dimension_(scs.dimension_),
        curr_dim_lower_(0),
        curr_dim_upper_(scs.dimension_),
        value_dim_lower_(value_dim_lower),
        value_dim_upper_(value_dim_upper)
    {
      if (!a_id.is_fixed(0)) {
        if (!stack_push(false, false)) {
          stack_pop();
          if (!stack_push(true, false)) {
            stack_pop();
            elementary_increment();
          }
        }
      }
      else if (!stack_push(true, false)) {
        stack_pop();
        elementary_increment();
      }
      update_value();
    }

    Face_iterator() : is_end_(true) {}
    
    
  protected:
    struct State_ {
      bool f, // is it fixed?
           p; // has the value different from the original?
      bool operator==(const State_& rhs) const {
        return f == rhs.f && p == rhs.p;
      }
    };
    Alcove_id a_id_; 
    Alcove_id value_;
    char family_;
    bool is_end_;
    std::size_t ambient_dimension_;
    std::size_t curr_dim_lower_, curr_dim_upper_;
    std::size_t value_dim_lower_, value_dim_upper_;
    std::vector<unsigned> basis_k_;
    std::stack<State_> stack_;
  };

  
  typedef boost::iterator_range<Face_iterator> Face_range;
  Face_range face_range(const Alcove_id& a_id,
                        std::size_t value_dim_lower,
                        std::size_t value_dim_upper) const {
    return Face_range(Face_iterator(a_id, *this, value_dim_lower, value_dim_upper),
                      Face_iterator());
  }

  Face_range face_range(const Alcove_id& a_id,
                        std::size_t value_dim) const {
    return face_range(a_id, value_dim, value_dim);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Coface range
  //////////////////////////////////////////////////////////////////////////////////////////////////

  class Coface_iterator : public boost::iterator_facade< Coface_iterator,
                                                         Alcove_id const,
                                                         boost::forward_traversal_tag> {
  protected:
    // using State_pair = std::pair<std::size_t, bool>;
    friend class boost::iterator_core_access;

    bool triplet_check(std::size_t i, std::size_t l, std::size_t j) const {
      int k  = (j*j+j-2)/2 - i;
      int k1 = (l*l+l-2)/2 - i;
      int k2 = (j*j+j-2)/2 - l;
      if (value_.is_fixed(k))
        return (value_.is_fixed(k1) && value_.is_fixed(k2) && value_[k1] + value_[k2] == value_[k]) ||
               (!value_.is_fixed(k1) && !value_.is_fixed(k2) && value_[k1] + value_[k2] + 1 == value_[k]);
      else
        return ((value_.is_fixed(k1) xor value_.is_fixed(k2)) && value_[k] == value_[k1] + value_[k2]) ||
               (!value_.is_fixed(k1) && !value_.is_fixed(k2) && (value_[k] == value_[k1] + value_[k2] ||
                                                                 value_[k] == value_[k1] + value_[k2] + 1));
    }
    
    void update_value() {
      if (is_end_)
        return;
      switch (family_) {
      case 'A': {
        // current position variables
        while (!is_end_) {
          std::size_t k = value_.size();
          std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
          std::size_t i = (j*j + j - 2)/2 - k;
          if (stack_.top().f)
            if (stack_.top().p)
              value_.push_back(a_id_[k]+1, true);
            else
              value_.push_back(a_id_[k], true);
          else
            if (stack_.top().p)
              value_.push_back(a_id_[k]-1, false);
            else
              value_.push_back(a_id_[k], false);
          bool curr_state_is_valid = true;
          if (basis_k_.size() == j && i != basis_k_.back())
            for (std::size_t l = basis_k_.back(); l < j && curr_state_is_valid; ++l)
              curr_state_is_valid = triplet_check(i,l,j);
          else
            for (std::size_t l = i + 1; l < j && curr_state_is_valid; ++l)
              curr_state_is_valid = triplet_check(i,l,j);
          if (curr_state_is_valid &&
              value_.size() == a_id_.size()) {
            value_.set_dimension(curr_dim_lower_);
            return;
          }
          if (!curr_state_is_valid) {
            elementary_increment();
            continue;
          }
          if (a_id_.is_fixed(k+1)) {
            if (!stack_push(false, false)) {              
              stack_pop();
              if (!stack_push(true, false)) {
                stack_pop();
                elementary_increment();
              }
            }
          }
          else 
            if (!stack_push(false, false)) {
              stack_pop();
              elementary_increment();
            }
        }
        break;
      }
      default:
        std::cerr << "Simple_coxeter_system::coface_iterator : The family " << family_ << " is not supported. "
                  << "Please use A family for the constructor (in capital).\n";
      }
    }
    
    bool equal(Coface_iterator const& other) const {
      return (is_end_ && other.is_end_) || (stack_ == other.stack_);
    }
    Alcove_id const& dereference() const {
      return value_;
    }
    void elementary_increment() {
      if (is_end_)
        return;
      while (!stack_.empty()) {
        bool f = stack_.top().f;
        bool p = stack_.top().p;
        stack_pop();
        value_.pop_back();
        if (a_id_.is_fixed(stack_.size()))
          if (f)
            continue;
          else
            if (p)
              if (stack_push(true, false)) {
                break;
              }
              else {
                stack_pop();
                continue;
              }
            else
              if (stack_push(false, true)) {
                break;
              }
              else {
                stack_pop();
                continue;
              }
        else
          continue;
      }
      if (stack_.empty())
        is_end_ = true;
    }
    
    void increment() {
      elementary_increment();
      update_value();
    }

    bool stack_push(bool fixed, bool plus_one) {
      std::size_t k = stack_.size();
      std::size_t j = std::floor(0.5*(1+std::sqrt(1+8*k)));
      std::size_t i = (j*j+j-2)/2 - k;
      if (basis_k_.size() == j - 1) {
        bool is_basis_vect = true;
        for (std::size_t l = i + 1; l < j && is_basis_vect; ++l) {
          std::size_t k1 = (l*l+l-2)/2 - i;
          is_basis_vect = !value_.is_fixed(k1);
        }
        if (!fixed)
          for (std::size_t l = 0; l < i && is_basis_vect; ++l) {
           is_basis_vect = false;
            for (std::size_t m = l+1; m < j && !is_basis_vect; ++m) {
              std::size_t k1 = (m*m+m-2)/2 - l;
              is_basis_vect = value_.is_fixed(k1);
            }
          }
        if (is_basis_vect) {
          basis_k_.push_back(i);
          if (fixed)
            curr_dim_upper_--;
          else
            curr_dim_lower_++;
        }
      }
      stack_.push({fixed, plus_one});        
      return curr_dim_lower_ <= curr_dim_upper_ &&
             curr_dim_lower_ <= value_dim_upper_ &&
             curr_dim_upper_ >= value_dim_lower_;
    }

    void stack_pop() {
      bool f = stack_.top().f;
      stack_.pop();
      std::size_t k = stack_.size();
      std::size_t j = std::floor(0.5*(1+std::sqrt(1+8*k)));
      std::size_t i = (j*j+j-2)/2 - k;
      if (basis_k_.size() == j && basis_k_.back() == i) {
        basis_k_.pop_back();
        if (f)
          curr_dim_upper_++;
        else
          curr_dim_lower_--;
      }
    }
    
  public:
    Coface_iterator(const Alcove_id& a_id,
                  const Simple_coxeter_system& scs,
                  std::size_t value_dim_lower,
                  std::size_t value_dim_upper)
      : a_id_(a_id),
        value_(a_id.level()),
        family_(scs.family()),
        is_end_(false),
        ambient_dimension_(scs.dimension_),
        curr_dim_lower_(0),
        curr_dim_upper_(scs.dimension_),
        value_dim_lower_(value_dim_lower),
        value_dim_upper_(value_dim_upper)
    {
      if (a_id.is_fixed(0)) {
        if (!stack_push(false, false)) {
          stack_pop();
          if (!stack_push(true, false)) {
            stack_pop();
            elementary_increment();
          }
        }
      }
      else if (!stack_push(false, false)) {
        stack_pop();
        elementary_increment();
      }
      update_value();
    }

    Coface_iterator() : is_end_(true) {}
    
    
  protected:
    struct State_ {
      bool f, // is it fixed?
           p; // has the value different from the original?
      bool operator==(const State_& rhs) const {
        return f == rhs.f && p == rhs.p;
      }
    };
    Alcove_id a_id_; 
    Alcove_id value_;
    char family_;
    bool is_end_;
    std::size_t ambient_dimension_;
    std::size_t curr_dim_lower_, curr_dim_upper_;
    std::size_t value_dim_lower_, value_dim_upper_;
    std::vector<unsigned> basis_k_;
    std::stack<State_> stack_;
  };

  
  typedef boost::iterator_range<Coface_iterator> Coface_range;
  Coface_range coface_range(const Alcove_id& a_id,
                            std::size_t value_dim_lower,
                            std::size_t value_dim_upper) const {
    return Coface_range(Coface_iterator(a_id, *this, value_dim_lower, value_dim_upper),
                        Coface_iterator());
  }

  Coface_range coface_range(const Alcove_id& a_id,
                            std::size_t value_dim) const {
    return coface_range(a_id, value_dim, value_dim);
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Vertex_range
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
public:  
  
  /** Add the vertices of the given simplex to a vertex-simplex map.
   * The size of si_it->first is d*(d+1)/2.
   */
  std::vector<Alcove_id> vertices_of_simplex(const Alcove_id& a_id) const
  {
    std::vector<Alcove_id> vertices;
    for (auto f_id: face_range(a_id, 0, 0)) {
      Alcove_id v_id(f_id.level(), 0);
      for (std::size_t i = 0; i < dimension_; ++i)
        v_id.push_back(f_id[i*(i+1)/2]);
      vertices.push_back(v_id);
    }
    return vertices;
  }

  /** Check if the given simplex and vertex are adjacent.
   */ 
  bool is_adjacent(const Vertex_id& v_id, const Alcove_id& a_id) const {
    auto alcove_it = a_id.begin();
    for (unsigned i = 1; i < v_id.size(); ++i) {
      int sum = 0; 
      for (unsigned j = i; j >= 1; --j) {
        sum += v_id[j-1];
        double
          v_plane = ((double)sum)/v_id.level()/vertex_level_,
          a_plane = ((double)*alcove_it++)/a_id.level()/vertex_level_;
        if (v_plane < a_plane || v_plane > a_plane + 1)
          return false;
      }
    }
    return true;
  }

  std::vector<double> barycenter(const Alcove_id& a_id) const {
    short d = dimension_;
    std::vector<double> result(d, 0);
    std::vector<Vertex_id> vertices = vertices_of_simplex(a_id);
    // std::cout << "Vertices of " << a_id << " are " << vertices << ". ";
    FT denom = a_id.level()*vertex_level_;
    // Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> chol(root_t_);
    Eigen::PartialPivLU<Matrix> chol(root_t_);
    for (auto v: vertices) {
      Eigen::VectorXd b(d);
      for (int i = 0; i < d; i++) {
        b(i) = v[i]/denom;
      }
      Eigen::VectorXd x = chol.solve(b);
      // std::cout << "Vertex " << v << "\n";
      // std::cout << x << "\n";
      for (int i = 0; i < d; i++)
        result[i] += x(i)/(d+1);
    }
    // std::cout << " The barycenter is " << result << ".\n";
    return result;
  }
  
  template <class VMap,
            class Simplex_range>
  void write_mesh(VMap& v_map, Simplex_range& range, std::string file_name = "toplex.mesh") const
  {
    short d = dimension_;
  
    std::ofstream ofs (file_name, std::ofstream::out);
    if (d <= 2)
      ofs << "MeshVersionFormatted 1\nDimension 2\n";
    else
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
    ofs << "Vertices\n" << v_map.size() << "\n";

    // std::vector<std::vector<double>> W;
    for (auto m: v_map) {
      FT denom = m.first.level()*vertex_level_;
      Eigen::VectorXd b(d);
      for (int i = 0; i < d; i++) {
        b(i) = m.first[i]/denom;
      }
      // Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> chol(root_t_);
      Eigen::PartialPivLU<Matrix> chol(root_t_);
      Eigen::VectorXd x = chol.solve(b);
      // if(chol.info()!=Eigen::Success) {
      //   std::cout << "solving failed\n";
      // }
      // std::cout << x << "\n\n";
      for (int i = 0; i < d; i++)
        ofs << x(i) << " ";
      ofs << "1" << std::endl;
    }

    struct Pointer_compare {
      typedef typename VMap::iterator Pointer;
      bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
        return lhs->first < rhs->first;
      }
    };    
    
    if (d == 2) {
      std::vector<typename Simplex_range::value_type> edges, triangles;
      for (auto s: range)
        if (s.size() == 2)
          edges.push_back(s);
        else if (s.size() == 3)
          triangles.push_back(s);
      ofs << "Edges " << edges.size() << "\n";
      for (auto s: edges) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "515" << std::endl;
      }
      
      ofs << "Triangles " << triangles.size() << "\n";
      for (auto s: triangles) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "516" << std::endl;
      }
    }
    else {
      std::vector<typename Simplex_range::value_type> edges, triangles, tetrahedra;
      // std::vector<std::vector<int> > facets;
      for (auto s: range)
        if (s.size() == 2)
          edges.push_back(s);
        else if (s.size() == 3)
          triangles.push_back(s);
        else if (s.size() == 4)
          tetrahedra.push_back(s);
      // for (auto t: tetrahedra) {
      //   for (unsigned i = 0; i < t.size(); i++) {
      //     std::vector<int> facet;
      //     auto it = t.begin();
      //     for (unsigned j = 0; it != t.end(); ++j, ++it)
      //       if (i != j)
      //         facet.push_back(*it);
      //     facets.push_back(facet);
      //   }
      // }
      // ofs << "Edges " << edges.size() << "\n";
      // for (auto s: edges) {
      //   for (auto v: s) {
      //     ofs << v+1 << " ";
      //   }
      //   ofs << "514" << std::endl;
      // }
      
      // ofs << "Triangles " << triangles.size() + facets.size() << "\n";
      // for (auto s: triangles) {
      //   for (auto v: s) {
      //     ofs << v+1 << " ";
      //   }
      //   ofs << "515" << std::endl;
      // }
      ofs << "Triangles " << triangles.size() << "\n";
      for (auto s: triangles) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "515" << std::endl;
      }
      // for (auto s: facets) {
      //   for (auto v: s) {
      //     ofs << v+1 << " ";
      //   }
      //   ofs << "320" << std::endl;
      // }

      ofs << "Tetrahedra " << tetrahedra.size() << "\n";
      for (auto s: tetrahedra) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "516" << std::endl;
      }
    }
  }
  
};



// Print the Simple_coxeter_system in os.
std::ostream& operator<<(std::ostream & os, const Simple_coxeter_system & scs) {
  os << scs.family() << scs.dimension();
  return os;
}

#endif
