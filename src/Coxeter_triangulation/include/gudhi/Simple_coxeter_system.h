#ifndef SIMPLE_COXETER_SYSTEM_H_
#define SIMPLE_COXETER_SYSTEM_H_

#include <iostream>
#include <vector>
#include <utility>
#include <exception>
#include <Eigen/Sparse>
// #include <Eigen/SPQRSupport>

class Simple_coxeter_system {

  typedef double FT;
  typedef Eigen::SparseMatrix<FT> Matrix;
  typedef Eigen::Triplet<FT> Triplet;
  typedef std::vector<int> Alcove_id;
  typedef Alcove_id Vertex_id;
  typedef std::vector<std::pair<int,int>> Change_range;
  typedef std::vector<Change_range> Local_changes;
  
  class wrong_family : public std::exception {  
  } wrong_family_exception_;

public:
  
  Matrix root_t_;
  char family_;
  unsigned short dimension_;
  unsigned short vertex_level_ = 1;
  
  Simple_coxeter_system() {
  }
  
  Simple_coxeter_system(char family, unsigned short dimension)
    : root_t_(dimension, dimension), family_(family), dimension_(dimension){
    switch (family) {
    case 'A': {
      unsigned short d = dimension;
      vertex_level_ = 1;
      std::vector<Triplet> cartan_triplets;
      cartan_triplets.reserve(3*d-2);
      for (unsigned i = 0; i < d; i++) {
        cartan_triplets.push_back(Triplet(i,i,2.0));
      }
      for (unsigned i = 1; i < d; i++) {
        cartan_triplets.push_back(Triplet(i-1,i,-1.0));
        cartan_triplets.push_back(Triplet(i,i-1,-1.0));
      }
      Matrix cartan(d,d);
      cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      // std::cout << "cartan =" << std::endl << cartan << std::endl;
      Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      root_t_ = chol.matrixL();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }
    case 'B': {
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 2);
      std::vector<Triplet> cartan_triplets;
      cartan_triplets.reserve(3*d-2);
      for (int i = 0; i < d-1; i++) {
        cartan_triplets.push_back(Triplet(i,i,4.0));
      }
      cartan_triplets.push_back(Triplet(d-1,d-1,2.0));
      for (int i = 1; i < d; i++) {
        cartan_triplets.push_back(Triplet(i-1,i,-2.0));
        cartan_triplets.push_back(Triplet(i,i-1,-2.0));
      }
      Matrix cartan(d,d);
      cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      // std::cout << "cartan =" << std::endl << cartan << std::endl;
      Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      root_t_ = chol.matrixL();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'C': {
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 2);
      std::vector<Triplet> cartan_triplets;
      cartan_triplets.reserve(3*d-2);
      for (int i = 0; i < d-1; i++) {
        cartan_triplets.push_back(Triplet(i,i,2.0));
      }
      cartan_triplets.push_back(Triplet(d-1,d-1,4.0));
      for (int i = 1; i < d-1; i++) {
        cartan_triplets.push_back(Triplet(i-1,i,-1.0));
        cartan_triplets.push_back(Triplet(i,i-1,-1.0));
      }
      cartan_triplets.push_back(Triplet(d-1,d-2,-2.0));
      cartan_triplets.push_back(Triplet(d-2,d-1,-2.0));
      Matrix cartan(d,d);
      cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      // std::cout << "cartan =" << std::endl << cartan << std::endl;
      Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      root_t_ = chol.matrixL();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'D': {   
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 3);
      std::vector<Triplet> cartan_triplets;
      cartan_triplets.reserve(3*d-2);
      for (unsigned i = 0; i < d; i++) {
        cartan_triplets.push_back(Triplet(i,i,2.0));
      }
      for (int i = 1; i < d-1; i++) {
        cartan_triplets.push_back(Triplet(i-1,i,-1.0));
        cartan_triplets.push_back(Triplet(i,i-1,-1.0));
      }
      cartan_triplets.push_back(Triplet(d-1,d-3,-1.0));
      cartan_triplets.push_back(Triplet(d-3,d-1,-1.0));
      Matrix cartan(d,d);
      cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      // std::cout << "cartan =" << std::endl << cartan << std::endl;
      Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      root_t_ = chol.matrixL();
      // std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      break;
    }  
    case 'E': {   
      unsigned short d = dimension;
      vertex_level_ = 2;
      assert(d >= 6 && d <= 8);
      std::vector<Triplet> cartan_triplets;
      cartan_triplets.reserve(3*d-2);
      for (unsigned i = 0; i < d; i++) {
        cartan_triplets.push_back(Triplet(i,i,2.0));
      }
      for (int i = 3; i < d; i++) {
        cartan_triplets.push_back(Triplet(i-1,i,-1.0));
        cartan_triplets.push_back(Triplet(i,i-1,-1.0));
      }
      cartan_triplets.push_back(Triplet(2,0,-1.0));
      cartan_triplets.push_back(Triplet(3,1,-1.0));
      cartan_triplets.push_back(Triplet(0,2,-1.0));
      cartan_triplets.push_back(Triplet(1,3,-1.0));
      Matrix cartan(d,d);
      cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
      std::cout << "cartan =" << std::endl << cartan << std::endl;
      Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
      root_t_ = chol.matrixL();
      std::cout << "root^t =" << std::endl << root_t_ << std::endl;
      Eigen::MatrixXf base(8,8);      
      std::vector<Triplet> base_triplets;
      base(0,0) = 0.5;
      base(7,0) = 0.5;
      base(1,0) = -0.5;
      base(2,0) = -0.5;
      base(3,0) = -0.5;
      base(4,0) = -0.5;
      base(5,0) = -0.5;
      base(6,0) = -0.5;
      base(0,1) = 1.0;
      base(1,1) = 1.0;
      for (int i = 0; i < d-2; ++i) {
        base(i, i+2) =  -1.0;
        base(i+1, i+2) =  1.0;
      }
      if (d <= 6) {
        base(5,6) = 1.0;
        base(6,6) = -1.0;
      }
      if (d <= 7) {
        base(6,7) = 1.0;
        base(7,7) = 1.0;
      }
      std::cout << "base = " << std::endl << base << std::endl;
      std::cout << "base^{-1} = " << std::endl << base.inverse() << std::endl;
      // for (int i = )
      break;
    }  
    default :
      std::cerr << "Simple_coxeter_system : The family " << family << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
  }

  unsigned short dimension() const {
    return dimension_;
  }

  char family() const {
    return family_;
  }

  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point,
            class OutputIterator>
  void alcove_coordinates(const Point& p, int level, OutputIterator output_it) const
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

  unsigned pos_root_count() {
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
  
private:  

  int gcd(int a, int b) const {
    return b == 0 ? a : gcd(b, a % b);
  }

  /** Common gcd simplification */
  template <class Id>
  Id reduced_id(Id& id) const {
    int common_gcd = 0;
    for (auto i: id) {
      common_gcd = gcd(i, common_gcd);
      if (common_gcd == 1)
        return id;
    }
    Id id_red(id);
    for (auto i_it = id_red.begin(); i_it != id_red.end(); ++i_it) {
      *i_it = *i_it / common_gcd;
    }
    return id_red;
  }

  template <class S_id_iterator>
  bool valid_coordinate(const Vertex_id& v_id, S_id_iterator& s_it, unsigned& integers, std::vector<Triplet>& triplets) {
    int k = v_id.size();
    unsigned short d = dimension_; 
    switch (family_) {
    case 'A': {
      // e_i - e_j
      int sum = 0;
      for (unsigned i = k-1; i >= 1; i--) {
        sum += v_id[i];
        if (sum < *s_it || sum > *s_it + 1)
          return false;
        s_it++;
      }
      return true;
    }
    case 'B': {
      if (k == d+1) {
        // e_i
        int sum = 0;
        for (unsigned i = d; i >= 1; i--) {
          sum += v_id[i];
          if (sum < 2*(*s_it) || sum > 2*(*s_it) + 2)
            return false;
          if (sum % 2 == 0) {
            triplets.push_back(Triplet(integers, i-1, 1.0));
            integers++;
          }
          s_it++;
        }
        // e_i + e_j
        int glob_sum = 0;
        for (unsigned i = d; i >= 1; i--) {
          glob_sum += 2*v_id[i];
          int sum = glob_sum;
          for (short j = i-1; j >= 1; j--) {
            sum += v_id[j];
            if (sum < 2*(*s_it) || sum > 2*(*s_it) + 2)
              return false;
            if (sum % 2 == 0) {
              triplets.push_back(Triplet(integers, i-1, 1.0));
              triplets.push_back(Triplet(integers, j-1, 1.0));
              integers++;
            }
            s_it++;
          }
        }
        return true;
      }
      // e_i - e_j
      int sum = 0;
      for (unsigned i = k-1; i >= 1; i--) {
        sum += v_id[i];
        if (sum < 2*(*s_it) || sum > 2*(*s_it + 1))
          return false;
        if (sum % 2 == 0) {
          triplets.push_back(Triplet(integers, i-1, 1.0));
          triplets.push_back(Triplet(integers, k-1, -1.0));
          integers++;
        }
        s_it++;
      }
      return true;
    }
    case 'C': {
      if (k == d+1) {
        // 2*e_i
        int sum = -v_id[d];
        for (unsigned i = d; i >= 1; i--) {
          sum += 2*v_id[i];
          if (sum < 2*(*s_it) || sum > 2*(*s_it) + 2)
            return false;
          if (sum % 2 == 0) {
            triplets.push_back(Triplet(integers, i-1, 2.0));
            integers++;
          }
          s_it++;
        }
        // e_i + e_j
        int glob_sum = -v_id[d];
        for (unsigned i = d; i >= 1; i--) {
          glob_sum += 2*v_id[i];
          int sum = glob_sum;
          for (short j = i-1; j >= 1; j--) {
            sum += v_id[j];
            if (sum < 2*(*s_it) || sum > 2*(*s_it) + 2)
              return false;
            if (sum % 2 == 0) {
              triplets.push_back(Triplet(integers, i-1, 1.0));
              triplets.push_back(Triplet(integers, j-1, 1.0));
              integers++;
            }
            s_it++;
          }
        }
        return true;
      }
      // e_i - e_j
      int sum = 0;
      for (unsigned i = k-1; i >= 1; i--) {
        sum += v_id[i];
        if (sum < 2*(*s_it) || sum > 2*(*s_it + 1))
          return false;
        if (sum % 2 == 0) {
          triplets.push_back(Triplet(integers, i-1, 1.0));
          triplets.push_back(Triplet(integers, k-1, -1.0));
          integers++;
        }
        s_it++;
      }
      return true;
    }
    case 'D': {
      if (k == d+1) {
        // e_i + e_j
        int glob_sum = -v_id[d]-v_id[d-1];
        for (int i = d; i >= 1; i--) {
          glob_sum += 2*v_id[i];
          int sum = glob_sum;
          for (short j = i-1; j >= 1; j--) {
            sum += v_id[j];
            if (sum < 2*(*s_it) || sum > 2*(*s_it) + 2)
              return false;
            if (sum % 2 == 0) {
              triplets.push_back(Triplet(integers, i-1, 1.0));
              triplets.push_back(Triplet(integers, j-1, 1.0));
              integers++;
            }
            s_it++;
          }
        }
        return true;
      }
      // e_i - e_j
      int sum = 0;
      for (unsigned i = k-1; i >= 1; i--) {
        sum += v_id[i];
        if (sum < 2*(*s_it) || sum > 2*(*s_it + 1))
          return false;
        if (sum % 2 == 0) {
          triplets.push_back(Triplet(integers, i-1, 1.0));
          triplets.push_back(Triplet(integers, k-1, -1.0));             
          integers++;
        }
        s_it++;
      }
      return true;
    }
    default :
      std::cerr << "Simple_coxeter_system::valid_coordinate : The family " << family_ << " is not supported. "
                << "Please use A, B, C or D family for the constructor (in capital).\n";
      throw wrong_family_exception_;
    }
    
  }
  
  /** Add the vertices of the given simplex to a vertex-simplex map.
   */
  template <class S_id_iterator>
  void rec_vertices_of_simplex(Vertex_id& v_id, S_id_iterator s_it, std::vector<Vertex_id>& vertices, unsigned& integers, std::vector<Triplet>& triplets)
  {
    unsigned short d = dimension_;
    int k = v_id.size();
    if (family_ == 'A') {
      if (k == d+1) {
        vertices.emplace_back(v_id);
        return;
      }
      for (unsigned i = 0; i <= vertex_level_; i++) {
        v_id.push_back(vertex_level_*(*s_it) + (int)i);
        S_id_iterator s_it_copy(s_it);
        if (valid_coordinate(v_id, s_it_copy, integers, triplets))
          rec_vertices_of_simplex(v_id, s_it_copy, vertices, integers, triplets);
        v_id.pop_back();
      }
    }
    else {
      if (k == d+1) {
        if (integers < d)
          return;
        Matrix int_roots(integers, d);
        int_roots.setFromTriplets(triplets.begin(), triplets.end());
        // Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> spQR(int_roots);
        Eigen::SparseQR<Matrix, Eigen::NaturalOrdering<int>> spQR(int_roots);
        // Eigen::SPQR<Matrix> spQR(int_roots);
        if (spQR.rank() == d)
          vertices.emplace_back(v_id);
        return;
      }
      for (unsigned i = 0; i <= vertex_level_; i++) {
        v_id.push_back(vertex_level_*(*s_it) + (int)i);
        S_id_iterator s_it_copy(s_it);
        unsigned integers_copy = integers;
        std::vector<Triplet> triplets_copy(triplets);
        if (valid_coordinate(v_id, s_it_copy, integers_copy, triplets_copy))
          rec_vertices_of_simplex(v_id, s_it_copy, vertices, integers_copy, triplets_copy);
        v_id.pop_back();
      }
    }
  }

public:  
  
  /** Add the vertices of the given simplex to a vertex-simplex map.
   * The size of si_it->first is d*(d+1)/2.
   */
  std::vector<Vertex_id> vertices_of_simplex(Alcove_id ai_id)
  {
    unsigned d = dimension_;
    Vertex_id v_id(1,*ai_id.begin() * vertex_level_);
    v_id.reserve(d+1);
    std::vector<Vertex_id> vertices;
    vertices.reserve(d+1);
    unsigned integers = 0;
    std::vector<Triplet> triplets;
    rec_vertices_of_simplex(v_id, ai_id.begin()+1, vertices, integers, triplets);
    return vertices;
  }

  /** Check if the given simplex and vertex are adjacent.
   */ 
  bool is_adjacent(const Vertex_id& v_id, const Alcove_id& a_id) const {
    auto alcove_it = a_id.begin()+1;
    for (unsigned i = 1; i < v_id.size(); ++i) {
      int sum = 0; 
      for (unsigned j = i; j >= 1; --j) {
        sum += v_id[j];
        double
          v_plane = ((double)sum)/v_id[0]/vertex_level_,
          a_plane = ((double)*alcove_it++)/a_id[0]/vertex_level_;
        if (v_plane < a_plane || v_plane > a_plane + 1)
          return false;
      }
    }
    return true;
  }

  template <class VMap,
            class AMap>
  void write_coxeter_mesh(VMap& v_map, AMap& a_map, std::string file_name = "coxeter.mesh") const
  {
    short d = dimension_;
    if (d > 3);
  
    std::ofstream ofs (file_name, std::ofstream::out);
    if (d <= 2)
      ofs << "MeshVersionFormatted 1\nDimension 2\n";
    else
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
    ofs << "Vertices\n" << v_map.size() << "\n";

    // std::vector<std::vector<double>> W;
    for (auto m: v_map) {
      FT denom = m.first[0]*vertex_level_;
      Eigen::VectorXd b(d);
      for (int i = 0; i < d; i++) {
        b(i) = m.first[i+1]/denom;
      }
      // std::cout << b << "\n\n";
      // Eigen::SimplicialLDLT<Matrix, Eigen::Lower> chol(root_t_);
      Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> chol(root_t_);
      Eigen::VectorXd x = chol.solve(b);
      if(chol.info()!=Eigen::Success) {
        std::cout << "solving failed\n";
      }
      // std::cout << x << "\n\n";
      for (int i = 0; i < d; i++)
        ofs << x(i) << " ";
      // std::vector<double> np;
      // for (int i = 0; i < d; i++)
      //   np.push_back(x(i));
      // W.push_back(np);
      ofs << "216" << std::endl;
    }
    // std::map<int,std::vector<int>> sv_map;
    // int j = 1;
    // for (auto m: v_map) {
    //   for (auto s: m.second) {
    //     auto find_it = sv_map.find(s);
    //     if (find_it == sv_map.end())
    //       sv_map.emplace(s, std::vector<int>(1,j));
    //     else
    //       find_it->second.push_back(j);
    //   }
    //   j++;
    // }

    struct Pointer_compare {
      typedef typename VMap::iterator Pointer;
      bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
        return lhs->first < rhs->first;
      }
    };    
    std::map<typename VMap::iterator, int, Pointer_compare> vi_map;
    unsigned index = 1;
    for (auto v_it = v_map.begin(); v_it != v_map.end(); ++v_it, ++index)
      vi_map.emplace(v_it, index);

    if (d == 2)
      ofs << "Triangles " << a_map.size() << "\n";
    else
      ofs << "Tetrahedra " << a_map.size() << "\n";
    
    for (auto a: a_map) {
      for (auto m_it: std::get<2>(a.second)) {
        ofs << vi_map[m_it] << " ";
      }
      ofs << "516" << std::endl;
    }
  }

  template <class VMap,
            class Fake_simplex_tree>
  void write_toplex_mesh(VMap& v_map, Fake_simplex_tree& stree, std::string file_name = "toplex.mesh") const
  {
    short d = dimension_;
    if (d > 3);
  
    std::ofstream ofs (file_name, std::ofstream::out);
    if (d <= 2)
      ofs << "MeshVersionFormatted 1\nDimension 2\n";
    else
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
    ofs << "Vertices\n" << v_map.size() << "\n";

    // std::vector<std::vector<double>> W;
    for (auto m: v_map) {
      FT denom = m.first[0]*vertex_level_;
      Eigen::VectorXd b(d);
      for (int i = 0; i < d; i++) {
        b(i) = m.first[i+1]/denom;
      }
      // std::cout << b << "\n\n";
      // Eigen::SimplicialLDLT<Matrix, Eigen::Lower> chol(root_t_);
      Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>> chol(root_t_);
      Eigen::VectorXd x = chol.solve(b);
      if(chol.info()!=Eigen::Success) {
        std::cout << "solving failed\n";
      }
      // std::cout << x << "\n\n";
      for (int i = 0; i < d; i++)
        ofs << x(i) << " ";
      // std::vector<double> np;
      // for (int i = 0; i < d; i++)
      //   np.push_back(x(i));
      // W.push_back(np);
      ofs << "1" << std::endl;
    }
    // std::map<int,std::vector<int>> sv_map;
    // int j = 1;
    // for (auto m: v_map) {
    //   for (auto s: m.second) {
    //     auto find_it = sv_map.find(s);
    //     if (find_it == sv_map.end())
    //       sv_map.emplace(s, std::vector<int>(1,j));
    //     else
    //       find_it->second.push_back(j);
    //   }
    //   j++;
    // }

    struct Pointer_compare {
      typedef typename VMap::iterator Pointer;
      bool operator()(const Pointer& lhs, const Pointer& rhs) const { 
        return lhs->first < rhs->first;
      }
    };    
    
    if (d == 2) {
      std::vector<typename Fake_simplex_tree::Simplex> edges, triangles;
      for (auto s: stree.max_simplices())
        if (s.size() == 2)
          edges.push_back(s);
        else
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
      std::vector<typename Fake_simplex_tree::Simplex> edges, triangles, tetrahedra;
      std::vector<std::vector<int>> facets;
      for (auto s: stree.max_simplices())
        if (s.size() == 2)
          edges.push_back(s);
        else if (s.size() == 3)
          triangles.push_back(s);
        else
          tetrahedra.push_back(s);
      for (auto t: tetrahedra) {
        for (unsigned i = 0; i < t.size(); i++) {
          std::vector<int> facet;
          auto it = t.begin();
          for (unsigned j = 0; it != t.end(); ++j, ++it)
            if (i != j)
              facet.push_back(*it);
          facets.push_back(facet);
        }
      }
      // ofs << "Edges " << edges.size() << "\n";
      // for (auto s: edges) {
      //   for (auto v: s) {
      //     ofs << v+1 << " ";
      //   }
      //   ofs << "514" << std::endl;
      // }
      
      ofs << "Triangles " << triangles.size() + facets.size() << "\n";
      for (auto s: triangles) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "515" << std::endl;
      }
      for (auto s: facets) {
        for (auto v: s) {
          ofs << v+1 << " ";
        }
        ofs << "320" << std::endl;
      }

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
std::ostream& operator<<(std::ostream & os, Simple_coxeter_system & scs) {
  os << scs.family() << scs.dimension();
  return os;
}

#endif
