#ifndef SIMPLE_COXETER_SYSTEM_H_
#define SIMPLE_COXETER_SYSTEM_H_

#include <iostream>
#include <vector>
#include <exception>
#include <Eigen/Sparse>

class Simple_coxeter_system {

  typedef double FT;
  typedef Eigen::SparseMatrix<FT> Matrix;
  typedef Eigen::Triplet<FT> Triplet;
  typedef std::vector<int> Alcove_id;
  typedef Alcove_id Vertex_id;
  
  class wrong_family : public std::exception {  
  } wrong_family_exception_;

public:
  
  Matrix root_t_;
  char family_;
  unsigned short dimension_;
  
  Simple_coxeter_system() {
  }
  
  Simple_coxeter_system(char family, unsigned short dimension)
    : root_t_(dimension, dimension), family_(family), dimension_(dimension){
    switch (family) {
    case 'A': {
      unsigned short d = dimension;
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
    default :
      std::cerr << "Simple_coxeter_system : The family " << family << " is not supported. "
                << "Please use A family for the constructor (in capital).\n";
      throw wrong_family_exception_;
      break;
    }
  }

  unsigned short dimension() const {
    return dimension_;
  }
  
  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point,
            class OutputIterator>
  void alcove_coordinates(Point& p, int level, OutputIterator output_it)
  {
    unsigned short d = p.size();
    assert(d == dimension_);
    Eigen::VectorXd p_vect(d);
    for (short i = 0; i < d; i++)
      p_vect(i) = p[i];
    Eigen::VectorXd scalprod_vect = root_t_ * p_vect;
    for (short i = 0; i < d; i++) {
      FT root_scalprod = 0;
      for (short j = i; j >= 0; j--) {
        root_scalprod += scalprod_vect(j);
       *output_it++ = std::floor(level * root_scalprod);
      }
    }
  }

  unsigned pos_root_count() {
    return dimension_*(dimension_ + 1)/2 ; 
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

  
  /** Add the vertices of the given simplex to a vertex-simplex map.
   * The size of si_it->first is d*(d+1)/2.
   */
  template <class S_id_iterator>
  void rec_vertices_of_simplex(Vertex_id& v_id, S_id_iterator s_it, unsigned d, std::vector<Vertex_id>& vertices)
  {
    unsigned k = v_id.size();
    if (k == d+1) {
      vertices.emplace_back(v_id);
      return;
    }
    int simplex_coord = *s_it; s_it++;
    v_id.push_back(simplex_coord);
    int sum = simplex_coord;
    S_id_iterator s_it_copy(s_it);
    bool valid = true;
    for (unsigned i = 1; i < k && valid; i++) {
      sum += v_id[k-i];
      if (sum < *s_it_copy || sum > *s_it_copy + 1)
        valid = false;
      s_it_copy++;
    }
    if (valid)
      rec_vertices_of_simplex(v_id, s_it_copy, d, vertices);
    v_id.pop_back();

    v_id.push_back(simplex_coord + 1);
    sum = simplex_coord + 1;
    s_it_copy = s_it;
    valid = true;
    for (unsigned i = 1; i < k && valid; i++) {
      sum += v_id[k-i];
      if (sum < *s_it_copy || sum > *s_it_copy + 1)
        valid = false;
      s_it_copy++;
    }
    if (valid)
      rec_vertices_of_simplex(v_id, s_it_copy, d, vertices);
    v_id.pop_back();
  }

public:  
  
  /** Add the vertices of the given simplex to a vertex-simplex map.
   * The size of si_it->first is d*(d+1)/2.
   */
  std::vector<Vertex_id> vertices_of_simplex(Alcove_id ai_id)
  {
    unsigned d = dimension_;
    Vertex_id v_id(1,*ai_id.begin());
    v_id.reserve(d+1);
    std::vector<Vertex_id> vertices;
    vertices.reserve(d+1);
    rec_vertices_of_simplex(v_id, ai_id.begin()+1, d, vertices);
    return vertices;
  }
  
};

#endif
