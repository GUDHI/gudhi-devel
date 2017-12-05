#ifndef COXETER_SYSTEM_H_
#define COXETER_SYSTEM_H_

#include <iostream>
#include <vector>
#include <exception>
#include <Eigen/Sparse>
#include <gudhi/Simple_coxeter_system.h>

class Coxeter_system  {

  typedef double FT;
  typedef Eigen::SparseMatrix<FT> Matrix;
  typedef Eigen::Triplet<FT> Triplet;
  typedef std::vector<int> Alcove_id;
  typedef Alcove_id Vertex_id;
  
  class wrong_family : public std::exception {  
  } wrong_family_exception_;

  std::vector<Simple_coxeter_system> simple_system_range_;
  unsigned short dimension_;
  
public:

  Coxeter_system()
    : dimension_(0) {
  }
  
  Coxeter_system(char family, short dimension)
    : simple_system_range_(1, Simple_coxeter_system(family, dimension)), dimension_(dimension) {
  }

  Coxeter_system(const std::vector<Simple_coxeter_system>& simple_system_range)
    : simple_system_range_(simple_system_range) {
    dimension_ = 0;
    for (auto scs: simple_system_range)
      dimension_ += scs.dimension();
  }
  
  void emplace_back(const Simple_coxeter_system& rhs) {
    simple_system_range_.emplace_back(rhs);
    dimension_ += rhs.dimension();
  }

  void emplace_back(char family, short dimension) {
    simple_system_range_.emplace_back(Simple_coxeter_system(family, dimension));
    dimension_ += dimension;
  }

  void pop_back() {
    dimension_ -= simple_system_range_.back().dimension();
    simple_system_range_.pop_back();
  }
  
  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point>
  Alcove_id alcove_coordinates(const Point& p, int level) const
  {
    Alcove_id a_id(1, level);
    auto p_it = p.begin();
    for (auto scs: simple_system_range_) {
      std::vector<FT> coordinate_segment;
      coordinate_segment.reserve(scs.dimension());
      for (unsigned i = 0; i < scs.dimension(); i++)
        coordinate_segment.push_back(*p_it++);
      scs.alcove_coordinates(coordinate_segment, level, std::back_inserter(a_id));
    }
    return a_id;
  }

private:

  int gcd(int a, int b) const {
    return b == 0 ? a : gcd(b, a % b);
  }

  /** Common gcd simplification */
  template <class Id>
  Id reduced_id(const Id& id) const {
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

  
  void rec_combine_chunks(std::vector<std::vector<Vertex_id>>::iterator chunks_it, std::vector<Vertex_id>& vertices, Vertex_id& v_id) const {
    int k = v_id.size()-1;
    if (k == dimension_) {
      vertices.push_back(reduced_id(v_id));
      return;
    }
    for (auto chunk: *chunks_it) {
      for (auto c_it = chunk.begin()+1; c_it != chunk.end(); ++c_it)
        v_id.push_back(*c_it);
      rec_combine_chunks(chunks_it+1, vertices, v_id);
      v_id.resize(v_id.size()-chunk.size()+1);
    }
  }

public:  
  
  std::vector<Vertex_id> vertices_of_alcove(const Alcove_id& ai_id) const
  {
    std::vector<Vertex_id> vertices;
    std::vector<std::vector<Vertex_id>> chunks;
    auto ai_it = ai_id.begin();
    int level = *ai_it++;
    for (auto scs: simple_system_range_) {
      std::vector<int> ai_id_part(1, level);
      unsigned pos_root_count = scs.pos_root_count();
      for (unsigned i = 0; i < pos_root_count; i++) {
        ai_id_part.push_back(*ai_it++);
      }
      chunks.emplace_back(scs.vertices_of_simplex(ai_id_part));
    }
    // std::vector<std::vector<Vertex_id>::iterator> iterators;
    // for (auto chunk: chunks)
    //   iterators.emplace_back(chunk.begin());
    Vertex_id v_id(1, level);
    rec_combine_chunks(chunks.begin(), vertices, v_id);
    return vertices;
  }
  
};

#endif
