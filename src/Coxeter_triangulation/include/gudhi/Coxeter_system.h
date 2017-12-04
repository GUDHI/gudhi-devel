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

  unsigned short dimension_;
  
public:
  
  std::vector<Simple_coxeter_system> simple_system_range_;
  
  Coxeter_system(char family, short dimension)
    : simple_system_range_(1, Simple_coxeter_system(family, dimension)), dimension_(dimension) {
  }

  Coxeter_system(const std::vector<Simple_coxeter_system>& simple_system_range)
    : simple_system_range_(simple_system_range) {
    dimension_ = 0;
    for (auto scs: simple_system_range)
      dimension_ += scs.dimension();
  }
  
  void emplace(const Simple_coxeter_system& rhs) {
    simple_system_range_.emplace_back(rhs);
    dimension_ += rhs.dimension();
  }

  void emplace(char family, short dimension) {
    simple_system_range_.emplace_back(Simple_coxeter_system(family, dimension));
    dimension_ += dimension;
  }
  
  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point>
  Alcove_id alcove_coordinates(Point& p, int level)
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

  std::vector<Vertex_id> vertices_of_simplex(Alcove_id ai_id)
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
      chunks.emplace_back(scs.vertices_of_simplex(ai_id));
    }
    std::vector<std::vector<Vertex_id>::iterator> iterators;
    for (auto chunk: chunks)
      iterators.emplace_back(chunk.begin());
    Vertex_id v_id(1, level);
    // while (iterator[0] != chunks[0].end()) {
    //   if (v_id.size() == dimension_+1)
    //     vertices.emplace_back(v_id);
    // }
    return vertices;
  }
  
};

#endif
