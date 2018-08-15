#ifndef COXETER_SYSTEM_H_
#define COXETER_SYSTEM_H_

#include <iostream>
#include <vector>
#include <exception>
#include <Eigen/Sparse>
#include <gudhi/Simple_coxeter_system_remastered.h>
// #include <gudhi/Simple_coxeter_system.h>
#include "../../example/cxx-prettyprint/prettyprint.hpp"

class Coxeter_system  {

  typedef double FT;
  typedef Eigen::SparseMatrix<FT> Matrix;
  typedef Eigen::Triplet<FT> Triplet;
  
  class wrong_family : public std::exception {  
  } wrong_family_exception_;

  std::vector<Simple_coxeter_system> simple_system_range_;
  unsigned short dimension_;

protected:
  
  
public:
  typedef typename Simple_coxeter_system::Alcove_id Alcove_id;
  typedef Alcove_id Vertex_id;
  typedef typename Simple_coxeter_system::Filtered_alcove Filtered_alcove;

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

  unsigned short dimension() const {
    return dimension_;
  }

  unsigned pos_root_count() const {
    unsigned value = 0;
    for (auto ss: simple_system_range_)
      value += ss.pos_root_count();
    return value;
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

  std::vector<Simple_coxeter_system>::const_iterator simple_coxeter_system_begin() const {
    return simple_system_range_.begin();
  }

  std::vector<Simple_coxeter_system>::const_iterator simple_coxeter_system_end() const {
    return simple_system_range_.end();
  }

  /** A conversion from Cartesian coordinates to the coordinates of the alcove containing the point.
   *  The matrix' rows are simple root vectors.
   */
  template <class Point>
  Alcove_id query_point_location(const Point& p, double level) const {
    Alcove_id a_id(level);
    auto p_it = p.begin();
    for (auto scs: simple_system_range_) {
      std::vector<FT> coordinate_segment;
      coordinate_segment.reserve(scs.dimension());
      for (unsigned i = 0; i < scs.dimension(); i++)
        coordinate_segment.push_back(*p_it++);
      scs.query_point_location(coordinate_segment, level, std::back_inserter(a_id));
    }
    return a_id;
  }

  template <class Point,
            class Visitor>
  void alcoves_of_ball(const Point& p,
                       double init_level,
                       double eps,
                       Visitor visitor,
                       bool root_coords = false) const
  {
    std::vector<std::vector<Filtered_alcove> > chunks;
    std::vector<std::vector<std::vector<Vertex_id> > > chunks_v;
    auto p_it = p.begin();
    for (auto scs: simple_system_range_) {
      std::vector<FT> p_part;
      unsigned dimension = scs.dimension();
      for (unsigned i = 0; i < dimension; i++) {
        p_part.push_back(*p_it++);
      }
      std::vector<Filtered_alcove> alcoves;
      std::vector<std::vector<Vertex_id> > vertices;
      scs.alcoves_of_ball(p_part, init_level, eps, alcoves, vertices, root_coords);
      chunks.emplace_back(alcoves);
      chunks_v.emplace_back(vertices);
    }
    Filtered_alcove a_id(Alcove_id(init_level, dimension_), 0);
    std::vector<std::vector<Vertex_id > > vertex_chunks;
    rec_combine_chunks_alcove(chunks.begin(),
                              chunks.end(),
                              chunks_v.begin(),
                              visitor,
                              eps,
                              a_id,
                              vertex_chunks);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Face range
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
private:

  class Face_iterator : public boost::iterator_facade< Face_iterator,
                                                        Alcove_id const,
                                                        boost::forward_traversal_tag> {
  protected:
    typedef typename std::vector<Simple_coxeter_system>::const_iterator Simple_coxeter_system_iterator;
    typedef typename Simple_coxeter_system::Face_iterator Scs_face_iterator;
    friend class boost::iterator_core_access;
    
    void update_value(std::size_t first_change_) {
      if (is_end_)
        return;
      for (std::size_t i = first_change_; i < decomposition_.size(); ++i) {
        Alcove_id face = *face_iterators_[i].first;
        for (std::size_t j = 0; j < face.size(); ++j)
          value_.push_back(face[j], face.is_fixed(j));
      }
      first_change_ = decomposition_.size();
    }

    std::size_t update_from(std::size_t pos, std::size_t rest) {
      for (std::size_t i = pos; i < decomposition_.size(); ++i) {
        decomposition_[i] =
          (rest < chunks_[i].dimension() ? rest : chunks_[i].dimension());
        auto face_range = scs_iterators_[i]->face_range(chunks_[i],
                                                        decomposition_[i]);
        face_iterators_[i] = std::make_pair(face_range.begin(), face_range.end());
        rest -= decomposition_[i];
      }
      return rest;
    }
    
    bool equal(Face_iterator const& other) const {
      return (is_end_ && other.is_end_) || (!is_end_ && !other.is_end_ && decomposition_ == other.decomposition_);
    }
    Alcove_id const& dereference() const {
      return value_;
    }
    void increment() {
      if (is_end_)
        return;
      std::size_t rest = 0;
      std::size_t pos = decomposition_.size() - 1;
      while (true) {
        value_.resize(value_.size() - chunks_[pos].size());
        face_iterators_[pos].first++;
        if (face_iterators_[pos].first == face_iterators_[pos].second) {
          if (decomposition_[pos] == 0) {
            if (pos == 0) {
              is_end_ = true;
              return;
            }
            rest += decomposition_[pos];
            pos--;
            continue;
          }
          decomposition_[pos]--;
          rest++;    
          auto face_range = scs_iterators_[pos]->face_range(chunks_[pos],
                                                            decomposition_[pos]);
          face_iterators_[pos] = std::make_pair(face_range.begin(), face_range.end());
        }
        if (update_from(pos + 1, rest)) {
          if (pos == 0) {
            is_end_ = true;
            return;
          }
          rest += decomposition_[pos];
          pos--;
          continue;
        }
        update_value(pos);
        return;
      }
    }
     
  public:
    Face_iterator(const Alcove_id& coface,
                   const Coxeter_system& cs,
                   std::size_t value_dimension)
      : value_(coface.level(), value_dimension),
        is_end_(false),
        decomposition_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin()),
        face_iterators_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin())
    {
      std::size_t pos = 0;
      for (auto scs_it = cs.simple_coxeter_system_begin();
           scs_it != cs.simple_coxeter_system_end();
           ++scs_it) {
        scs_iterators_.push_back(scs_it);
        Alcove_id chunk(coface.level(), scs_it->dimension());
        for (std::size_t i = pos; i < pos + scs_it->pos_root_count(); ++i)
          chunk.push_back(coface[i], coface.is_fixed(i));
        chunk.set_dimension(scs_it->alcove_dimension(chunk));
        chunks_.push_back(chunk);
        pos += scs_it->pos_root_count();
      }
      if (update_from(0, value_dimension)) {
        is_end_ = true;
        return;
      }
      update_value(0);
    }

  protected:
    Alcove_id value_;
    bool is_end_;
    std::vector<Alcove_id> chunks_;
    std::vector<std::size_t> decomposition_;
    std::vector<std::pair<Scs_face_iterator, Scs_face_iterator> > face_iterators_;
    std::vector<Simple_coxeter_system_iterator> scs_iterators_;
  };

  
public:
  typedef boost::iterator_range<Face_iterator> Face_range;
  
  Face_range face_range(const Alcove_id& a_id, std::size_t k) const {
    return Face_range(Face_iterator(a_id, *this, k),
                      Face_iterator(a_id, *this, dimension_ + 1));
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Coface range
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
private:

  class Coface_iterator : public boost::iterator_facade< Coface_iterator,
                                                        Alcove_id const,
                                                        boost::forward_traversal_tag> {
  protected:
    typedef typename std::vector<Simple_coxeter_system>::const_iterator Simple_coxeter_system_iterator;
    typedef typename Simple_coxeter_system::Coface_iterator Scs_coface_iterator;
    friend class boost::iterator_core_access;
    
    void update_value(std::size_t first_change_) {
      if (is_end_)
        return;
      for (std::size_t i = first_change_; i < decomposition_.size(); ++i) {
        Alcove_id coface = *coface_iterators_[i].first;
        for (std::size_t j = 0; j < coface.size(); ++j)
          value_.push_back(coface[j], coface.is_fixed(j));
      }
      first_change_ = decomposition_.size();
    }

    std::size_t update_from(std::size_t pos, std::size_t rest_codim) {
      for (std::size_t i = pos; i < decomposition_.size(); ++i) {
        unsigned chunk_codim = scs_iterators_[i]->dimension() - chunks_[i].dimension();
        decomposition_[i] = (rest_codim < chunk_codim ? rest_codim : chunk_codim);
        unsigned decomp_dim = scs_iterators_[i]->dimension() - decomposition_[i];
        auto coface_range = scs_iterators_[i]->coface_range(chunks_[i], decomp_dim);
        coface_iterators_[i] = std::make_pair(coface_range.begin(), coface_range.end());
        rest_codim -= decomposition_[i];
      }
      return rest_codim;
    }
    
    bool equal(Coface_iterator const& other) const {
      return (is_end_ && other.is_end_) ||
        (!is_end_ && !other.is_end_ && decomposition_ == other.decomposition_);
    }
    Alcove_id const& dereference() const {
      return value_;
    }
    void increment() {
      if (is_end_)
        return;
      std::size_t rest_codim = 0;
      std::size_t pos = decomposition_.size() - 1;
      while (true) {
        value_.resize(value_.size() - chunks_[pos].size());
        coface_iterators_[pos].first++;
        if (coface_iterators_[pos].first == coface_iterators_[pos].second) {
          if (decomposition_[pos] == 0) {
            if (pos == 0) {
              is_end_ = true;
              return;
            }
            rest_codim += decomposition_[pos];
            pos--;
            continue;
          }
          decomposition_[pos]--;
          rest_codim++;    
          unsigned decomp_dim = scs_iterators_[pos]->dimension() - decomposition_[pos];
          auto coface_range = scs_iterators_[pos]->coface_range(chunks_[pos], decomp_dim);
          coface_iterators_[pos] = std::make_pair(coface_range.begin(), coface_range.end());
        }
        if (update_from(pos + 1, rest_codim)) {
          if (pos == 0) {
            is_end_ = true;
            return;
          }
          rest_codim += decomposition_[pos];
          pos--;
          continue;
        }
        update_value(pos);
        return;
      }
    }
     
  public:
    Coface_iterator(const Alcove_id& cell,
                   const Coxeter_system& cs,
                   std::size_t value_dimension)
      : value_(cell.level(), value_dimension),
        is_end_(false),
        decomposition_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin()),
        coface_iterators_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin())
    {
      std::size_t pos = 0;
      for (auto scs_it = cs.simple_coxeter_system_begin();
           scs_it != cs.simple_coxeter_system_end();
           ++scs_it) {
        scs_iterators_.push_back(scs_it);
        Alcove_id chunk(cell.level(), scs_it->dimension());
        for (std::size_t i = pos; i < pos + scs_it->pos_root_count(); ++i)
          chunk.push_back(cell[i], cell.is_fixed(i));
        chunk.set_dimension(scs_it->alcove_dimension(chunk));
        chunks_.push_back(chunk);
        pos += scs_it->pos_root_count();
      }
      if (cs.dimension() < value_dimension || update_from(0, cs.dimension() - value_dimension)) {
        is_end_ = true;
        return;
      }
      update_value(0);
    }

  protected:
    Alcove_id value_;
    bool is_end_;
    std::vector<Alcove_id> chunks_;
    std::vector<std::size_t> decomposition_;
    std::vector<std::pair<Scs_coface_iterator, Scs_coface_iterator> > coface_iterators_;
    std::vector<Simple_coxeter_system_iterator> scs_iterators_;
  };
  
public:
  typedef boost::iterator_range<Coface_iterator> Coface_range;
  
  Coface_range coface_range(const Alcove_id& a_id, std::size_t k) const {
    return Coface_range(Coface_iterator(a_id, *this, k),
                        Coface_iterator(a_id, *this, dimension_ + 1));
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Face range experimental
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
private:

  class Face2_iterator : public boost::iterator_facade< Face2_iterator,
                                                        Alcove_id const,
                                                        boost::forward_traversal_tag> {
  protected:
    typedef typename std::vector<Simple_coxeter_system>::const_iterator Simple_coxeter_system_iterator;
    typedef typename Simple_coxeter_system::Face2_iterator Scs_face_iterator;
    friend class boost::iterator_core_access;
    
    void update_value(std::size_t first_change_) {
      if (is_end_)
        return;
      for (std::size_t i = first_change_; i < decomposition_.size(); ++i) {
        Alcove_id face = *face_iterators_[i].first;
        for (std::size_t j = 0; j < face.size(); ++j)
          value_.push_back(face[j], face.is_fixed(j));
      }
      first_change_ = decomposition_.size();
    }

    std::size_t update_from(std::size_t pos, std::size_t rest) {
      for (std::size_t i = pos; i < decomposition_.size(); ++i) {
        decomposition_[i] =
          (rest < chunks_[i].dimension() ? rest : chunks_[i].dimension());
        auto face_range = scs_iterators_[i]->face2_range(chunks_[i],
                                                         decomposition_[i]);
        face_iterators_[i] = std::make_pair(face_range.begin(), face_range.end());
        rest -= decomposition_[i];
      }
      return rest;
    }
    
    bool equal(Face2_iterator const& other) const {
      return (is_end_ && other.is_end_) || (!is_end_ && !other.is_end_ && decomposition_ == other.decomposition_);
    }
    Alcove_id const& dereference() const {
      return value_;
    }
    void increment() {
      if (is_end_)
        return;
      std::size_t rest = 0;
      std::size_t pos = decomposition_.size() - 1;
      while (true) {
        value_.resize(value_.size() - chunks_[pos].size());
        face_iterators_[pos].first++;
        if (face_iterators_[pos].first == face_iterators_[pos].second) {
          if (decomposition_[pos] == 0) {
            if (pos == 0) {
              is_end_ = true;
              return;
            }
            rest += decomposition_[pos];
            pos--;
            continue;
          }
          decomposition_[pos]--;
          rest++;    
          auto face_range = scs_iterators_[pos]->face2_range(chunks_[pos],
                                                             decomposition_[pos]);
          face_iterators_[pos] = std::make_pair(face_range.begin(), face_range.end());
        }
        if (update_from(pos + 1, rest)) {
          if (pos == 0) {
            is_end_ = true;
            return;
          }
          rest += decomposition_[pos];
          pos--;
          continue;
        }
        update_value(pos);
        return;
      }
    }
     
  public:
    Face2_iterator(const Alcove_id& coface,
                   const Coxeter_system& cs,
                   std::size_t value_dimension)
      : value_(coface.level(), value_dimension),
        is_end_(false),
        decomposition_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin()),
        face_iterators_(cs.simple_coxeter_system_end() - cs.simple_coxeter_system_begin())
    {
      std::size_t pos = 0;
      for (auto scs_it = cs.simple_coxeter_system_begin();
           scs_it != cs.simple_coxeter_system_end();
           ++scs_it) {
        scs_iterators_.push_back(scs_it);
        Alcove_id chunk(coface.level(), scs_it->dimension());
        for (std::size_t i = pos; i < pos + scs_it->pos_root_count(); ++i)
          chunk.push_back(coface[i], coface.is_fixed(i));
        chunk.set_dimension(scs_it->alcove_dimension(chunk));
        chunks_.push_back(chunk);
        pos += scs_it->pos_root_count();
      }
      if (update_from(0, value_dimension)) {
        is_end_ = true;
        return;
      }
      update_value(0);
    }

  protected:
    Alcove_id value_;
    bool is_end_;
    std::vector<Alcove_id> chunks_;
    std::vector<std::size_t> decomposition_;
    std::vector<std::pair<Scs_face_iterator, Scs_face_iterator> > face_iterators_;
    std::vector<Simple_coxeter_system_iterator> scs_iterators_;
  };

  
public:
  typedef boost::iterator_range<Face2_iterator> Face2_range;
  
  Face2_range face2_range(const Alcove_id& a_id, std::size_t k) const {
    return Face2_range(Face2_iterator(a_id, *this, k),
                       Face2_iterator(a_id, *this, dimension_ + 1));
  }  


  
private:

  int gcd(int a, int b) const {
    return b == 0 ? std::abs(a) : gcd(b, a % b);
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

  template <class Visitor>
  void rec_combine_chunks_alcove(std::vector<std::vector<Filtered_alcove > >::iterator chunks_it,
                                 std::vector<std::vector<Filtered_alcove > >::iterator chunks_end,
                                 std::vector<std::vector<std::vector<Vertex_id > > >::iterator chunks_v_it,
                                 Visitor& visitor,
                                 double eps,
                                 Filtered_alcove& a,
                                 std::vector<std::vector<Vertex_id > >& vertex_chunks) const {
    if (chunks_it == chunks_end) {
      std::vector<Vertex_id> vertices;
#ifdef CC_A_V_VISITORS
      Vertex_id v_id(a.id.level());
      rec_combine_chunks(vertex_chunks.begin(),
                         vertices,
                         v_id);
#endif
      visitor(a, vertices);
      return;
    }
    auto chunk_ = chunks_it->begin();
    auto chunk_v_ = chunks_v_it->begin();
    for (; chunk_ != chunks_it->end(); ++chunk_, ++chunk_v_) {
      // Filtered_alcove::iterator         c_it = chunk_->id.begin();
      // std::vector<Vertex_id >::iterator c_v_it = chunk_v_->begin();
      for (auto c: chunk_->id)
        a.id.push_back(c);
#ifdef CC_A_V_VISITORS
      vertex_chunks.push_back(*chunk_v_);
#endif
      a.f += chunk_->f;
      if (a.f <= eps*eps)
        rec_combine_chunks_alcove(chunks_it+1,
                                  chunks_end,
                                  chunks_v_it+1,
                                  visitor,
                                  eps,
                                  a,
                                  vertex_chunks);
      a.id.resize(a.id.size()-chunk_->id.size());
      a.f -= chunk_->f;
#ifdef CC_A_V_VISITORS
      vertex_chunks.pop_back();
#endif
    }
  }

  
  void rec_combine_chunks(std::vector<std::vector<Vertex_id>>::iterator chunks_it,
                          std::vector<Vertex_id>& vertices,
                          Vertex_id& v_id) const {
    int k = v_id.size();
    if (k == dimension_) {
      //      vertices.push_back(reduced_id(v_id));
      vertices.push_back(v_id);
      return;
    }
    for (auto chunk: *chunks_it) {
      for (auto c_it = chunk.begin(); c_it != chunk.end(); ++c_it)
        v_id.push_back(*c_it);
      rec_combine_chunks(chunks_it+1, vertices, v_id);
      v_id.resize(v_id.size()-chunk.size());
    }
  }
  
public:  
  
  std::vector<Vertex_id> vertices_of_alcove(const Alcove_id& ai_id) const
  {
    std::vector<Vertex_id> vertices;
    std::vector<std::vector<Vertex_id>> chunks;
    auto ai_it = ai_id.begin();
    double level = ai_id.level();
    for (auto scs: simple_system_range_) {
      Alcove_id ai_id_part(level);
      unsigned pos_root_count = scs.pos_root_count();
      for (unsigned i = 0; i < pos_root_count; i++) {
        ai_id_part.push_back(*ai_it++);
      }
      chunks.emplace_back(scs.vertices_of_simplex(ai_id_part));
    }
    // std::vector<std::vector<Vertex_id>::iterator> iterators;
    // for (auto chunk: chunks)
    //   iterators.emplace_back(chunk.begin());
    Vertex_id v_id(level);
    rec_combine_chunks(chunks.begin(), vertices, v_id);
    return vertices;
  }

  std::vector<double> barycenter(const Alcove_id& a_id) const {
    std::vector<double> result;
    auto a_it = a_id.begin();
    for (auto scs: simple_system_range_) {
      Alcove_id coordinate_segment(a_id.level());
      unsigned pos_root_count = scs.pos_root_count();
      coordinate_segment.reserve(pos_root_count);
      for (unsigned i = 0; i < pos_root_count; i++)
        coordinate_segment.push_back(*a_it++);
      std::vector<double> barycenter = scs.barycenter(coordinate_segment);
      result.insert(result.end(), barycenter.begin(), barycenter.end());
    }
    return result;
  }
  
  bool is_adjacent(const Vertex_id& v_id, const Alcove_id& a_id) const {
    int i = 1, j = 1;
    for (auto scs: simple_system_range_) {
      Vertex_id chunk_v(v_id.level());
      for (unsigned ii = 0; ii < scs.dimension(); i++, ii++)
        chunk_v.push_back(v_id[i]);
      Alcove_id chunk_a(a_id.level());
      for (unsigned jj = 0; jj < scs.pos_root_count(); j++, jj++)
        chunk_a.push_back(a_id[j]);
      if (!scs.is_adjacent(chunk_v, chunk_a))
        return false;
     }
    return true;
  }

  template <class VMap,
            class Simplex_range>
  void write_mesh(VMap& v_map, Simplex_range& range, std::string file_name = "toplex.mesh") const {
    if (simple_system_range_.size() == 1)
      if (simple_system_range_.begin()->dimension() == 2 || simple_system_range_.begin()->dimension() == 3)
        simple_system_range_.begin()->write_mesh(v_map, range, file_name);
  }

};


  // Print the Coxeter_system in os.
std::ostream& operator<<(std::ostream & os, Coxeter_system& cs) {
  os << "[";
  auto scs_it = cs.simple_coxeter_system_begin();
  if (scs_it != cs.simple_coxeter_system_end())
    os << *scs_it++;
  for (; scs_it != cs.simple_coxeter_system_end(); ++scs_it)
    os << ", " << *scs_it;
  os << "]" << std::endl;
  return os;
}


#endif
