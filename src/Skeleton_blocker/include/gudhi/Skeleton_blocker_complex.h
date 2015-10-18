/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SKELETON_BLOCKER_COMPLEX_H_
#define SKELETON_BLOCKER_COMPLEX_H_

#include <gudhi/Skeleton_blocker/iterators/Skeleton_blockers_iterators.h>
#include <gudhi/Skeleton_blocker_link_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_link_superior.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_simplex.h>

#include <gudhi/Skeleton_blocker/Skeleton_blocker_complex_visitor.h>
#include <gudhi/Skeleton_blocker/internal/Top_faces.h>
#include <gudhi/Skeleton_blocker/internal/Trie.h>

#include <gudhi/Utils.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/map.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

namespace Gudhi {

namespace skbl {

/**
 *@class Skeleton_blocker_complex
 *@brief Abstract Simplicial Complex represented with a skeleton/blockers pair.
 *@ingroup skbl
 */
template<class SkeletonBlockerDS>
class Skeleton_blocker_complex {
  template<class ComplexType> friend class Vertex_iterator;
  template<class ComplexType> friend class Neighbors_vertices_iterator;
  template<class ComplexType> friend class Edge_iterator;
  template<class ComplexType> friend class Edge_around_vertex_iterator;

  template<class ComplexType> friend class Skeleton_blocker_link_complex;
  template<class ComplexType> friend class Skeleton_blocker_link_superior;
  template<class ComplexType> friend class Skeleton_blocker_sub_complex;

 public:
  /**
   * @brief The type of stored vertex node, specified by the template SkeletonBlockerDS
   */
  typedef typename SkeletonBlockerDS::Graph_vertex Graph_vertex;

  /**
   * @brief The type of stored  edge node, specified by the template SkeletonBlockerDS
   */
  typedef typename SkeletonBlockerDS::Graph_edge Graph_edge;

  typedef typename SkeletonBlockerDS::Root_vertex_handle Root_vertex_handle;

  /**
   * @brief The type of an handle to a vertex of the complex.
   */
  typedef typename SkeletonBlockerDS::Vertex_handle Vertex_handle;
  typedef typename Root_vertex_handle::boost_vertex_handle boost_vertex_handle;

  /**
   * @brief A ordered set of integers that represents a simplex.
   */
  typedef Skeleton_blocker_simplex<Vertex_handle> Simplex;
  typedef Skeleton_blocker_simplex<Root_vertex_handle> Root_simplex_handle;

  /**
   * @brief Handle to a blocker of the complex.
   */
  typedef Simplex* Blocker_handle;

  typedef typename Root_simplex_handle::Simplex_vertex_const_iterator Root_simplex_iterator;
  typedef typename Simplex::Simplex_vertex_const_iterator Simplex_handle_iterator;

 protected:
  typedef typename boost::adjacency_list<boost::setS,  // edges
  boost::vecS,  // vertices
  boost::undirectedS, Graph_vertex, Graph_edge> Graph;
  // todo/remark : edges are not sorted, it heavily penalizes computation for SuperiorLink
  // (eg Link with greater vertices)
  // that burdens simplex iteration / complex initialization via list of simplices.
  // to avoid that, one should modify the graph by storing two lists of adjacency for every
  // vertex, the one with superior and the one with lower vertices, that way there is
  // no more extra cost for computation of SuperiorLink
  typedef typename boost::graph_traits<Graph>::vertex_iterator boost_vertex_iterator;
  typedef typename boost::graph_traits<Graph>::edge_iterator boost_edge_iterator;

 protected:
  typedef typename boost::graph_traits<Graph>::adjacency_iterator boost_adjacency_iterator;

 public:
  /**
   * @brief Handle to an edge of the complex.
   */
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_handle;

 protected:
  typedef std::multimap<Vertex_handle, Simplex *> BlockerMap;
  typedef typename std::multimap<Vertex_handle, Simplex *>::value_type BlockerPair;
  typedef typename std::multimap<Vertex_handle, Simplex *>::iterator BlockerMapIterator;
  typedef typename std::multimap<Vertex_handle, Simplex *>::const_iterator BlockerMapConstIterator;

 protected:
  size_t num_vertices_;
  size_t num_blockers_;

  typedef Skeleton_blocker_complex_visitor<Vertex_handle> Visitor;
  // typedef Visitor* Visitor_ptr;
  Visitor* visitor;

  /**
   * @details If 'x' is a Vertex_handle of a vertex in the complex then degree[x] = d is its degree.
   *
   * This quantity is updated when adding/removing edge.
   *
   * This is useful because the operation
   * list.size() is done in linear time.
   */
  std::vector<boost_vertex_handle> degree_;
  Graph skeleton; /** 1-skeleton of the simplicial complex. */

  /** Each vertex can access to the blockers passing through it. */
  BlockerMap blocker_map_;

 public:
  /////////////////////////////////////////////////////////////////////////////
  /** @name Constructors, Destructors
   */
  //@{

  /**
   *@brief constructs a simplicial complex with a given number of vertices and a visitor.
   */
  explicit Skeleton_blocker_complex(size_t num_vertices_ = 0, Visitor* visitor_ = NULL)
      : visitor(visitor_) {
    clear();
    for (size_t i = 0; i < num_vertices_; ++i) {
      add_vertex();
    }
  }

 private:
  // typedef Trie<Skeleton_blocker_complex<SkeletonBlockerDS>> STrie;
  typedef Trie<Simplex> STrie;

 public:
  /**
   * @brief Constructor with a list of simplices.
   * @details is_flag_complex indicates if the complex is a flag complex or not (to know if blockers have to be computed or not).
   */
  template<typename SimpleHandleOutputIterator>
  Skeleton_blocker_complex(SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end,
                           bool is_flag_complex = false, Visitor* visitor_ = NULL)
      : num_vertices_(0),
      num_blockers_(0),
      visitor(visitor_) {
    add_vertex_and_edges(simplex_begin, simplex_end);

    if (!is_flag_complex)
      // need to compute blockers
      add_blockers(simplex_begin, simplex_end);
  }

 private:
  template<typename SimpleHandleOutputIterator>
  void add_vertex_and_edges(SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end) {
    std::vector<std::pair<Vertex_handle, Vertex_handle>> edges;
    // first pass to add vertices and edges
    int num_vertex = -1;
    for (auto s_it = simplex_begin; s_it != simplex_end; ++s_it) {
      if (s_it->dimension() == 0) num_vertex = (std::max)(num_vertex, s_it->first_vertex().vertex);
      if (s_it->dimension() == 1) edges.emplace_back(s_it->first_vertex(), s_it->last_vertex());
    }
    while (num_vertex-- >= 0) add_vertex();

    for (const auto& e : edges)
      add_edge_without_blockers(e.first, e.second);
  }

  template<typename SimpleHandleOutputIterator>
  void add_blockers(SimpleHandleOutputIterator simplex_begin, SimpleHandleOutputIterator simplex_end) {
    Tries<Simplex> tries(num_vertices(), simplex_begin, simplex_end);
    tries.init_next_dimension();
    auto simplices(tries.next_dimension_simplices());

    while (!simplices.empty()) {
      simplices = tries.next_dimension_simplices();
      for (auto& sigma : simplices) {
        // common_positive_neighbors is the set of vertices u such that
        // for all s in sigma, us is an edge and u>s
        Simplex common_positive_neighbors(tries.positive_neighbors(sigma.last_vertex()));
        for (auto sigma_it = sigma.rbegin(); sigma_it != sigma.rend(); ++sigma_it)
          if (sigma_it != sigma.rbegin())
            common_positive_neighbors.intersection(tries.positive_neighbors(*sigma_it));

        for (auto x : common_positive_neighbors) {
          // first test that all edges sx are here for all s in sigma
          bool all_edges_here = true;
          for (auto s : sigma)
            if (!contains_edge(x, s)) {
              all_edges_here = false;
              break;
            }
          if (!all_edges_here) continue;

          // all edges of sigma \cup x are here
          // we have a blocker if all proper faces of sigma \cup x
          // are in the complex and if sigma \cup x is not in the complex
          // the first is equivalent at checking if blocks(sigma \cup x) is true
          // as all blockers of lower dimension have already been computed
          sigma.add_vertex(x);
          if (!tries.contains(sigma) && !blocks(sigma))
            add_blocker(sigma);
          sigma.remove_vertex(x);
        }
      }
    }
  }

 public:
  // We cannot use the default copy constructor since we need
  // to make a copy of each of the blockers

  Skeleton_blocker_complex(const Skeleton_blocker_complex& copy) {
    visitor = NULL;
    degree_ = copy.degree_;
    skeleton = Graph(copy.skeleton);
    num_vertices_ = copy.num_vertices_;

    num_blockers_ = 0;
    // we copy the blockers
    for (auto blocker : copy.const_blocker_range()) {
      add_blocker(*blocker);
    }
  }

  /**
   */
  Skeleton_blocker_complex& operator=(const Skeleton_blocker_complex& copy) {
    clear();
    visitor = NULL;
    degree_ = copy.degree_;
    skeleton = Graph(copy.skeleton);
    num_vertices_ = copy.num_vertices_;

    num_blockers_ = 0;
    // we copy the blockers
    for (auto blocker : copy.const_blocker_range())
      add_blocker(*blocker);
    return *this;
  }

  /**
   * return true if both complexes have the same simplices.
   */
  bool operator==(const Skeleton_blocker_complex& other) const {
    if (other.num_vertices() != num_vertices()) return false;
    if (other.num_edges() != num_edges()) return false;
    if (other.num_blockers() != num_blockers()) return false;

    for (auto v : vertex_range())
      if (!other.contains_vertex(v)) return false;

    for (auto e : edge_range())
      if (!other.contains_edge(first_vertex(e), second_vertex(e))) return false;

    for (const auto b : const_blocker_range())
      if (!other.contains_blocker(*b)) return false;

    return true;
  }

  bool operator!=(const Skeleton_blocker_complex& other) const {
    return !(*this == other);
  }

  /**
   * The destructor delete all blockers allocated.
   */
  virtual ~Skeleton_blocker_complex() {
    clear();
  }

  /**
   * @details Clears the simplicial complex. After a call to this function,
   * blockers are destroyed. The 1-skeleton and the set of blockers
   * are both empty.
   */
  virtual void clear() {
    // xxx for now the responsabilty of freeing the visitor is for
    // the user
    visitor = NULL;

    degree_.clear();
    num_vertices_ = 0;

    remove_blockers();

    skeleton.clear();
  }

  /**
   *@brief allows to change the visitor.
   */
  void set_visitor(Visitor* other_visitor) {
    visitor = other_visitor;
  }

  //@}

  /////////////////////////////////////////////////////////////////////////////
  /** @name Vertices operations
   */
  //@{
 public:
  /**
   * @brief Return a local Vertex_handle of a vertex given a global one.
   * @remark Assume that the vertex is present in the complex.
   */
  Vertex_handle operator[](Root_vertex_handle global) const {
    auto local(get_address(global));
    assert(local);
    return *local;
  }

  /**
   * @brief Return the vertex node associated to local Vertex_handle.
   * @remark Assume that the vertex is present in the complex.
   */
  Graph_vertex& operator[](Vertex_handle address) {
    assert(
           0 <= address.vertex && address.vertex < boost::num_vertices(skeleton));
    return skeleton[address.vertex];
  }

  /**
   * @brief Return the vertex node associated to local Vertex_handle.
   * @remark Assume that the vertex is present in the complex.
   */
  const Graph_vertex& operator[](Vertex_handle address) const {
    assert(
           0 <= address.vertex && address.vertex < boost::num_vertices(skeleton));
    return skeleton[address.vertex];
  }

  /**
   * @brief Adds a vertex to the simplicial complex and returns its Vertex_handle.
   */
  Vertex_handle add_vertex() {
    Vertex_handle address(boost::add_vertex(skeleton));
    num_vertices_++;
    (*this)[address].activate();
    // safe since we now that we are in the root complex and the field 'address' and 'id'
    // are identical for every vertices
    (*this)[address].set_id(Root_vertex_handle(address.vertex));
    degree_.push_back(0);
    if (visitor)
      visitor->on_add_vertex(address);
    return address;
  }

  /**
   * @brief Remove a vertex from the simplicial complex
   * @remark It just deactivates the vertex with a boolean flag but does not
   * remove it from vertices from complexity issues.
   */
  void remove_vertex(Vertex_handle address) {
    assert(contains_vertex(address));
    // We remove b
    boost::clear_vertex(address.vertex, skeleton);
    (*this)[address].deactivate();
    num_vertices_--;
    degree_[address.vertex] = -1;
    if (visitor)
      visitor->on_remove_vertex(address);
  }

  /**
   */
  bool contains_vertex(Vertex_handle u) const {
    if (u.vertex < 0 || u.vertex >= boost::num_vertices(skeleton))
      return false;
    return (*this)[u].is_active();
  }

  /**
   */
  bool contains_vertex(Root_vertex_handle u) const {
    boost::optional<Vertex_handle> address = get_address(u);
    return address && (*this)[*address].is_active();
  }

  /**
   * @return true iff the simplicial complex contains all vertices
   * of simplex sigma
   */
  bool contains_vertices(const Simplex & sigma) const {
    for (auto vertex : sigma)
      if (!contains_vertex(vertex))
        return false;
    return true;
  }

  /**
   * @brief Given an Id return the address of the vertex having this Id in the complex.
   * @remark For a simplicial complex, the address is the id but it may not be the case for a SubComplex.
   */
  virtual boost::optional<Vertex_handle> get_address(
                                                     Root_vertex_handle id) const {
    boost::optional<Vertex_handle> res;
    if (id.vertex < boost::num_vertices(skeleton))
      res = Vertex_handle(id.vertex);  // xxx
    return res;
  }

  /**
   * return the id of a vertex of adress local present in the graph
   */
  Root_vertex_handle get_id(Vertex_handle local) const {
    assert(0 <= local.vertex && local.vertex < boost::num_vertices(skeleton));
    return (*this)[local].get_id();
  }

  /**
   * @brief Convert an address of a vertex of a complex to the address in
   * the current complex.
   * @details
   * If the current complex is a sub (or sup) complex of 'other', it converts
   * the address of a vertex v expressed in 'other' to the address of the vertex
   * v in the current one.
   * @remark this methods uses Root_vertex_handle to identify the vertex and
   * assumes the vertex is present in the current complex.
   */
  Vertex_handle convert_handle_from_another_complex(const Skeleton_blocker_complex& other,
                                                    Vertex_handle vh_in_other) const {
    auto vh_in_current_complex = get_address(other.get_id(vh_in_other));
    assert(vh_in_current_complex);
    return *vh_in_current_complex;
  }

  /**
   * @brief return the graph degree of a vertex.
   */
  int degree(Vertex_handle local) const {
    assert(0 <= local.vertex && local.vertex < boost::num_vertices(skeleton));
    return degree_[local.vertex];
  }

  //@}

  /////////////////////////////////////////////////////////////////////////////
  /** @name Edges operations
   */
  //@{
 public:
  /**
   * @brief return an edge handle if the two vertices forms
   * an edge in the complex
   */
  boost::optional<Edge_handle> operator[](
                                          const std::pair<Vertex_handle, Vertex_handle>& ab) const {
    boost::optional<Edge_handle> res;
    std::pair<Edge_handle, bool> edge_pair(
                                           boost::edge(ab.first.vertex, ab.second.vertex, skeleton));
    if (edge_pair.second)
      res = edge_pair.first;
    return res;
  }

  /**
   * @brief returns the stored node associated to an edge
   */
  Graph_edge& operator[](Edge_handle edge_handle) {
    return skeleton[edge_handle];
  }

  /**
   * @brief returns the stored node associated to an edge
   */
  const Graph_edge& operator[](Edge_handle edge_handle) const {
    return skeleton[edge_handle];
  }

  /**
   * @brief returns the first vertex of an edge
   * @details it assumes that the edge is present in the complex
   */
  Vertex_handle first_vertex(Edge_handle edge_handle) const {
    return static_cast<Vertex_handle> (source(edge_handle, skeleton));
  }

  /**
   * @brief returns the first vertex of an edge
   * @details it assumes that the edge is present in the complex
   */
  Vertex_handle second_vertex(Edge_handle edge_handle) const {
    return static_cast<Vertex_handle> (target(edge_handle, skeleton));
  }

  /**
   * @brief returns the simplex made with the two vertices of the edge
   * @details it assumes that the edge is present in the complex

   */
  Simplex get_vertices(Edge_handle edge_handle) const {
    auto edge((*this)[edge_handle]);
    return Simplex((*this)[edge.first()], (*this)[edge.second()]);
  }

  /**
   * @brief Adds an edge between vertices a and b.
   * @details For instance, the complex contains edge 01 and 12, then calling
   * add_edge with vertex 0 and 2 will create a complex containing
   * the edges 01, 12, 20 but not the triangle 012 (and hence this complex
   * will contains a blocker 012).
   */
  Edge_handle add_edge(Vertex_handle a, Vertex_handle b) {    
    //if the edge is already there we musnt go further
    //as we may add blockers that should not be here
    if(contains_edge(a,b)) 
      return *((*this)[std::make_pair(a,b)]);
    auto res = add_edge_without_blockers(a,b);
    add_blockers_after_simplex_insertion(Simplex(a,b));
    return res;
  }

    /**
   * @brief Adds all edges of s in the complex.
   */
  void add_edge(const Simplex& s) {
    for(auto i = s.begin(); i != s.end(); ++i)
      for(auto j = i; ++j != s.end(); /**/)
        add_edge(*i,*j);
  }

  /**
   * @brief Adds an edge between vertices a and b without blockers.
   * @details For instance, the complex contains edge 01 and 12, then calling
   * add_edge with vertex 0 and 2 will create a complex containing
   * the triangle 012.
   */
  Edge_handle add_edge_without_blockers(Vertex_handle a, Vertex_handle b) {
    assert(contains_vertex(a) && contains_vertex(b));
    assert(a != b);

    auto edge_handle((*this)[std::make_pair(a, b)]);
    if (!edge_handle) {
      edge_handle = boost::add_edge(a.vertex, b.vertex, skeleton).first;
      (*this)[*edge_handle].setId(get_id(a), get_id(b));
      degree_[a.vertex]++;
      degree_[b.vertex]++;
      if (visitor)
        visitor->on_add_edge_without_blockers(a, b);
    }
    return *edge_handle;
  }


  /**
   * @brief Adds all edges of s in the complex without adding blockers.
   */
  void add_edge_without_blockers(Simplex s) {
    for(auto i = s.begin(); i != s.end(); ++i){
      for(auto j = i; ++j != s.end(); /**/)
        add_edge_without_blockers(*i,*j);
    }
  }

  /**
   * @brief Removes an edge from the simplicial complex and all its cofaces.
   * @details returns the former Edge_handle representing the edge
   */
  virtual Edge_handle remove_edge(Vertex_handle a, Vertex_handle b) {
    bool found;
    Edge_handle edge;
    tie(edge, found) = boost::edge(a.vertex, b.vertex, skeleton);
    if (found) {
      if (visitor)
        visitor->on_remove_edge(a, b);
      boost::remove_edge(a.vertex, b.vertex, skeleton);
      degree_[a.vertex]--;
      degree_[b.vertex]--;
    }
    return edge;
  }

  /**
   * @brief Removes edge and its cofaces from the simplicial complex.
   */
  void remove_edge(Edge_handle edge) {
    assert(contains_vertex(first_vertex(edge)));
    assert(contains_vertex(second_vertex(edge)));
    remove_edge(first_vertex(edge), second_vertex(edge));
  }

  /**
   * @brief The complex is reduced to its set of vertices.
   * All the edges and blockers are removed.
   */
  void keep_only_vertices() {
    remove_blockers();

    for (auto u : vertex_range()) {
      while (this->degree(u) > 0) {
        Vertex_handle v(*(adjacent_vertices(u.vertex, this->skeleton).first));
        this->remove_edge(u, v);
      }
    }
  }

  /**
   * @return true iff the simplicial complex contains an edge between
   * vertices a and b
   */
  bool contains_edge(Vertex_handle a, Vertex_handle b) const {
    // if (a.vertex<0 || b.vertex <0) return false;
    return boost::edge(a.vertex, b.vertex, skeleton).second;
  }

  /**
   * @return true iff the simplicial complex contains all vertices
   * and all edges of simplex sigma
   */
  bool contains_edges(const Simplex & sigma) const {
    for (auto i = sigma.begin(); i != sigma.end(); ++i) {
      if (!contains_vertex(*i))
        return false;
      for (auto j = i; ++j != sigma.end();) {
        if (!contains_edge(*i, *j))
          return false;
      }
    }
    return true;
  }
  //@}

  /////////////////////////////////////////////////////////////////////////////
  /** @name Blockers operations
   */
  //@{

  /**
   * @brief Adds the simplex to the set of blockers and
   * returns a Blocker_handle toward it if was not present before and 0 otherwise.
   */
  Blocker_handle add_blocker(const Simplex& blocker) {
    assert(blocker.dimension() > 1);
    if (contains_blocker(blocker)) {
      return 0;
    } else {
      if (visitor)
        visitor->on_add_blocker(blocker);
      Blocker_handle blocker_pt = new Simplex(blocker);
      num_blockers_++;
      auto vertex = blocker_pt->begin();
      while (vertex != blocker_pt->end()) {
        blocker_map_.insert(BlockerPair(*vertex, blocker_pt));
        ++vertex;
      }
      return blocker_pt;
    }
  }

 protected:
  /**
   * @brief Adds the simplex to the set of blockers
   */
  void add_blocker(Blocker_handle blocker) {
    if (contains_blocker(*blocker)) {
      // std::cerr << "ATTEMPT TO ADD A BLOCKER ALREADY THERE ---> BLOCKER IGNORED" << endl;
      return;
    } else {
      if (visitor)
        visitor->on_add_blocker(*blocker);
      num_blockers_++;
      auto vertex = blocker->begin();
      while (vertex != blocker->end()) {
        blocker_map_.insert(BlockerPair(*vertex, blocker));
        ++vertex;
      }
    }
  }

 protected:
  /**
   * Removes sigma from the blocker map of vertex v
   */
  void remove_blocker(const Blocker_handle sigma, Vertex_handle v) {
    Complex_blocker_around_vertex_iterator blocker;
    for (blocker = blocker_range(v).begin(); blocker != blocker_range(v).end();
         ++blocker) {
      if (*blocker == sigma)
        break;
    }
    if (*blocker != sigma) {
      std::cerr
          << "bug ((*blocker).second == sigma) ie try to remove a blocker not present\n";
      assert(false);
    } else {
      blocker_map_.erase(blocker.current_position());
    }
  }

 public:
  /**
   * @brief Removes the simplex from the set of blockers.
   * @remark sigma has to belongs to the set of blockers
   */
  void remove_blocker(const Blocker_handle sigma) {
    for (auto vertex : *sigma)
      remove_blocker(sigma, vertex);
    num_blockers_--;
  }

  /**
   * @brief Remove all blockers, in other words, it expand the simplicial
   * complex to the smallest flag complex that contains it.
   */
  void remove_blockers() {
    // Desallocate the blockers
    while (!blocker_map_.empty()) {
      delete_blocker(blocker_map_.begin()->second);
    }
    num_blockers_ = 0;
    blocker_map_.clear();
  }

 protected:
  /**
   * Removes the simplex sigma from the set of blockers.
   * sigma has to belongs to the set of blockers
   *
   * @remark contrarily to delete_blockers does not call the destructor
   */
  void remove_blocker(const Simplex& sigma) {
    assert(contains_blocker(sigma));
    for (auto vertex : sigma)
      remove_blocker(sigma, vertex);
    num_blockers_--;
  }

 public:
  /**
   * Removes the simplex s from the set of blockers
   * and desallocate s.
   */
  void delete_blocker(Blocker_handle sigma) {
    if (visitor)
      visitor->on_delete_blocker(sigma);
    remove_blocker(sigma);
    delete sigma;
  }

  /**
   * @return true iff s is a blocker of the simplicial complex
   */
  bool contains_blocker(const Blocker_handle s) const {
    if (s->dimension() < 2)
      return false;

    Vertex_handle a = s->first_vertex();

    for (const auto blocker : const_blocker_range(a)) {
      if (s == *blocker)
        return true;
    }
    return false;
  }

  /**
   * @return true iff s is a blocker of the simplicial complex
   */
  bool contains_blocker(const Simplex & s) const {
    if (s.dimension() < 2)
      return false;

    Vertex_handle a = s.first_vertex();

    for (auto blocker : const_blocker_range(a)) {
      if (s == *blocker)
        return true;
    }
    return false;
  }

 private:
  /**
   * @return true iff a blocker of the simplicial complex
   * is a face of sigma.
   */
  bool blocks(const Simplex & sigma) const {
    for (auto s : sigma)
      for (auto blocker : const_blocker_range(s))
        if (sigma.contains(*blocker))
          return true;
    return false;
  }

  //@}

 protected:
  /**
   * @details Adds to simplex the neighbours of v e.g. \f$ n \leftarrow n \cup N(v) \f$.
   * If keep_only_superior is true then only vertices that are greater than v are added.
   */
  virtual void add_neighbours(Vertex_handle v, Simplex & n,
                              bool keep_only_superior = false) const {
    boost_adjacency_iterator ai, ai_end;
    for (tie(ai, ai_end) = adjacent_vertices(v.vertex, skeleton); ai != ai_end;
         ++ai) {
      if (keep_only_superior) {
        if (*ai > v.vertex) {
          n.add_vertex(Vertex_handle(*ai));
        }
      } else {
        n.add_vertex(Vertex_handle(*ai));
      }
    }
  }

  /**
   * @details Add to simplex res all vertices which are
   * neighbours of alpha: ie \f$ res \leftarrow res \cup N(alpha) \f$.
   *
   * If 'keep_only_superior' is true then only vertices that are greater than alpha are added.
   * todo revoir
   *
   */
  virtual void add_neighbours(const Simplex &alpha, Simplex & res,
                              bool keep_only_superior = false) const {
    res.clear();
    auto alpha_vertex = alpha.begin();
    add_neighbours(*alpha_vertex, res, keep_only_superior);
    for (alpha_vertex = (alpha.begin())++; alpha_vertex != alpha.end();
         ++alpha_vertex)
      keep_neighbours(*alpha_vertex, res, keep_only_superior);
  }

  /**
   * @details Remove from simplex n all vertices which are
   * not neighbours of v e.g. \f$ res \leftarrow res \cap N(v) \f$.
   * If 'keep_only_superior' is true then only vertices that are greater than v are keeped.
   */
  virtual void keep_neighbours(Vertex_handle v, Simplex& res,
                               bool keep_only_superior = false) const {
    Simplex nv;
    add_neighbours(v, nv, keep_only_superior);
    res.intersection(nv);
  }

  /**
   * @details Remove from simplex all vertices which are
   * neighbours of v eg \f$ res \leftarrow res \setminus N(v) \f$.
   * If 'keep_only_superior' is true then only vertices that are greater than v are added.
   */
  virtual void remove_neighbours(Vertex_handle v, Simplex & res,
                                 bool keep_only_superior = false) const {
    Simplex nv;
    add_neighbours(v, nv, keep_only_superior);
    res.difference(nv);
  }

 public:
  typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex> Link_complex;

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Link_complex link(Vertex_handle v) const {
    return Link_complex(*this, Simplex(v));
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Link_complex link(Edge_handle edge) const {
    return Link_complex(*this, edge);
  }

  /**
   * Constructs the link of 'simplex' with points coordinates.
   */
  Link_complex link(const Simplex& simplex) const {
    return Link_complex(*this, simplex);
  }

  /**
   * @brief Compute the local vertices of 's' in the current complex
   * If one of them is not present in the complex then the return value is uninitialized.
   *
   *
   */
  // xxx rename get_address et place un using dans sub_complex

  boost::optional<Simplex> get_simplex_address(
                                                      const Root_simplex_handle& s) const {
    boost::optional<Simplex> res;

    Simplex s_address;
    // Root_simplex_const_iterator i;
    for (auto i = s.begin(); i != s.end(); ++i) {
      boost::optional<Vertex_handle> address = get_address(*i);
      if (!address)
        return res;
      else
        s_address.add_vertex(*address);
    }
    res = s_address;
    return res;
  }

  /**
   * @brief returns a simplex with vertices which are the id of vertices of the
   * argument.
   */
  Root_simplex_handle get_id(const Simplex& local_simplex) const {
    Root_simplex_handle global_simplex;
    for (auto x = local_simplex.begin(); x != local_simplex.end(); ++x) {
      global_simplex.add_vertex(get_id(*x));
    }
    return global_simplex;
  }

  /**
   * @brief returns true iff the simplex s belongs to the simplicial
   * complex.
   */
  virtual bool contains(const Simplex & s) const {
    if (s.dimension() == -1) {
      return false;
    } else if (s.dimension() == 0) {
      return contains_vertex(s.first_vertex());
    } else {
      return (contains_edges(s) && !blocks(s));
    }
  }

  /*
   * @brief returnrs true iff the complex is empty.
   */
  bool empty() const {
    return num_vertices() == 0;
  }

  /*
   * @brief returns the number of vertices in the complex.
   */
  int num_vertices() const {
    // remark boost::num_vertices(skeleton) counts deactivated vertices
    return num_vertices_;
  }

  /*
   * @brief returns the number of edges in the complex.
   * @details currently in O(n)
   */
  // todo cache the value

  int num_edges() const {
    return boost::num_edges(skeleton);
  }

  int num_triangles() const {
    auto triangles = triangle_range();
    return std::distance(triangles.begin(), triangles.end());
  }


  /*
   * @brief returns the number of simplices of a given dimension in the complex.
   */  
  size_t num_simplices() const {
    auto simplices = complex_simplex_range();
    return std::distance(simplices.begin(), simplices.end());
  }

  /*
   * @brief returns the number of simplices of a given dimension in the complex.
   */  
  size_t num_simplices(unsigned dimension) const {
    //todo iterator on k-simplices
    size_t res = 0;
    for(const auto& s: complex_simplex_range()) 
      if(s.dimension() == dimension) 
        ++res;
    return res;
  }

  /*
   * @brief returns the number of blockers in the complex.
   */
  size_t num_blockers() const {
    return num_blockers_;
  }

  /*
   * @brief returns true iff the graph of the 1-skeleton of the complex is complete.
   */
  bool complete() const {
    return (num_vertices() * (num_vertices() - 1)) / 2 == num_edges();
  }

  /**
   * @brief returns the number of connected components in the graph of the 1-skeleton.
   */
  int num_connected_components() const {
    int num_vert_collapsed = skeleton.vertex_set().size() - num_vertices();
    std::vector<int> component(skeleton.vertex_set().size());
    return boost::connected_components(this->skeleton, &component[0])
        - num_vert_collapsed;
  }

  /**
   * @brief %Test if the complex is a cone.
   * @details Runs in O(n) where n is the number of vertices.
   */
  bool is_cone() const {
    if (num_vertices() == 0)
      return false;
    if (num_vertices() == 1)
      return true;
    for (auto vi : vertex_range()) {
      // xxx todo faire une methode bool is_in_blocker(Vertex_handle)
      if (blocker_map_.find(vi) == blocker_map_.end()) {
        // no blocker passes through the vertex, we just need to
        // check if the current vertex is linked to all others vertices of the complex
        if (degree_[vi.vertex] == num_vertices() - 1)
          return true;
      }
    }
    return false;
  }

  //@}
  /** @Simplification operations
   */
  //@{

  /**
   * Returns true iff the blocker 'sigma' is popable.
   * To define popable, let us call 'L' the complex that
   * consists in the current complex without the blocker 'sigma'.
   * A blocker 'sigma' is then "popable" if the link of 'sigma'
   * in L is reducible.
   *
   */
  bool is_popable_blocker(Blocker_handle sigma) const;

  /**
   * Removes all the popable blockers of the complex and delete them.
   * @returns the number of popable blockers deleted
   */
  void remove_popable_blockers();

  /**
   * Removes all the popable blockers of the complex passing through v and delete them.
   */
  void remove_popable_blockers(Vertex_handle v);

  /**
   * @brief Removes all the popable blockers of the complex passing through v and delete them.
   * Also remove popable blockers in the neighborhood if they became popable.
   *
   */
  void remove_all_popable_blockers(Vertex_handle v);

  /**
   * Remove the star of the vertex 'v'
   */
  void remove_star(Vertex_handle v);

 private:
  /**
   * after removing the star of a simplex, blockers sigma that contains this simplex must be removed.
   * Furthermore, all simplices tau of the form sigma \setminus simplex_to_be_removed must be added
   * whenever the dimension of tau is at least 2.
   */
  void update_blockers_after_remove_star_of_vertex_or_edge(const Simplex& simplex_to_be_removed);

 public:
  /**
   * Remove the star of the edge connecting vertices a and b.
   * @returns the number of blocker that have been removed
   */
  void remove_star(Vertex_handle a, Vertex_handle b);

  /**
   * Remove the star of the edge 'e'.
   */
  void remove_star(Edge_handle e);

  /**
   * Remove the star of the simplex 'sigma' which needs to belong to the complex
   */
  void remove_star(const Simplex& sigma);

  /**
   * @brief add a simplex.
   * @details the simplex must have dimension greater than one (otherwise use add_vertex or add_edge_without_blockers).
   * and all vertices lower than the higher vertex of sigma must already be in the complex.
   * if some edges of sigma are not in the complex, then insert_edges_of_sigma flag must be 
   * set to true.
   */
  void add_simplex(const Simplex& sigma, bool insert_edges_of_sigma = false);

 private:
  void add_blockers_after_simplex_insertion(Simplex s);

  /**
   * remove all blockers that contains sigma
   */
  void remove_blocker_containing_simplex(const Simplex& sigma);

  /**
   * remove all blockers that contains sigma
   */
  void remove_blocker_include_in_simplex(const Simplex& sigma);

 public:
  enum simplifiable_status {
    NOT_HOMOTOPY_EQ, MAYBE_HOMOTOPY_EQ, HOMOTOPY_EQ
  };

  simplifiable_status is_remove_star_homotopy_preserving(const Simplex& simplex) {
    // todo write a virtual method 'link' in Skeleton_blocker_complex which will be overloaded by the current one of
    // Skeleton_blocker_geometric_complex
    // then call it there to build the link and return the value of link.is_contractible()
    return MAYBE_HOMOTOPY_EQ;
  }

  enum contractible_status {
    NOT_CONTRACTIBLE, MAYBE_CONTRACTIBLE, CONTRACTIBLE
  };

  /**
   * @brief %Test if the complex is reducible using a strategy defined in the class
   * (by default it tests if the complex is a cone)
   * @details Note that NO could be returned if some invariant ensures that the complex
   * is not a point (for instance if the euler characteristic is different from 1).
   * This function will surely have to return MAYBE in some case because the
   * associated problem is undecidable but it in practice, it can often
   * be solved with the help of geometry.
   */
  virtual contractible_status is_contractible() const {
    if (this->is_cone()) {
      return CONTRACTIBLE;
    } else {
      return MAYBE_CONTRACTIBLE;
    }
  }
  //@}

  /** @Edge contraction operations
   */
  //@{

  /**
   * @return If ignore_popable_blockers is true
   * then the result is true iff the link condition at edge ab is satisfied
   * or equivalently iff no blocker contains ab.
   * If ignore_popable_blockers is false then the
   * result is true iff all blocker containing ab are popable.
   */
  bool link_condition(Vertex_handle a, Vertex_handle b, bool ignore_popable_blockers = false) const {
    for (auto blocker : this->const_blocker_range(a))
      if (blocker->contains(b)) {
        // false if ignore_popable_blockers is false
        // otherwise the blocker has to be popable
        return ignore_popable_blockers && is_popable_blocker(blocker);
      }
    return true;
  }

  /**
   * @return If ignore_popable_blockers is true
   * then the result is true iff the link condition at edge ab is satisfied
   * or equivalently iff no blocker contains ab.
   * If ignore_popable_blockers is false then the
   * result is true iff all blocker containing ab are popable.
   */
  bool link_condition(Edge_handle e, bool ignore_popable_blockers = false) const {
    const Graph_edge& edge = (*this)[e];
    assert(this->get_address(edge.first()));
    assert(this->get_address(edge.second()));
    Vertex_handle a(*this->get_address(edge.first()));
    Vertex_handle b(*this->get_address(edge.second()));
    return link_condition(a, b, ignore_popable_blockers);
  }

 protected:
  /**
   * Compute simplices beta such that a.beta is an order 0 blocker
   * that may be used to construct a new blocker after contracting ab.
   * It requires that the link condition is satisfied.
   */
  void tip_blockers(Vertex_handle a, Vertex_handle b, std::vector<Simplex> & buffer) const;

 private:
  /**
   * @brief "Replace" the edge 'bx' by the edge 'ax'.
   * Assume that the edge 'bx' was present whereas 'ax' was not.
   * Precisely, it does not replace edges, but remove 'bx' and then add 'ax'.
   * The visitor 'on_swaped_edge' is called just after edge 'ax' had been added
   * and just before edge 'bx' had been removed. That way, it can
   * eventually access to information of 'ax'.
   */
  void swap_edge(Vertex_handle a, Vertex_handle b, Vertex_handle x);

 private:
  /**
   * @brief removes all blockers passing through the edge 'ab'
   */
  void delete_blockers_around_vertex(Vertex_handle v);

  /**
   * @brief removes all blockers passing through the edge 'ab'
   */
  void delete_blockers_around_edge(Vertex_handle a, Vertex_handle b);

 public:
  /**
   * Contracts the edge.
   * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
   * it removes first all blockers passing through 'ab'
   */
  void contract_edge(Edge_handle edge) {
    contract_edge(this->first_vertex(edge), this->second_vertex(edge));
  }

  /**
   * Contracts the edge connecting vertices a and b.
   * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
   * it removes first all blockers passing through 'ab'
   */
  void contract_edge(Vertex_handle a, Vertex_handle b);

 private:
  void get_blockers_to_be_added_after_contraction(Vertex_handle a, Vertex_handle b,
                                                  std::set<Simplex>& blockers_to_add);
  /**
   * delete all blockers that passes through a or b
   */
  void delete_blockers_around_vertices(Vertex_handle a, Vertex_handle b);
  void update_edges_after_contraction(Vertex_handle a, Vertex_handle b);
  void notify_changed_edges(Vertex_handle a);
  //@}

 public:
  /////////////////////////////////////////////////////////////////////////////
  /** @name Vertex iterators
   */
  //@{
  typedef Vertex_iterator<Skeleton_blocker_complex> Complex_vertex_iterator;

  //
  // Range over the vertices of the simplicial complex.
  // Methods .begin() and .end() return a Complex_vertex_iterator.
  //
  typedef boost::iterator_range<Complex_vertex_iterator> Complex_vertex_range;

  /**
   * @brief Returns a Complex_vertex_range over all vertices of the complex
   */
  Complex_vertex_range vertex_range() const {
    auto begin = Complex_vertex_iterator(this);
    auto end = Complex_vertex_iterator(this, 0);
    return Complex_vertex_range(begin, end);
  }

  typedef Neighbors_vertices_iterator<Skeleton_blocker_complex> Complex_neighbors_vertices_iterator;


  typedef boost::iterator_range<Complex_neighbors_vertices_iterator> Complex_neighbors_vertices_range;

  /**
   * @brief Returns a Complex_edge_range over all edges of the simplicial complex that passes trough v
   */
  Complex_neighbors_vertices_range vertex_range(Vertex_handle v) const {
    auto begin = Complex_neighbors_vertices_iterator(this, v);
    auto end = Complex_neighbors_vertices_iterator(this, v, 0);
    return Complex_neighbors_vertices_range(begin, end);
  }

  //@}

  /** @name Edge iterators
   */
  //@{

  typedef Edge_iterator<Skeleton_blocker_complex> Complex_edge_iterator;


  typedef boost::iterator_range<Complex_edge_iterator> Complex_edge_range;

  /**
   * @brief Returns a Complex_edge_range over all edges of the simplicial complex
   */
  Complex_edge_range edge_range() const {
    auto begin = Complex_edge_iterator(this);
    auto end = Complex_edge_iterator(this, 0);
    return Complex_edge_range(begin, end);
  }


  typedef Edge_around_vertex_iterator<Skeleton_blocker_complex> Complex_edge_around_vertex_iterator;


  typedef boost::iterator_range <Complex_edge_around_vertex_iterator> Complex_edge_around_vertex_range;

  /**
   * @brief Returns a Complex_edge_range over all edges of the simplicial complex that passes
   * through 'v'
   */
  Complex_edge_around_vertex_range edge_range(Vertex_handle v) const {
    auto begin = Complex_edge_around_vertex_iterator(this, v);
    auto end = Complex_edge_around_vertex_iterator(this, v, 0);
    return Complex_edge_around_vertex_range(begin, end);
  }

  //@}

  /** @name Triangles iterators
   */
  //@{
 private:
  typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<SkeletonBlockerDS> > Link;
  typedef Skeleton_blocker_link_superior<Skeleton_blocker_complex<SkeletonBlockerDS> > Superior_link;

 public:
  typedef Triangle_around_vertex_iterator<Skeleton_blocker_complex, Superior_link>
  Superior_triangle_around_vertex_iterator;
  typedef boost::iterator_range < Triangle_around_vertex_iterator<Skeleton_blocker_complex, Link> >
  Complex_triangle_around_vertex_range;

  /**
   * @brief Range over triangles around a vertex of the simplicial complex.
   * Methods .begin() and .end() return a Triangle_around_vertex_iterator.
   *
   */
  Complex_triangle_around_vertex_range triangle_range(Vertex_handle v) const {
    auto begin = Triangle_around_vertex_iterator<Skeleton_blocker_complex, Link>(this, v);
    auto end = Triangle_around_vertex_iterator<Skeleton_blocker_complex, Link>(this, v, 0);
    return Complex_triangle_around_vertex_range(begin, end);
  }

  typedef boost::iterator_range<Triangle_iterator<Skeleton_blocker_complex> > Complex_triangle_range;
  typedef Triangle_iterator<Skeleton_blocker_complex> Complex_triangle_iterator;

  /**
   * @brief Range over triangles of the simplicial complex.
   * Methods .begin() and .end() return a Triangle_around_vertex_iterator.
   *
   */
  Complex_triangle_range triangle_range() const {
    auto end = Triangle_iterator<Skeleton_blocker_complex>(this, 0);
    if (empty()) {
      return Complex_triangle_range(end, end);
    } else {
      auto begin = Triangle_iterator<Skeleton_blocker_complex>(this);
      return Complex_triangle_range(begin, end);
    }
  }

  //@}

  /** @name Simplices iterators
   */
  //@{
  typedef Simplex_around_vertex_iterator<Skeleton_blocker_complex, Link> Complex_simplex_around_vertex_iterator;

  /**
   * @brief Range over the simplices of the simplicial complex around a vertex.
   * Methods .begin() and .end() return a Complex_simplex_around_vertex_iterator.
   */
  typedef boost::iterator_range < Complex_simplex_around_vertex_iterator > Complex_simplex_around_vertex_range;

  /**
   * @brief Returns a Complex_simplex_around_vertex_range over all the simplices around a vertex of the complex
   */
  Complex_simplex_around_vertex_range star_simplex_range(Vertex_handle v) const {
    assert(contains_vertex(v));
    return Complex_simplex_around_vertex_range(
                                               Complex_simplex_around_vertex_iterator(this, v),
                                               Complex_simplex_around_vertex_iterator(this, v, true));
  }
  typedef Simplex_coboundary_iterator<Skeleton_blocker_complex, Link> Complex_simplex_coboundary_iterator;

  /**
   * @brief Range over the simplices of the coboundary of a simplex.
   * Methods .begin() and .end() return a Complex_simplex_coboundary_iterator.
   */
  typedef boost::iterator_range < Complex_simplex_coboundary_iterator > Complex_coboundary_range;

  /**
   * @brief Returns a Complex_simplex_coboundary_iterator over the simplices of the coboundary of a simplex.
   */
  Complex_coboundary_range coboundary_range(const Simplex& s) const {
    assert(contains(s));
    return Complex_coboundary_range(Complex_simplex_coboundary_iterator(this, s),
                                    Complex_simplex_coboundary_iterator(this, s, true));
  }

  // typedef Simplex_iterator<Skeleton_blocker_complex,Superior_link> Complex_simplex_iterator;
  typedef Simplex_iterator<Skeleton_blocker_complex> Complex_simplex_iterator;

  typedef boost::iterator_range < Complex_simplex_iterator > Complex_simplex_range;

  /**
   * @brief Returns a Complex_simplex_range over all the simplices of the complex
   */
  Complex_simplex_range complex_simplex_range() const {
    Complex_simplex_iterator end(this, true);
    if (empty()) {
      return Complex_simplex_range(end, end);
    } else {
      Complex_simplex_iterator begin(this);
      return Complex_simplex_range(begin, end);
    }
  }

  //@}

  /** @name Blockers iterators
   */
  //@{
 private:
  /**
   * @brief Iterator over the blockers adjacent to a vertex
   */
  typedef Blocker_iterator_around_vertex_internal<
  typename std::multimap<Vertex_handle, Simplex *>::iterator,
  Blocker_handle>
  Complex_blocker_around_vertex_iterator;

  /**
   * @brief Iterator over (constant) blockers adjacent to a vertex
   */
  typedef Blocker_iterator_around_vertex_internal<
  typename std::multimap<Vertex_handle, Simplex *>::const_iterator,
  const Blocker_handle>
  Const_complex_blocker_around_vertex_iterator;

  typedef boost::iterator_range <Complex_blocker_around_vertex_iterator> Complex_blocker_around_vertex_range;
  typedef boost::iterator_range <Const_complex_blocker_around_vertex_iterator>
  Const_complex_blocker_around_vertex_range;

 public:
  /**
   * @brief Returns a range of the blockers of the complex passing through a vertex
   */
  Complex_blocker_around_vertex_range blocker_range(Vertex_handle v) {
    auto begin = Complex_blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
    auto end = Complex_blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
    return Complex_blocker_around_vertex_range(begin, end);
  }

  /**
   * @brief Returns a range of the blockers of the complex passing through a vertex
   */
  Const_complex_blocker_around_vertex_range const_blocker_range(Vertex_handle v) const {
    auto begin = Const_complex_blocker_around_vertex_iterator(blocker_map_.lower_bound(v));
    auto end = Const_complex_blocker_around_vertex_iterator(blocker_map_.upper_bound(v));
    return Const_complex_blocker_around_vertex_range(begin, end);
  }

 private:
  /**
   * @brief Iterator over the blockers.
   */
  typedef Blocker_iterator_internal<
  typename std::multimap<Vertex_handle, Simplex *>::iterator,
  Blocker_handle>
  Complex_blocker_iterator;

  /**
   * @brief Iterator over the (constant) blockers.
   */
  typedef Blocker_iterator_internal<
  typename std::multimap<Vertex_handle, Simplex *>::const_iterator,
  const Blocker_handle>
  Const_complex_blocker_iterator;

  typedef boost::iterator_range <Complex_blocker_iterator> Complex_blocker_range;
  typedef boost::iterator_range <Const_complex_blocker_iterator> Const_complex_blocker_range;

 public:
  /**
   * @brief Returns a range of the blockers of the complex
   */
  Complex_blocker_range blocker_range() {
    auto begin = Complex_blocker_iterator(blocker_map_.begin(), blocker_map_.end());
    auto end = Complex_blocker_iterator(blocker_map_.end(), blocker_map_.end());
    return Complex_blocker_range(begin, end);
  }

  /**
   * @brief Returns a range of the blockers of the complex
   */
  Const_complex_blocker_range const_blocker_range() const {
    auto begin = Const_complex_blocker_iterator(blocker_map_.begin(), blocker_map_.end());
    auto end = Const_complex_blocker_iterator(blocker_map_.end(), blocker_map_.end());
    return Const_complex_blocker_range(begin, end);
  }

  //@}

  /////////////////////////////////////////////////////////////////////////////
  /** @name Print and IO methods
   */
  //@{
 public:
  std::string to_string() const {
    std::ostringstream stream;
    stream << num_vertices() << " vertices:\n" << vertices_to_string() << std::endl;
    stream << num_edges() << " edges:\n" << edges_to_string() << std::endl;
    stream << num_blockers() << " blockers:\n" << blockers_to_string() << std::endl;
    return stream.str();
  }

  std::string vertices_to_string() const {
    std::ostringstream stream;
    for (auto vertex : vertex_range()) {
      stream << "{" << (*this)[vertex].get_id() << "} ";
    }
    stream << std::endl;
    return stream.str();
  }

  std::string edges_to_string() const {
    std::ostringstream stream;
    for (auto edge : edge_range())
      stream << "{" << (*this)[edge].first() << "," << (*this)[edge].second() << "} ";
    stream << std::endl;
    return stream.str();
  }

  std::string blockers_to_string() const {
    std::ostringstream stream;

    for (auto b : const_blocker_range())
      stream << *b << std::endl;
    return stream.str();
  }
  //@}
};

/**
 * build a simplicial complex from a collection
 * of top faces.
 * return the total number of simplices
 */
template<typename Complex, typename SimplexHandleIterator>
Complex make_complex_from_top_faces(SimplexHandleIterator simplex_begin, SimplexHandleIterator simplex_end,
                                    bool is_flag_complex = false) {
  //todo use add_simplex instead! should be more efficient and more elegant :)
  typedef typename Complex::Simplex Simplex;
  std::vector<Simplex> simplices;
  for (auto top_face = simplex_begin; top_face != simplex_end; ++top_face) {
    auto subfaces_topface = subfaces(*top_face);
    simplices.insert(simplices.end(), subfaces_topface.begin(), subfaces_topface.end());
  }
  return Complex(simplices.begin(), simplices.end(), is_flag_complex);
}

}  // namespace skbl

}  // namespace Gudhi

#include "Skeleton_blocker_simplifiable_complex.h"

#endif  // SKELETON_BLOCKER_COMPLEX_H_
