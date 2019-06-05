/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       François Godi, Vincent Rouvreau
 *
 *    Copyright (C) 2018  INRIA
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

#ifndef LAZY_TOPLEX_MAP_H
#define LAZY_TOPLEX_MAP_H

#include <gudhi/Toplex_map.h>
#include <boost/heap/fibonacci_heap.hpp>

namespace Gudhi {

/**
 * \brief Lazy toplex map data structure for representing unfiltered simplicial complexes.
 *
 * \details A Toplex_map is an unordered map from vertices to maximal simplices (aka. toplices).
 * The lazy version is not always up to date as it requires clean operation in order to be.
 *
 * \ingroup toplex_map */
class Lazy_toplex_map {
 public:
  /** Vertex is the type of vertices. */
  using Vertex = Toplex_map::Vertex;

  /** Simplex is the type of simplices. */
  using Simplex = Toplex_map::Simplex;

  /** The type of the pointers to maximal simplices. */
  using Simplex_ptr = Toplex_map::Simplex_ptr;

  /** The type of the sets of Simplex_ptr. */
  using Simplex_ptr_set = Toplex_map::Simplex_ptr_set;

  /** Adds the given simplex to the complex.
   * The simplex must not be in the complex already, and it must not contain one of the current toplices. */
  template <typename Input_vertex_range>
  void insert_independent_simplex(const Input_vertex_range &vertex_range);

  /** \brief Adds the given simplex to the complex.
   * Nothing happens if the simplex is already in the complex (i.e. it is a face of one of the toplices). */
  template <typename Input_vertex_range>
  bool insert_simplex(const Input_vertex_range &vertex_range);

  /** \brief Removes the given simplex and its cofaces from the complex.
   * Its faces are kept inside. */
  template <typename Input_vertex_range>
  void remove_simplex(const Input_vertex_range &vertex_range);

  /** Does a simplex belong to the complex ? */
  template <typename Input_vertex_range>
  bool membership(const Input_vertex_range &vertex_range);

  /** Do all the facets of a simplex belong to the complex ? */
  template <typename Input_vertex_range>
  bool all_facets_inside(const Input_vertex_range &vertex_range);

  /** Contracts one edge in the complex.
   * The edge has to verify the link condition if you want to preserve topology.
   * Returns the remaining vertex. */
  Vertex contraction(const Vertex x, const Vertex y);

  /** \brief Number of maximal simplices. */
  std::size_t num_maximal_simplices() const { return size; }

  /** \brief Number of vertices. */
  std::size_t num_vertices() const { return t0.size(); }

 private:
  template <typename Input_vertex_range>
  void erase_max(const Input_vertex_range &vertex_range);
  template <typename Input_vertex_range>
  Vertex best_index(const Input_vertex_range &vertex_range);
  void clean(const Vertex v);

  std::unordered_map<Vertex, std::size_t> gamma0_lbounds;

  std::unordered_map<Vertex, Simplex_ptr_set> t0;
  bool empty_toplex = true;  // Is the empty simplex a toplex ?

  typedef boost::heap::fibonacci_heap<std::pair<std::size_t, Vertex>> PriorityQueue;
  PriorityQueue cleaning_priority;
  std::unordered_map<Vertex, PriorityQueue::handle_type> cp_handles;

  std::size_t get_gamma0_lbound(const Vertex v) const;

  std::size_t size_lbound = 0;
  std::size_t size = 0;

  const double ALPHA = 4;  // time
  const double BETTA = 8;  // memory
};

template <typename Input_vertex_range>
void Lazy_toplex_map::insert_independent_simplex(const Input_vertex_range &vertex_range) {
  for (const Vertex &v : vertex_range)
    if (!gamma0_lbounds.count(v))
      gamma0_lbounds.emplace(v, 1);
    else
      gamma0_lbounds[v]++;
  size_lbound++;
  insert_simplex(vertex_range);
}

template <typename Input_vertex_range>
bool Lazy_toplex_map::insert_simplex(const Input_vertex_range &vertex_range) {
  Simplex sigma(vertex_range.begin(), vertex_range.end());
  // Check empty face management
  empty_toplex = (sigma.size() == 0);
  Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
  bool inserted = false;
  for (const Vertex &v : sigma) {
    if (!t0.count(v)) {
      t0.emplace(v, Simplex_ptr_set());
      auto v_handle = cleaning_priority.push(std::make_pair(0, v));
      cp_handles.emplace(v, v_handle);
    }
    inserted = t0.at(v).emplace(sptr).second;
    cleaning_priority.update(cp_handles.at(v), std::make_pair(t0.at(v).size() - get_gamma0_lbound(v), v));
  }
  if (inserted) size++;
  if (size > (size_lbound + 1) * BETTA) clean(cleaning_priority.top().second);
  return inserted;
}

template <typename Input_vertex_range>
void Lazy_toplex_map::remove_simplex(const Input_vertex_range &vertex_range) {
  if (vertex_range.begin() == vertex_range.end()) {
    t0.clear();
    gamma0_lbounds.clear();
    cleaning_priority.clear();
    size_lbound = 0;
    size = 0;
    empty_toplex = false;
  } else {
    const Vertex &v = best_index(vertex_range);
    // Copy constructor needed because the set is modified
    if (t0.count(v))
      for (const Simplex_ptr &sptr : Simplex_ptr_set(t0.at(v)))
        if (included(vertex_range, *sptr)) {
          erase_max(*sptr);
          for (const Simplex &f : facets(vertex_range)) insert_independent_simplex(f);
        }
  }
}

template <typename Input_vertex_range>
bool Lazy_toplex_map::membership(const Input_vertex_range &vertex_range) {
  if (t0.size() == 0 && !empty_toplex) return false;            // empty complex
  if (vertex_range.begin() == vertex_range.end()) return true;  // empty query simplex
  Vertex v = best_index(vertex_range);
  if (!t0.count(v)) return false;
  for (const Simplex_ptr &sptr : t0.at(v))
    if (included(vertex_range, *sptr)) return true;
  return false;
}

template <typename Input_vertex_range>
bool Lazy_toplex_map::all_facets_inside(const Input_vertex_range &vertex_range) {
  Simplex sigma(vertex_range.begin(), vertex_range.end());
  Vertex v = best_index(sigma);
  if (!t0.count(v)) return false;
  Simplex f = sigma;
  f.erase(v);
  if (!membership(f)) return false;
  std::unordered_set<Vertex> facets_inside;
  for (const Simplex_ptr &sptr : t0.at(v))
    for (const Vertex &w : sigma) {
      f = sigma;
      f.erase(w);
      if (included(f, *sptr)) facets_inside.insert(w);
    }
  return facets_inside.size() == sigma.size() - 1;
}

/* Returns the remaining vertex */
Lazy_toplex_map::Vertex Lazy_toplex_map::contraction(const Vertex x, const Vertex y) {
  if (!t0.count(x)) return y;
  if (!t0.count(y)) return x;
  Vertex k, d;
  if (t0.at(x).size() > t0.at(y).size())
    k = x, d = y;
  else
    k = y, d = x;
  // Copy constructor needed because the set is modified
  for (const Simplex_ptr &sptr : Simplex_ptr_set(t0.at(d))) {
    Simplex sigma(*sptr);
    erase_max(sigma);
    sigma.erase(d);
    sigma.insert(k);
    insert_simplex(sigma);
  }
  t0.erase(d);
  return k;
}

/* No facets insert_simplexed */
template <typename Input_vertex_range>
inline void Lazy_toplex_map::erase_max(const Input_vertex_range &vertex_range) {
  Simplex sigma(vertex_range.begin(), vertex_range.end());
  empty_toplex = false;
  Simplex_ptr sptr = std::make_shared<Simplex>(sigma);
  bool erased = false;
  for (const Vertex &v : sigma) {
    erased = t0.at(v).erase(sptr) > 0;
    if (t0.at(v).size() == 0) t0.erase(v);
  }
  if (erased) size--;
}

template <typename Input_vertex_range>
Lazy_toplex_map::Vertex Lazy_toplex_map::best_index(const Input_vertex_range &vertex_range) {
  Simplex tau(vertex_range.begin(), vertex_range.end());
  std::size_t min = std::numeric_limits<size_t>::max();
  Vertex arg_min = -1;
  for (const Vertex &v : tau)
    if (!t0.count(v))
      return v;
    else if (t0.at(v).size() < min)
      min = t0.at(v).size(), arg_min = v;
  if (min > ALPHA * get_gamma0_lbound(arg_min)) clean(arg_min);
  return arg_min;
}

std::size_t Lazy_toplex_map::get_gamma0_lbound(const Vertex v) const {
  return gamma0_lbounds.count(v) ? gamma0_lbounds.at(v) : 0;
}

void Lazy_toplex_map::clean(const Vertex v) {
  Toplex_map toplices;
  std::unordered_map<int, std::vector<Simplex>> dsorted_simplices;
  std::size_t max_dim = 0;
  for (const Simplex_ptr &sptr : Simplex_ptr_set(t0.at(v))) {
    if (sptr->size() > max_dim) {
      for (std::size_t d = max_dim + 1; d <= sptr->size(); d++) dsorted_simplices.emplace(d, std::vector<Simplex>());
      max_dim = sptr->size();
    }
    dsorted_simplices[sptr->size()].emplace_back(*sptr);
    erase_max(*sptr);
  }
  for (std::size_t d = max_dim; d >= 1; d--)
    for (const Simplex &s : dsorted_simplices.at(d))
      if (!toplices.membership(s)) toplices.insert_independent_simplex(s);
  Simplex sv;
  sv.insert(v);
  auto clean_cofaces = toplices.maximal_cofaces(sv);
  size_lbound = size_lbound - get_gamma0_lbound(v) + clean_cofaces.size();
  gamma0_lbounds[v] = clean_cofaces.size();
  for (const Simplex_ptr &sptr : clean_cofaces) insert_simplex(*sptr);
}

}  // namespace Gudhi

#endif /* LAZY_TOPLEX_MAP_H */
