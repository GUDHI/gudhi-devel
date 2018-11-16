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

#ifndef TOPLEX_MAP_H
#define TOPLEX_MAP_H

#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <limits>

namespace Gudhi {

/**
 * \brief Toplex map data structure for representing unfiltered simplicial complexes.
 *
 * \details A Toplex_map is an unordered map from vertices to maximal simplices (aka. toplices).
 *
 * \ingroup toplex_map */
class Toplex_map {
 public:
  /** Vertex is the type of vertices. */
  using Vertex = std::size_t;

  /** Simplex is the type of simplices. */
  using Simplex = std::set<Toplex_map::Vertex>;

  /** The type of the pointers to maximal simplices. */
  using Simplex_ptr = std::shared_ptr<Toplex_map::Simplex>;

  struct Sptr_hash {
    std::size_t operator()(const Toplex_map::Simplex_ptr& s) const;
  };

  struct Sptr_equal {
    std::size_t operator()(const Toplex_map::Simplex_ptr& a, const Toplex_map::Simplex_ptr& b) const;
  };

  /** The type of the sets of Toplex_map::Simplex_ptr. */
  using Simplex_ptr_set = std::unordered_set<Toplex_map::Simplex_ptr, Sptr_hash, Sptr_equal>;

  /** \brief Adds the given simplex to the complex.
   * Nothing happens if the simplex has a coface in the complex (i.e. it is a face of one of the toplices). */
  template <typename Input_vertex_range>
  void insert_simplex(const Input_vertex_range& vertex_range);

  /** \brief Removes the given simplex and its cofaces from the complex.
   * Its faces are kept inside. */
  template <typename Input_vertex_range>
  void remove_simplex(const Input_vertex_range& vertex_range);

  /** Does a simplex belong to the complex ? */
  template <typename Input_vertex_range>
  bool membership(const Input_vertex_range& vertex_range) const;

  /** Does a simplex is a toplex ? */
  template <typename Input_vertex_range>
  bool maximality(const Input_vertex_range& vertex_range) const;

  /** Gives a set of pointers to the maximal cofaces of a simplex.
   * Gives all the toplices if given the empty simplex.
   * Gives not more than max_number maximal cofaces if max_number is strictly positive. */
  template <typename Input_vertex_range>
  Toplex_map::Simplex_ptr_set maximal_cofaces(const Input_vertex_range& vertex_range,
                                              const std::size_t max_number = 0) const;

  /** Gives a set of pointers to the maximal simplices.
   * Gives not more than max_number maximal cofaces if max_number is strictly positive. */
  Toplex_map::Simplex_ptr_set maximal_simplices(const std::size_t max_number = 0) const {
    return maximal_cofaces(Simplex(), max_number);
  }

  /** Contracts one edge in the complex.
   * The edge has to verify the link condition if you want to preserve topology.
   * Returns the remaining vertex. */
  Toplex_map::Vertex contraction(const Toplex_map::Vertex x, const Toplex_map::Vertex y);

  /** Removes a vertex from any simplex containing it. */
  void remove_vertex(const Toplex_map::Vertex x);

  /** \brief Number of maximal simplices. */
  std::size_t num_maximal_simplices() const { return maximal_simplices().size(); }

  /** \brief Number of vertices. */
  std::size_t num_vertices() const { return t0.size(); }

  std::set<Toplex_map::Vertex> unitary_collapse(const Toplex_map::Vertex k, const Toplex_map::Vertex d);

  /** Adds the given simplex to the complex.
   * The simplex must not be in the complex already, and it must not contain one of the current toplices. */
  template <typename Input_vertex_range>
  void insert_independent_simplex(const Input_vertex_range& vertex_range);

 protected:
  /** \internal Gives an index in order to look for a simplex quickly. */
  template <typename Input_vertex_range>
  Toplex_map::Vertex best_index(const Input_vertex_range& vertex_range) const;

  /** \internal The map from vertices to toplices */
  std::unordered_map<Toplex_map::Vertex, Toplex_map::Simplex_ptr_set> t0;

  const Toplex_map::Vertex VERTEX_UPPER_BOUND = std::numeric_limits<Toplex_map::Vertex>::max();

  /** \internal Removes a toplex without adding facets after. */
  void erase_maximal(const Toplex_map::Simplex_ptr& sptr);
};

// Pointers are also used as key in the hash sets.
template <typename Input_vertex_range>
Toplex_map::Simplex_ptr get_key(const Input_vertex_range& vertex_range);

// Is the first simplex a face of the second ?
template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1& vertex_range1, const Input_vertex_range2& vertex_range2);

// All the facets of the given simplex.
template <typename Input_vertex_range>
std::vector<Toplex_map::Simplex> facets(const Input_vertex_range& vertex_range);

template <typename Input_vertex_range>
void Toplex_map::insert_simplex(const Input_vertex_range& vertex_range) {
  if (membership(vertex_range)) return;
  bool replace_facets = true;
  for (const Toplex_map::Simplex& facet : facets(vertex_range))
    if (!maximality(facet)) {
      replace_facets = false;
      break;
    }
  if (replace_facets)
    for (const Toplex_map::Simplex& facet : facets(vertex_range)) erase_maximal(get_key(facet));
  else
    for (const Toplex_map::Vertex& v : vertex_range)
      if (t0.count(v))
        for (const Toplex_map::Simplex_ptr& fptr : Simplex_ptr_set(t0.at(v)))
          // Copy constructor needed because the set is modified
          if (included(*fptr, vertex_range)) erase_maximal(fptr);
  // We erase all the maximal faces of the simplex
  insert_independent_simplex(vertex_range);
}

template <typename Input_vertex_range>
void Toplex_map::remove_simplex(const Input_vertex_range& vertex_range) {
  if (vertex_range.begin() == vertex_range.end()) t0.clear();
  // Removal of the empty simplex means cleaning everything
  else {
    const Toplex_map::Vertex& v = best_index(vertex_range);
    if (t0.count(v))
      for (const Toplex_map::Simplex_ptr& sptr : Simplex_ptr_set(t0.at(v)))
        // Copy constructor needed because the set is modified
        if (included(vertex_range, *sptr)) {
          erase_maximal(sptr);
          for (const Toplex_map::Simplex& f : facets(vertex_range))
            if (!membership(f)) insert_independent_simplex(f);
          // We add the facets which are new maximal simplices
        }
  }
}

template <typename Input_vertex_range>
bool Toplex_map::membership(const Input_vertex_range& vertex_range) const {
  if (t0.size() == 0) return false;
  const Toplex_map::Vertex& v = best_index(vertex_range);
  if (!t0.count(v)) return false;
  if (maximality(vertex_range)) return true;
  for (const Toplex_map::Simplex_ptr& sptr : t0.at(v))
    if (included(vertex_range, *sptr)) return true;
  return false;
}

template <typename Input_vertex_range>
bool Toplex_map::maximality(const Input_vertex_range& vertex_range) const {
  const Toplex_map::Vertex& v = best_index(vertex_range);
  if (!t0.count(v)) return false;
  return t0.at(v).count(get_key(vertex_range));
}

template <typename Input_vertex_range>
Toplex_map::Simplex_ptr_set Toplex_map::maximal_cofaces(const Input_vertex_range& vertex_range,
                                                        const std::size_t max_number) const {
  Simplex_ptr_set cofaces;
  if (maximality(vertex_range))
    cofaces.emplace(get_key(vertex_range));
  else if (vertex_range.begin() == vertex_range.end())
    for (const auto& kv : t0)
      for (const Toplex_map::Simplex_ptr& sptr : kv.second) {
        // kv.second is a Simplex_ptr_set
        cofaces.emplace(sptr);
        if (cofaces.size() == max_number) return cofaces;
      }
  else {
    const Toplex_map::Vertex& v = best_index(vertex_range);
    if (t0.count(v))
      for (const Toplex_map::Simplex_ptr& sptr : t0.at(v))
        if (included(vertex_range, *sptr)) {
          cofaces.emplace(sptr);
          if (cofaces.size() == max_number) return cofaces;
        }
  }
  return cofaces;
}

Toplex_map::Vertex Toplex_map::contraction(const Toplex_map::Vertex x, const Toplex_map::Vertex y) {
  if (!t0.count(x)) return y;
  if (!t0.count(y)) return x;
  int k, d;
  if (t0.at(x).size() > t0.at(y).size())
    k = x, d = y;
  else
    k = y, d = x;
  for (const Toplex_map::Simplex_ptr& sptr : Simplex_ptr_set(t0.at(d))) {
    // Copy constructor needed because the set is modified
    Simplex sigma(*sptr);
    erase_maximal(sptr);
    sigma.erase(d);
    sigma.insert(k);
    insert_simplex(sigma);
  }
  return k;
}

std::set<Toplex_map::Vertex> Toplex_map::unitary_collapse(const Toplex_map::Vertex k, const Toplex_map::Vertex d) {
  std::set<Toplex_map::Vertex> r;
  for (const Toplex_map::Simplex_ptr& sptr : Simplex_ptr_set(t0.at(d))) {
    // Copy constructor needed because the set is modified
    Simplex sigma(*sptr);
    erase_maximal(sptr);
    sigma.erase(d);
    for (const Toplex_map::Vertex v : sigma) r.insert(v);
    sigma.insert(k);
    insert_simplex(sigma);
  }
  return r;
}

template <typename Input_vertex_range>
void Toplex_map::insert_independent_simplex(const Input_vertex_range& vertex_range) {
  auto key = get_key(vertex_range);
  for (const Toplex_map::Vertex& v : vertex_range) {
    if (!t0.count(v)) t0.emplace(v, Simplex_ptr_set());
    t0.at(v).emplace(key);
  }
}

void Toplex_map::remove_vertex(const Toplex_map::Vertex x) {
  for (const Toplex_map::Simplex_ptr& sptr : Simplex_ptr_set(t0.at(x))) {
    Simplex sigma(*sptr);
    erase_maximal(sptr);
    sigma.erase(x);
    insert_simplex(sigma);
  }
}

inline void Toplex_map::erase_maximal(const Toplex_map::Simplex_ptr& sptr) {
  Simplex sigma(*sptr);
  if (sptr->size() == 0) sigma.insert(VERTEX_UPPER_BOUND);
  for (const Toplex_map::Vertex& v : sigma) {
    t0.at(v).erase(sptr);
    if (t0.at(v).size() == 0) t0.erase(v);
  }
}

template <typename Input_vertex_range>
Toplex_map::Vertex Toplex_map::best_index(const Input_vertex_range& vertex_range) const {
  std::size_t min = std::numeric_limits<size_t>::max();
  Vertex arg_min = VERTEX_UPPER_BOUND;
  for (const Toplex_map::Vertex& v : vertex_range)
    if (!t0.count(v))
      return v;
    else if (t0.at(v).size() < min)
      min = t0.at(v).size(), arg_min = v;
  return arg_min;
}

std::size_t Toplex_map::Sptr_equal::operator()(const Toplex_map::Simplex_ptr& s1,
                                               const Toplex_map::Simplex_ptr& s2) const {
  if (s1->size() != s2->size()) return false;
  return included(*s1, *s2);
  // inclusion tests equality for same size simplices
}

std::size_t Toplex_map::Sptr_hash::operator()(const Toplex_map::Simplex_ptr& s) const {
  std::hash<double> h_f;
  // double hash works better than int hash
  size_t h = 0;
  for (const Toplex_map::Vertex& v : *s) h += h_f(static_cast<double>(v));
  return h;
}

template <typename Input_vertex_range>
Toplex_map::Simplex_ptr get_key(const Input_vertex_range& vertex_range) {
  Toplex_map::Simplex s(vertex_range.begin(), vertex_range.end());
  return std::make_shared<Toplex_map::Simplex>(s);
}

template <typename Input_vertex_range1, typename Input_vertex_range2>
bool included(const Input_vertex_range1& vertex_range1, const Input_vertex_range2& vertex_range2) {
  Toplex_map::Simplex s2(vertex_range2.begin(), vertex_range2.end());
  for (const Toplex_map::Vertex& v : vertex_range1)
    if (!s2.count(v)) return false;
  return true;
}

template <typename Input_vertex_range>
std::vector<Toplex_map::Simplex> facets(const Input_vertex_range& vertex_range) {
  std::vector<Toplex_map::Simplex> facets;
  Toplex_map::Simplex f(vertex_range.begin(), vertex_range.end());
  for (const Toplex_map::Vertex& v : vertex_range) {
    f.erase(v);
    facets.emplace_back(f);
    f.insert(v);
  }
  return facets;
}

}  // namespace Gudhi

#endif /* TOPLEX_MAP_H */
