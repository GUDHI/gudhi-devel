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
#ifndef SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SIMPLIFIABLE_COMPLEX_H_
#define SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SIMPLIFIABLE_COMPLEX_H_

#include <list>
#include <vector>
#include <set>

#include "gudhi/Skeleton_blocker/Skeleton_blocker_sub_complex.h"

namespace Gudhi {

namespace skbl {

/**
 * Returns true iff the blocker 'sigma' is popable.
 * To define popable, let us call 'L' the complex that
 * consists in the current complex without the blocker 'sigma'.
 * A blocker 'sigma' is then "popable" if the link of 'sigma'
 * in L is reducible.
 *
 */
template<typename SkeletonBlockerDS>
bool Skeleton_blocker_complex<SkeletonBlockerDS>::is_popable_blocker(Blocker_handle sigma) const {
  assert(this->contains_blocker(*sigma));
  Skeleton_blocker_link_complex<Skeleton_blocker_complex> link_blocker_sigma;
  build_link_of_blocker(*this, *sigma, link_blocker_sigma);

  bool res = link_blocker_sigma.is_contractible() == CONTRACTIBLE;
  return res;
}

/**
 * Removes all the popable blockers of the complex and delete them.
 * @returns the number of popable blockers deleted
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_popable_blockers() {
  std::list<Vertex_handle> vertex_to_check;
  for (auto v : this->vertex_range())
    vertex_to_check.push_front(v);

  while (!vertex_to_check.empty()) {
    Vertex_handle v = vertex_to_check.front();
    vertex_to_check.pop_front();

    bool blocker_popable_found = true;
    while (blocker_popable_found) {
      blocker_popable_found = false;
      for (auto block : this->blocker_range(v)) {
        if (this->is_popable_blocker(block)) {
          for (Vertex_handle nv : *block)
            if (nv != v) vertex_to_check.push_back(nv);
          this->delete_blocker(block);
          blocker_popable_found = true;
          break;
        }
      }
    }
  }
}

/**
 * Removes all the popable blockers of the complex passing through v and delete them.
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_popable_blockers(Vertex_handle v) {
  bool blocker_popable_found = true;
  while (blocker_popable_found) {
    blocker_popable_found = false;
    for (auto block : this->blocker_range(v)) {
      if (is_popable_blocker(block)) {
        this->delete_blocker(block);
        blocker_popable_found = true;
      }
    }
  }
}

/**
 * @brief Removes all the popable blockers of the complex passing through v and delete them.
 * Also remove popable blockers in the neighborhood if they became popable.
 *
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_all_popable_blockers(Vertex_handle v) {
  std::list<Vertex_handle> vertex_to_check;
  vertex_to_check.push_front(v);

  while (!vertex_to_check.empty()) {
    Vertex_handle v = vertex_to_check.front();
    vertex_to_check.pop_front();

    bool blocker_popable_found = true;
    while (blocker_popable_found) {
      blocker_popable_found = false;
      for (auto block : this->blocker_range(v)) {
        if (this->is_popable_blocker(block)) {
          for (Vertex_handle nv : *block)
            if (nv != v) vertex_to_check.push_back(nv);
          this->delete_blocker(block);
          blocker_popable_found = true;
          break;
        }
      }
    }
  }
}

/**
 * Remove the star of the vertice 'v'
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_star(Vertex_handle v) {
  // we remove the blockers that are not consistent anymore
  update_blockers_after_remove_star_of_vertex_or_edge(Simplex_handle(v));
  while (this->degree(v) > 0) {
    Vertex_handle w(* (adjacent_vertices(v.vertex, this->skeleton).first));
    this->remove_edge(v, w);
  }
  this->remove_vertex(v);
}

/**
 * after removing the star of a simplex, blockers sigma that contains this simplex must be removed.
 * Furthermore, all simplices tau of the form sigma \setminus simplex_to_be_removed must be added
 * whenever the dimension of tau is at least 2.
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::update_blockers_after_remove_star_of_vertex_or_edge(
    const Simplex_handle& simplex_to_be_removed) {
  std::list <Blocker_handle> blockers_to_update;
  if (simplex_to_be_removed.empty()) return;

  auto v0 = simplex_to_be_removed.first_vertex();
  for (auto blocker : this->blocker_range(v0)) {
    if (blocker->contains(simplex_to_be_removed))
      blockers_to_update.push_back(blocker);
  }

  for (auto blocker_to_update : blockers_to_update) {
    Simplex_handle sub_blocker_to_be_added;
    bool sub_blocker_need_to_be_added =
        (blocker_to_update->dimension() - simplex_to_be_removed.dimension()) >= 2;
    if (sub_blocker_need_to_be_added) {
      sub_blocker_to_be_added = *blocker_to_update;
      sub_blocker_to_be_added.difference(simplex_to_be_removed);
    }
    this->delete_blocker(blocker_to_update);
    if (sub_blocker_need_to_be_added)
      this->add_blocker(sub_blocker_to_be_added);
  }
}

/**
 * Remove the star of the edge connecting vertices a and b.
 * @returns the number of blocker that have been removed
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_star(Vertex_handle a, Vertex_handle b) {
  update_blockers_after_remove_star_of_vertex_or_edge(Simplex_handle(a, b));
  // we remove the edge
  this->remove_edge(a, b);
}

/**
 * Remove the star of the edge 'e'.
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_star(Edge_handle e) {
  return remove_star(this->first_vertex(e), this->second_vertex(e));
}

/**
 * Remove the star of the simplex 'sigma' which needs to belong to the complex
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_star(const Simplex_handle& sigma) {
  assert(this->contains(sigma));
  if (sigma.dimension() == 0) {
    remove_star(sigma.first_vertex());
  } else if (sigma.dimension() == 1) {
    remove_star(sigma.first_vertex(), sigma.last_vertex());
  } else {
    remove_blocker_containing_simplex(sigma);
    this->add_blocker(sigma);
  }
}

/**
 * @brief add a maximal simplex plus all its cofaces. All vertices lower than the higher vertex of
 * sigma must already be present.
 * @details the simplex must have dimension greater than one (otherwise use add_vertex or add_edge).
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::add_simplex(const Simplex_handle& sigma) {
  assert(!this->contains(sigma));
  assert(sigma.dimension() > 1);

  int num_vertex_to_add = 0;
  for (auto v : sigma)
    if (!contains_vertex(v)) ++num_vertex_to_add;
  while (num_vertex_to_add--) add_vertex();

  for (auto u_it = sigma.begin(); u_it != sigma.end(); ++u_it)
    for (auto v_it = u_it; ++v_it != sigma.end(); /**/) {
      std::cout << "add edge" << *u_it << " " << *v_it << std::endl;
      add_edge(*u_it, *v_it);
    }
  remove_blocker_include_in_simplex(sigma);
}

/**
 * remove all blockers that contains sigma
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_blocker_containing_simplex(const Simplex_handle& sigma) {
  std::vector <Blocker_handle> blockers_to_remove;
  for (auto blocker : this->blocker_range(sigma.first_vertex())) {
    if (blocker->contains(sigma))
      blockers_to_remove.push_back(blocker);
  }
  for (auto blocker_to_update : blockers_to_remove)
    this->delete_blocker(blocker_to_update);
}

/**
 * remove all blockers that contains sigma
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::remove_blocker_include_in_simplex(const Simplex_handle& sigma) {
  std::vector <Blocker_handle> blockers_to_remove;
  for (auto blocker : this->blocker_range(sigma.first_vertex())) {
    if (sigma.contains(*blocker))
      blockers_to_remove.push_back(blocker);
  }
  for (auto blocker_to_update : blockers_to_remove)
    this->delete_blocker(blocker_to_update);
}

/**
 * Compute simplices beta such that a.beta is an order 0 blocker
 * that may be used to construct a new blocker after contracting ab.
 * It requires that the link condition is satisfied.
 */
template<typename SkeletonBlockerDS>
void Skeleton_blocker_complex<SkeletonBlockerDS>::tip_blockers(Vertex_handle a, Vertex_handle b,
                                                               std::vector<Simplex_handle> & buffer) const {
  for (auto const & blocker : this->const_blocker_range(a)) {
    Simplex_handle beta = (*blocker);
    beta.remove_vertex(a);
    buffer.push_back(beta);
  }

  Simplex_handle n;
  this->add_neighbours(b, n);
  this->remove_neighbours(a, n);
  n.remove_vertex(a);


  for (Vertex_handle y : n) {
    Simplex_handle beta;
    beta.add_vertex(y);
    buffer.push_back(beta);
  }
}

/**
 * @brief "Replace" the edge 'bx' by the edge 'ax'.
 * Assume that the edge 'bx' was present whereas 'ax' was not.
 * Precisely, it does not replace edges, but remove 'bx' and then add 'ax'.
 * The visitor 'on_swaped_edge' is called just after edge 'ax' had been added
 * and just before edge 'bx' had been removed. That way, it can
 * eventually access to information of 'ax'.
 */
template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::swap_edge(Vertex_handle a, Vertex_handle b, Vertex_handle x) {
  this->add_edge(a, x);
  if (this->visitor) this->visitor->on_swaped_edge(a, b, x);
  this->remove_edge(b, x);
}

template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::delete_blockers_around_vertex(Vertex_handle v) {
  std::list <Blocker_handle> blockers_to_delete;
  for (auto blocker : this->blocker_range(v)) {
    blockers_to_delete.push_back(blocker);
  }
  while (!blockers_to_delete.empty()) {
    this->remove_blocker(blockers_to_delete.back());
    blockers_to_delete.pop_back();
  }
}

/**
 * @brief removes all blockers passing through the edge 'ab'
 */
template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::delete_blockers_around_edge(Vertex_handle a, Vertex_handle b) {
  std::list<Blocker_handle> blocker_to_delete;
  for (auto blocker : this->blocker_range(a))
    if (blocker->contains(b)) blocker_to_delete.push_back(blocker);
  while (!blocker_to_delete.empty()) {
    this->delete_blocker(blocker_to_delete.back());
    blocker_to_delete.pop_back();
  }
}

/**
 * Contracts the edge connecting vertices a and b.
 * @remark If the link condition Link(ab) = Link(a) inter Link(b) is not satisfied,
 * it removes first all blockers passing through 'ab'
 */
template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::contract_edge(Vertex_handle a, Vertex_handle b) {
  assert(this->contains_vertex(a));
  assert(this->contains_vertex(b));

  if(this->contains_edge(a, b))
    this->add_edge(a, b);

  // if some blockers passes through 'ab', we need to remove them.
  if (!link_condition(a, b))
    delete_blockers_around_edge(a, b);

  std::set<Simplex_handle> blockers_to_add;

  get_blockers_to_be_added_after_contraction(a, b, blockers_to_add);

  delete_blockers_around_vertices(a, b);

  update_edges_after_contraction(a, b);

  this->remove_vertex(b);

  notify_changed_edges(a);

  for (auto block : blockers_to_add)
    this->add_blocker(block);

  assert(this->contains_vertex(a));
  assert(!this->contains_vertex(b));
}

template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::get_blockers_to_be_added_after_contraction(Vertex_handle a, Vertex_handle b, 
    std::set<Simplex_handle>& blockers_to_add) {
  blockers_to_add.clear();

  typedef Skeleton_blocker_link_complex<Skeleton_blocker_complex<SkeletonBlockerDS> > LinkComplexType;

  LinkComplexType link_a(*this, a);
  LinkComplexType link_b(*this, b);

  std::vector<Simplex_handle> vector_alpha, vector_beta;

  tip_blockers(a, b, vector_alpha);
  tip_blockers(b, a, vector_beta);

  for (auto alpha = vector_alpha.begin(); alpha != vector_alpha.end(); ++alpha) {
    for (auto beta = vector_beta.begin(); beta != vector_beta.end(); ++beta) {
      Simplex_handle sigma = *alpha;
      sigma.union_vertices(*beta);
      Root_simplex_handle sigma_id = this->get_id(sigma);
      if (this->contains(sigma) &&
          proper_faces_in_union<Skeleton_blocker_complex < SkeletonBlockerDS >> (sigma_id, link_a, link_b)) {
        // Blocker_handle blocker = new Simplex_handle(sigma);
        sigma.add_vertex(a);
        blockers_to_add.insert(sigma);
      }
    }
  }
}

/**
 * delete all blockers that passes through a or b
 */
template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::delete_blockers_around_vertices(Vertex_handle a, Vertex_handle b) {
  std::vector<Blocker_handle> blocker_to_delete;
  for (auto bl : this->blocker_range(a))
    blocker_to_delete.push_back(bl);
  for (auto bl : this->blocker_range(b))
    blocker_to_delete.push_back(bl);
  while (!blocker_to_delete.empty()) {
    this->delete_blocker(blocker_to_delete.back());
    blocker_to_delete.pop_back();
  }
}

template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::update_edges_after_contraction(Vertex_handle a, Vertex_handle b) {
  // We update the set of edges
  this->remove_edge(a, b);

  // For all edges {b,x} incident to b,
  // we remove {b,x} and add {a,x} if not already there.
  while (this->degree(b) > 0) {
    Vertex_handle x(*(adjacent_vertices(b.vertex, this->skeleton).first));
    if (!this->contains_edge(a, x))
      // we 'replace' the edge 'bx' by the edge 'ax'
      this->swap_edge(a, b, x);
    else
      this->remove_edge(b, x);
  }
}

template<typename SkeletonBlockerDS>
void
Skeleton_blocker_complex<SkeletonBlockerDS>::notify_changed_edges(Vertex_handle a) {
  // We notify the visitor that all edges incident to 'a' had changed
  boost_adjacency_iterator v, v_end;

  for (tie(v, v_end) = adjacent_vertices(a.vertex, this->skeleton); v != v_end; ++v)
    if (this->visitor) this->visitor->on_changed_edge(a, Vertex_handle(*v));
}


}  // namespace skbl

}  // namespace Gudhi

#endif  // SRC_SKELETON_BLOCKER_INCLUDE_GUDHI_SKELETON_BLOCKER_SIMPLIFIABLE_COMPLEX_H_
