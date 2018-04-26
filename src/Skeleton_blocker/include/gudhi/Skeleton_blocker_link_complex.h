/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
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

#ifndef SKELETON_BLOCKER_LINK_COMPLEX_H_
#define SKELETON_BLOCKER_LINK_COMPLEX_H_

#include <gudhi/Skeleton_blocker_complex.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {

namespace skeleton_blocker {

template<class ComplexType> class Skeleton_blocker_sub_complex;

/**
 *  \brief Class representing the link of a simplicial complex encoded by a skeleton/blockers pair.
 *  It inherits from Skeleton_blocker_sub_complex because such complex is a sub complex of a
 *  root complex.
 * \ingroup skbl
 */
template<typename ComplexType>
class Skeleton_blocker_link_complex : public Skeleton_blocker_sub_complex<
ComplexType> {
  template<typename T> friend class Skeleton_blocker_link_superior;
  typedef typename ComplexType::Edge_handle Edge_handle;

  typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

 private:
  bool only_superior_vertices_;

 public:
  typedef typename ComplexType::Vertex_handle Vertex_handle;
  typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;

  typedef typename ComplexType::Simplex Simplex;
  typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;

  typedef typename ComplexType::Blocker_handle Blocker_handle;

  typedef typename ComplexType::Root_simplex_handle::Simplex_vertex_const_iterator Root_simplex_handle_iterator;

  explicit Skeleton_blocker_link_complex(bool only_superior_vertices = false)
      : only_superior_vertices_(only_superior_vertices) { }

  /**
   * If the parameter only_superior_vertices is true,
   * only vertices greater than the one of alpha are added.
   * Only vertices are computed if only_vertices is true.
   */
  Skeleton_blocker_link_complex(const ComplexType & parent_complex,
                                const Simplex& alpha_parent_adress,
                                bool only_superior_vertices = false,
                                bool only_vertices = false)
      : only_superior_vertices_(only_superior_vertices) {
    if (!alpha_parent_adress.empty())
      build_link(parent_complex, alpha_parent_adress, only_vertices);
  }

  /**
   * If the parameter only_superior_vertices is true,
   * only vertices greater than the one of the vertex are added.
   */
  Skeleton_blocker_link_complex(const ComplexType & parent_complex,
                                Vertex_handle a_parent_adress,
                                bool only_superior_vertices = false)
      : only_superior_vertices_(only_superior_vertices) {
    Simplex alpha_simplex(a_parent_adress);
    build_link(parent_complex, alpha_simplex);
  }

  /**
   * If the parameter only_superior_vertices is true,
   * only vertices greater than the one of the edge are added.
   */
  Skeleton_blocker_link_complex(const ComplexType & parent_complex,
                                Edge_handle edge, bool only_superior_vertices =
                                false)
      : only_superior_vertices_(only_superior_vertices) {
    Simplex alpha_simplex(parent_complex.first_vertex(edge),
                          parent_complex.second_vertex(edge));
    build_link(parent_complex, alpha_simplex);
  }

 protected:
  /**
   * @brief compute vertices of the link.
   * If the boolean only_superior_vertices is true, then only the vertices
   * are greater than  vertices of alpha_parent_adress are added.
   */
  void compute_link_vertices(const ComplexType & parent_complex,
                             const Simplex& alpha_parent_adress,
                             bool only_superior_vertices,
                             bool is_alpha_blocker = false) {
    if (alpha_parent_adress.dimension() == 0) {
      // for a vertex we know exactly the number of vertices of the link (and the size of the corresponding vector)
      // thus we call a specific function that will reserve a vector with appropriate size
      this->compute_link_vertices(parent_complex,
                                  alpha_parent_adress.first_vertex(),
                                  only_superior_vertices_);
    } else {
      // we compute the intersection of neighbors of alpha and store it in link_vertices
      Simplex link_vertices_parent;
      parent_complex.add_neighbours(alpha_parent_adress, link_vertices_parent,
                                    only_superior_vertices);
      // For all vertex 'v' in this intersection, we go through all its adjacent blockers.
      // If one blocker minus 'v' is included in alpha then the vertex is not in the link complex.
      for (auto v_parent : link_vertices_parent) {
        bool new_vertex = true;
        for (auto beta : parent_complex.const_blocker_range(v_parent)) {
          if (!is_alpha_blocker || *beta != alpha_parent_adress) {
            new_vertex = !(alpha_parent_adress.contains_difference(*beta,
                                                                   v_parent));
            if (!new_vertex)
              break;
          }
        }
        if (new_vertex)
          this->add_vertex(parent_complex.get_id(v_parent));
      }
    }
  }

  /**
   * @brief compute vertices of the link.
   * If the boolean only_superior_vertices is true, then only the vertices
   * are greater than  vertices of alpha_parent_adress are added.
   */
  void compute_link_vertices(const ComplexType & parent_complex,
                             Vertex_handle alpha_parent_adress,
                             bool only_superior_vertices) {
    // for a vertex we know exactly the number of vertices of the link (and the size of the corresponding vector
    this->skeleton.m_vertices.reserve(
                                      parent_complex.degree(alpha_parent_adress));

    // For all vertex 'v' in this intersection, we go through all its adjacent blockers.
    // If one blocker minus 'v' is included in alpha then the vertex is not in the link complex.
    for (auto v_parent : parent_complex.vertex_range(alpha_parent_adress)) {
      if (!only_superior_vertices
          || v_parent.vertex > alpha_parent_adress.vertex)
        this->add_vertex(parent_complex.get_id(v_parent));
    }
  }

  void compute_link_edges(const ComplexType & parent_complex,
                          const Simplex& alpha_parent_adress,
                          bool is_alpha_blocker = false) {
    if (this->num_vertices() <= 1)
      return;

    for (auto x_link = this->vertex_range().begin();
         x_link != this->vertex_range().end(); ++x_link) {
      for (auto y_link = x_link; ++y_link != this->vertex_range().end();) {
        Vertex_handle x_parent = *parent_complex.get_address(
                                                             this->get_id(*x_link));
        Vertex_handle y_parent = *parent_complex.get_address(
                                                             this->get_id(*y_link));
        if (parent_complex.contains_edge(x_parent, y_parent)) {
          // we check that there is no blocker subset of alpha passing trough x and y
          bool new_edge = true;
          for (auto blocker_parent : parent_complex.const_blocker_range(
                                                                        x_parent)) {
            if (!is_alpha_blocker || *blocker_parent != alpha_parent_adress) {
              if (blocker_parent->contains(y_parent)) {
                new_edge = !(alpha_parent_adress.contains_difference(
                                                                     *blocker_parent, x_parent, y_parent));
                if (!new_edge)
                  break;
              }
            }
          }
          if (new_edge)
            this->add_edge_without_blockers(*x_link, *y_link);
        }
      }
    }
  }

  /**
   * @brief : Given an address in the current complex, it returns the
   * corresponding address in 'other_complex'.
   * It assumes that other_complex have a vertex 'this.get_id(address)'
   */
  boost::optional<Vertex_handle> give_equivalent_vertex(const ComplexType & other_complex,
                                                        Vertex_handle address) const {
    Root_vertex_handle id((*this)[address].get_id());
    return other_complex.get_address(id);
  }

  /*
   * compute the blockers of the link if is_alpha_blocker is false.
   * Otherwise, alpha is a blocker, and the link is computed in the complex where
   * the blocker is anticollapsed.
   */
  void compute_link_blockers(const ComplexType & parent_complex,
                             const Simplex& alpha_parent,
                             bool is_alpha_blocker = false) {
    for (auto x_link : this->vertex_range()) {
      Vertex_handle x_parent = *this->give_equivalent_vertex(parent_complex,
                                                             x_link);

      for (auto blocker_parent : parent_complex.const_blocker_range(x_parent)) {
        if (!is_alpha_blocker || *blocker_parent != alpha_parent) {
          Simplex sigma_parent(*blocker_parent);

          sigma_parent.difference(alpha_parent);

          if (sigma_parent.dimension() >= 2
              && sigma_parent.first_vertex() == x_parent) {
            Root_simplex_handle sigma_id(parent_complex.get_id(sigma_parent));
            auto sigma_link = this->get_simplex_address(sigma_id);
            // ie if the vertices of sigma are vertices of the link
            if (sigma_link) {
              bool is_new_blocker = true;
              for (auto a : alpha_parent) {
                for (auto eta_parent : parent_complex.const_blocker_range(a)) {
                  if (!is_alpha_blocker || *eta_parent != alpha_parent) {
                    Simplex eta_minus_alpha(*eta_parent);
                    eta_minus_alpha.difference(alpha_parent);
                    if (eta_minus_alpha != sigma_parent
                        && sigma_parent.contains_difference(*eta_parent,
                                                            alpha_parent)) {
                      is_new_blocker = false;
                      break;
                    }
                  }
                }
                if (!is_new_blocker)
                  break;
              }
              if (is_new_blocker)
                this->add_blocker(new Simplex(*sigma_link));
            }
          }
        }
      }
    }
  }

 public:
  /**
   * @brief compute vertices, edges and blockers of the link.
   * @details If the boolean only_superior_vertices is true, then the link is computed only
   * with vertices that are greater than  vertices of alpha_parent_adress.
   */
  void build_link(const ComplexType & parent_complex,
                  const Simplex& alpha_parent_adress,
                  bool is_alpha_blocker = false,
                  bool only_vertices = false) {
    assert(is_alpha_blocker || parent_complex.contains(alpha_parent_adress));
    compute_link_vertices(parent_complex, alpha_parent_adress, only_superior_vertices_);
    if (!only_vertices) {
      compute_link_edges(parent_complex, alpha_parent_adress, is_alpha_blocker);
      compute_link_blockers(parent_complex, alpha_parent_adress, is_alpha_blocker);
    }
  }

  /**
   * @brief build the link of a blocker which is the link
   * of the blocker's simplex if this simplex had been
   * removed from the blockers of the complex.
   */
  friend void build_link_of_blocker(const ComplexType & parent_complex,
                                    Simplex& blocker,
                                    Skeleton_blocker_link_complex & result) {
    assert(blocker.dimension() >= 2);
    assert(parent_complex.contains_blocker(blocker));
    result.clear();
    result.build_link(parent_complex, blocker, true);
  }
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_LINK_COMPLEX_H_
