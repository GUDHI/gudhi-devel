/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SKELETON_BLOCKER_SKELETON_BLOCKER_SUB_COMPLEX_H_
#define SKELETON_BLOCKER_SKELETON_BLOCKER_SUB_COMPLEX_H_

#include <gudhi/Skeleton_blocker_complex.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_simplex.h>
#include <gudhi/Debug_utils.h>

#include <map>
#include <vector>

namespace Gudhi {

namespace skeleton_blocker {

/**
 * @brief Simplicial subcomplex of a complex represented by a skeleton/blockers pair.
 * @extends Skeleton_blocker_complex
 * @details Stores a subcomplex of a simplicial complex.
 * To simplify explanations below, we will suppose that :
 * - K is the root simplicial complex
 * - L is a subcomplex of K.
 *
 * One vertex of K may exists in L but with a different address.
 * To be able to locate the vertices in K from vertices of L, the class
 * stores a map 'addresses' between vertices of K and vertices of L.
 *
 * Note that the type for handle of vertices of L is 'Vertex_handle' and
 * the type for handle of vertices of K is 'Root_vertex_handle'.
 *
 * The template ComplexType is type of the root complex. It allows to know
 * if the subcomplex is geometric or not.
 * It has to be either 'Skeleton_blockers_complex' or 'Skeleton_blockers_geometric_complex'.
 *
 */
template<typename ComplexType>
class Skeleton_blocker_sub_complex : public ComplexType {
 protected:
  template<class T> friend class Skeleton_blocker_link_complex;

  typedef typename ComplexType::Graph Graph;
  typedef typename ComplexType::Edge_handle Edge_handle;

  typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

 public:
  using ComplexType::add_vertex;
  using ComplexType::add_edge_without_blockers;
  using ComplexType::add_blocker;

  typedef typename ComplexType::Vertex_handle Vertex_handle;
  typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;
  typedef typename ComplexType::Simplex Simplex;
  typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;

 protected:
  /**
   * @brief Determines whether all proper faces of simplex 'sigma' belong to 'link1' \cup 'link2'
   * where 'link1' and 'link2' are subcomplexes of the same complex of type ComplexType
   */
  typedef std::map<Root_vertex_handle, Vertex_handle> IdAddressMap;
  typedef typename IdAddressMap::value_type AddressPair;
  typedef typename IdAddressMap::iterator IdAddressMapIterator;
  typedef typename IdAddressMap::const_iterator IdAddressMapConstIterator;
  std::map<Root_vertex_handle, Vertex_handle> addresses;

 public:
  /**
   * Add a vertex 'global' of K to L. When added to L, this vertex will receive
   * another number, addresses(global), its local address.
   * return the address where the vertex lay on L.
   * The vertex corresponding to 'global' must not be already present
   * in the complex.
   */
  Vertex_handle add_vertex(Root_vertex_handle global) {
    assert(!this->contains_vertex(global));
    Vertex_handle address(boost::add_vertex(this->skeleton));
    this->num_vertices_++;
    (*this)[address].activate();
    (*this)[address].set_id(global);
    addresses.insert(AddressPair(global, address));
    this->degree_.push_back(0);
    return address;
  }

  /**
   * Add an edge (v1_root,v2_root) to the sub-complex.
   * It assumes that both vertices corresponding to v1_root and v2_root are present
   * in the sub-complex.
   */
  void add_edge_without_blockers(Root_vertex_handle v1_root, Root_vertex_handle v2_root) {
    auto v1_sub(this->get_address(v1_root));
    auto v2_sub(this->get_address(v2_root));
    assert(v1_sub && v2_sub);
    this->ComplexType::add_edge_without_blockers(*v1_sub, *v2_sub);
  }

  /**
   * Add a blocker to the sub-complex.
   * It assumes that all vertices of blocker_root are present
   * in the sub-complex.
   */
  void add_blocker(const Root_simplex_handle& blocker_root) {
    auto blocker_sub = this->get_address(blocker_root);
    assert(blocker_sub);
    this->add_blocker(new Simplex(*blocker_sub));
  }

 public:
  /**
   * Constructs the restricted complex of 'parent_complex' to
   * vertices of 'simplex'.
   */
  void make_restricted_complex(const ComplexType & parent_complex,
                               const Simplex& simplex) {
    this->clear();
    // add vertices to the sub complex
    for (auto x : simplex) {
      assert(parent_complex.contains_vertex(x));
      auto x_local = this->add_vertex(parent_complex[x].get_id());
      (*this)[x_local] = parent_complex[x];
    }

    // add edges to the sub complex
    for (auto x : simplex) {
      // x_neigh is the neighbor of x intersected with vertices_simplex
      Simplex x_neigh;
      parent_complex.add_neighbours(x, x_neigh, true);
      x_neigh.intersection(simplex);
      for (auto y : x_neigh) {
        this->add_edge_without_blockers(parent_complex[x].get_id(), parent_complex[y].get_id());
      }
    }

    // add blockers to the sub complex
    for (auto blocker : parent_complex.const_blocker_range()) {
      // check if it is the first time we encounter the blocker
      if (simplex.contains(*blocker)) {
        Root_simplex_handle blocker_root(parent_complex.get_id(*(blocker)));
        Simplex blocker_restr(
                              *(this->get_simplex_address(blocker_root)));
        this->add_blocker(new Simplex(blocker_restr));
      }
    }
  }

  void clear() {
    addresses.clear();
    ComplexType::clear();
  }

  /**
   * Compute the local vertex in L corresponding to the vertex global in K.
   * runs in O(log n) if n = num_vertices()
   */
  boost::optional<Vertex_handle> get_address(Root_vertex_handle global) const {
    boost::optional < Vertex_handle > res;
    IdAddressMapConstIterator it = addresses.find(global);
    if (it == addresses.end())
      res.reset();
    else
      res = (*it).second;
    return res;
  }

  // /**
  //  * Allocates a simplex in L corresponding to the simplex s in K
  //  * with its local addresses and returns an AddressSimplex.
  //  */
  // boost::optional<Simplex> get_address(const Root_simplex_handle & s) const;

  // private:
  /**
   *  same as get_address except that it will return a simplex in any case.
   *  The vertices that were not found are not added.
   */
  // @remark should be private but problem with VS

  std::vector<boost::optional<Vertex_handle> > get_addresses(
                                                             const Root_simplex_handle & s) const {
    std::vector < boost::optional<Vertex_handle> > res;
    for (auto i : s) {
      res.push_back(get_address(i));
    }
    return res;
  }
};

/**
 * @remark waste of time to create a new simplex each time when we could use instead of addresses_sigma_in_link a
 * simplex with special values (ComplexDS::null_vertex e.g.) to indicate that a vertex does not belong to the complex.
 */
template<typename ComplexType>
bool proper_face_in_union(
                          Skeleton_blocker_sub_complex<ComplexType> & link,
                          std::vector<boost::optional<typename ComplexType::Vertex_handle> > & addresses_sigma_in_link,
                          std::size_t vertex_to_be_ignored) {
  // we test that all vertices of 'addresses_sigma_in_link' but 'vertex_to_be_ignored'
  // are in link1 if it is the case we construct the corresponding simplex
  bool vertices_sigma_are_in_link = true;
  typename ComplexType::Simplex sigma_in_link;
  for (std::size_t i = 0; i < addresses_sigma_in_link.size(); ++i) {
    if (i != vertex_to_be_ignored) {
      if (!addresses_sigma_in_link[i]) {
        vertices_sigma_are_in_link = false;
        break;
      } else {
        sigma_in_link.add_vertex(*addresses_sigma_in_link[i]);
      }
    }
  }
  // If one of vertices of the simplex is not in the complex then it returns false
  // Otherwise, it tests if the simplex is in the complex
  return vertices_sigma_are_in_link && link.contains(sigma_in_link);
}

// Remark: this function should be friend in order to leave get_addresses private
// however doing so seems currently not possible due to a visual studio bug c2668
// "the compiler does not support partial ordering of template functions as specified in the C++ Standard"
// http://www.serkey.com/error-c2668-ambiguous-call-to-overloaded-function-bb45ft.html

template<typename ComplexType>
bool proper_faces_in_union(
                           Skeleton_blocker_simplex<typename ComplexType::Root_vertex_handle> & sigma,
                           Skeleton_blocker_sub_complex<ComplexType> & link1,
                           Skeleton_blocker_sub_complex<ComplexType> & link2) {
  typedef typename ComplexType::Vertex_handle Vertex_handle;
  std::vector < boost::optional<Vertex_handle> > addresses_sigma_in_link1 =
      link1.get_addresses(sigma);
  std::vector < boost::optional<Vertex_handle> > addresses_sigma_in_link2 =
      link2.get_addresses(sigma);

  for (std::size_t current_index = 0; current_index < addresses_sigma_in_link1.size();
       ++current_index) {
    if (!proper_face_in_union(link1, addresses_sigma_in_link1, current_index)
        && !proper_face_in_union(link2, addresses_sigma_in_link2,
                                 current_index)) {
      return false;
    }
  }
  return true;
}

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_SKELETON_BLOCKER_SUB_COMPLEX_H_
