/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA
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

#ifndef SKELETON_BLOCKER_SKELETON_BLOCKER_LINK_SUPERIOR_H_
#define SKELETON_BLOCKER_SKELETON_BLOCKER_LINK_SUPERIOR_H_

#include <gudhi/Skeleton_blocker_link_complex.h>

namespace Gudhi {

namespace skeleton_blocker {

template<class ComplexType> class Skeleton_blocker_sub_complex;

/**
 *  \brief Class representing the link of a simplicial complex encoded by a skeleton/blockers pair.
 *  It computes only vertices greater than the simplex used to build the link.
 */
template<typename ComplexType>
class Skeleton_blocker_link_superior : public Skeleton_blocker_link_complex<
ComplexType> {
  typedef typename ComplexType::Edge_handle Edge_handle;

  typedef typename ComplexType::boost_vertex_handle boost_vertex_handle;

 public:
  typedef typename ComplexType::Vertex_handle Vertex_handle;
  typedef typename ComplexType::Root_vertex_handle Root_vertex_handle;
  typedef typename ComplexType::Simplex Simplex;
  typedef typename ComplexType::Root_simplex_handle Root_simplex_handle;
  typedef typename ComplexType::BlockerMap BlockerMap;
  typedef typename ComplexType::BlockerPair BlockerPair;
  typedef typename ComplexType::BlockerMapIterator BlockerMapIterator;
  typedef typename ComplexType::BlockerMapConstIterator BlockerMapConstIterator;
  typedef typename ComplexType::Simplex::Simplex_vertex_const_iterator AddressSimplexConstIterator;
  typedef typename ComplexType::Root_simplex_handle::Simplex_vertex_const_iterator IdSimplexConstIterator;

  Skeleton_blocker_link_superior()
      : Skeleton_blocker_link_complex<ComplexType>(true) { }

  Skeleton_blocker_link_superior(const ComplexType & parent_complex,
                                 Simplex& alpha_parent_adress)
      : Skeleton_blocker_link_complex<ComplexType>(parent_complex,
                                                   alpha_parent_adress, true) { }

  Skeleton_blocker_link_superior(const ComplexType & parent_complex,
                                 Vertex_handle a_parent_adress)
      : Skeleton_blocker_link_complex<ComplexType>(parent_complex,
                                                   a_parent_adress, true) { }
};

}  // namespace skeleton_blocker

namespace skbl = skeleton_blocker;

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_SKELETON_BLOCKER_LINK_SUPERIOR_H_
