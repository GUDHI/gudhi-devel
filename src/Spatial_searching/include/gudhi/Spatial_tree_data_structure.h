/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016  INRIA Sophia-Antipolis (France)
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

#ifndef GUDHI_POINT_CLOUD_H
#define GUDHI_POINT_CLOUD_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>

#include <boost/iterator/counting_iterator.hpp>

#include <cstddef>
#include <vector>

namespace Gudhi {

template <typename K, typename Point_container_>
class Spatial_tree_data_structure
{
public:
  typedef typename Point_container_::value_type             Point;
  typedef K                                                 Kernel;
  typedef typename Kernel::FT                               FT;

  typedef CGAL::Search_traits<
    FT, Point,
    typename Kernel::Cartesian_const_iterator_d, 
    typename Kernel::Construct_cartesian_const_iterator_d>  Traits_base;
  // using a pointer as a special property map type
  typedef CGAL::Search_traits_adapter<
    std::ptrdiff_t, Point*, Traits_base>                    STraits;
  
  typedef CGAL::Orthogonal_k_neighbor_search<STraits>       K_neighbor_search;
  typedef typename K_neighbor_search::Tree                  Tree;
  typedef typename K_neighbor_search::Distance              Distance;
  typedef typename K_neighbor_search::iterator              KNS_iterator;
  typedef K_neighbor_search                                 KNS_range;

  typedef CGAL::Orthogonal_incremental_neighbor_search<
    STraits, Distance, CGAL::Sliding_midpoint<STraits>, Tree>
                                                   Incremental_neighbor_search;
  typedef typename Incremental_neighbor_search::iterator    INS_iterator;
  typedef Incremental_neighbor_search                       INS_range;

  /// Constructor
  Spatial_tree_data_structure(Point_container_ const& points)
  : m_points(points),
    m_tree(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(points.size()),
           typename Tree::Splitter(),
           STraits((Point*)&(points[0])) )
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// Constructor
  template <typename Point_indices_range>
  Spatial_tree_data_structure(
    Point_container_ const& points,
    Point_indices_range const& only_these_points)
    : m_points(points),
    m_tree(
    only_these_points.begin(), only_these_points.end(),
    typename Tree::Splitter(),
    STraits((Point*)&(points[0])))
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// Constructor
  Spatial_tree_data_structure(
    Point_container_ const& points,
    std::size_t begin_idx, std::size_t past_the_end_idx)
  : m_points(points),
    m_tree(
      boost::counting_iterator<std::ptrdiff_t>(begin_idx),
      boost::counting_iterator<std::ptrdiff_t>(past_the_end_idx),
      typename Tree::Splitter(),
      STraits((Point*)&(points[0])) )
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /*Point_container_ &points()
  {
    return m_points;
  }

  const Point_container_ &points() const
  {
    return m_points;
  }*/

  // Be careful, this function invalidates the tree,
  // which will be recomputed at the next query
  void insert(std::ptrdiff_t point_idx)
  {
    m_tree.insert(point_idx);
  }

  KNS_range query_ANN(const
    Point &sp,
    unsigned int k,
    bool sorted = true) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      sp,
      k,
      FT(0),
      true,
      CGAL::Distance_adapter<std::ptrdiff_t,Point*,CGAL::Euclidean_distance<Traits_base> >(
        (Point*)&(m_points[0])),
      sorted);

    return search;
  }

  INS_range query_incremental_ANN(const Point &sp) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      sp,
      FT(0),
      true,
      CGAL::Distance_adapter<std::ptrdiff_t, Point*, CGAL::Euclidean_distance<Traits_base> >(
        (Point*)&(m_points[0])) );

    return search;
  }

protected:
  Point_container_ const& m_points;
  Tree m_tree;
};

} //namespace Gudhi

#endif // GUDHI_POINT_CLOUD_H
