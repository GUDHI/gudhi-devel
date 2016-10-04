/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 INRIA
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

#ifndef GUDHI_SPATIAL_TREE_DS_H_
#define GUDHI_SPATIAL_TREE_DS_H_

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <cstddef>
#include <vector>

namespace Gudhi {
namespace spatial_searching {


  /**
  * \class Kd_tree_search Kd_tree_search.h gudhi/Kd_tree_search.h
  * \brief Spatial tree data structure to perform (approximate) nearest and farthest neighbor search.
  *
  * \ingroup spatial_searching
  *
  * \details
  * The class Kd_tree_search is a tree-based data structure, based on
  * <a target="_blank" href="http://doc.cgal.org/latest/Spatial_searching/index.html">CGAL dD spatial searching data structures</a>.
  * It provides a simplified API to perform (approximate) nearest and farthest neighbor searches. Contrary to CGAL default behavior, the tree
  * does not store the points themselves, but stores indices.
  *
  * There are two types of queries: the <i>k-nearest or k-farthest neighbor query</i>, where <i>k</i> is fixed and the <i>k</i> nearest 
  * or farthest points are computed right away,
  * and the <i>incremental nearest or farthest neighbor query</i>, where no number of neighbors is provided during the call, as the
  * neighbors will be computed incrementally when the iterator on the range is incremented.
  *
  * \tparam K must be a model of the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Spatial_searching/classSearchTraits.html">SearchTraits</a>
  *   concept, such as the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class, which
  *   can be static if you know the ambiant dimension at compile-time, or dynamic if you don't.
  * \tparam Point_range is the type of the range that provides the points.
  *   It must be a range whose iterator type is a `RandomAccessIterator`.
  */
template <typename K, typename Point_range>
class Kd_tree_search
{
  typedef boost::iterator_property_map<
    typename Point_range::const_iterator,
    CGAL::Identity_property_map<std::ptrdiff_t> >           Point_property_map;

public:
  /// The kernel.
  typedef K                                                 Kernel;
  /// Number type used for distances.
  typedef typename Kernel::FT                               FT;
  /// The point type.
  typedef typename Point_range::value_type                  Point;

  typedef CGAL::Search_traits<
    FT, Point,
    typename Kernel::Cartesian_const_iterator_d, 
    typename Kernel::Construct_cartesian_const_iterator_d>  Traits_base;
  
  typedef CGAL::Search_traits_adapter<
    std::ptrdiff_t,
    Point_property_map,
    Traits_base>                                            STraits;
  
  typedef CGAL::Orthogonal_k_neighbor_search<STraits>       K_neighbor_search;
  typedef typename K_neighbor_search::Tree                  Tree;
  typedef typename K_neighbor_search::Distance              Distance;
  /// \brief The range returned by a k-nearest or k-farthest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef K_neighbor_search                                 KNS_range;

  typedef CGAL::Orthogonal_incremental_neighbor_search<
    STraits, Distance, CGAL::Sliding_midpoint<STraits>, Tree>
                                                   Incremental_neighbor_search;
  /// \brief The range returned by an incremental nearest or farthest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef Incremental_neighbor_search                       INS_range;

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  Kd_tree_search(Point_range const& points)
  : m_points(points),
    m_tree(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(points.size()),
           typename Tree::Splitter(),
           STraits(std::begin(points)) )
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  /// @param[in] only_these_points Specifies the indices of the points that
  /// should be actually inserted into the tree. The other points are ignored.
  template <typename Point_indices_range>
  Kd_tree_search(
    Point_range const& points,
    Point_indices_range const& only_these_points)
    : m_points(points),
      m_tree(
        only_these_points.begin(), only_these_points.end(),
        typename Tree::Splitter(),
        STraits(std::begin(points)))
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  /// @param[in] begin_idx, past_the_end_idx Define the subset of the points that
  /// should be actually inserted into the tree. The other points are ignored.
  Kd_tree_search(
    Point_range const& points,
    std::size_t begin_idx, std::size_t past_the_end_idx)
  : m_points(points),
    m_tree(
      boost::counting_iterator<std::ptrdiff_t>(begin_idx),
      boost::counting_iterator<std::ptrdiff_t>(past_the_end_idx),
      typename Tree::Splitter(),
      STraits(std::begin(points)) )
  {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  // Be careful, this function invalidates the tree,
  // which will be recomputed at the next query
  void insert(std::ptrdiff_t point_idx)
  {
    m_tree.insert(point_idx);
  }

  /// \brief Search for the k-nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] k Number of nearest points to search.
  /// @param[in] sorted Indicates if the computed sequence of k-nearest neighbors needs to be sorted.
  /// @param[in] eps Approximation factor.
  /// @return A range containing the k-nearest neighbors.
  KNS_range query_k_nearest_neighbors(const
    Point &p,
    unsigned int k,
    bool sorted = true,
    FT eps = FT(0)) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      p,
      k,
      eps,
      true,
      CGAL::Distance_adapter<std::ptrdiff_t,Point_property_map,CGAL::Euclidean_distance<Traits_base> >(
        std::begin(m_points)),
      sorted);

    return search;
  }

  /// \brief Search incrementally for the nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] eps Approximation factor.
  /// @return A range containing the neighbors sorted by their distance to p. 
  /// All the neighbors are not computed by this function, but they will be
  /// computed incrementally when the iterator on the range is incremented.
  INS_range query_incremental_nearest_neighbors(const Point &p, FT eps = FT(0)) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      p,
      eps,
      true,
      CGAL::Distance_adapter<std::ptrdiff_t, Point_property_map, CGAL::Euclidean_distance<Traits_base> >(
        std::begin(m_points)) );

    return search;
  }

  /// \brief Search for the k-farthest points from a query point.
  /// @param[in] p The query point.
  /// @param[in] k Number of farthest points to search.
  /// @param[in] sorted Indicates if the computed sequence of k-farthest neighbors needs to be sorted.
  /// @param[in] eps Approximation factor.
  /// @return A range containing the k-farthest neighbors.
  KNS_range query_k_farthest_neighbors(const
    Point &p,
    unsigned int k,
    bool sorted = true,
    FT eps = FT(0)) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      p,
      k,
      eps,
      false,
      CGAL::Distance_adapter<std::ptrdiff_t,Point_property_map,CGAL::Euclidean_distance<Traits_base> >(
        std::begin(m_points)),
      sorted);

    return search;
  }

  /// \brief Search incrementally for the farthest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] eps Approximation factor.
  /// @return A range containing the neighbors sorted by their distance to p. 
  /// All the neighbors are not computed by this function, but they will be
  /// computed incrementally when the iterator on the range is incremented.
  INS_range query_incremental_farthest_neighbors(const Point &p, FT eps = FT(0)) const
  {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      p,
      eps,
      false,
      CGAL::Distance_adapter<std::ptrdiff_t, Point_property_map, CGAL::Euclidean_distance<Traits_base> >(
        std::begin(m_points)) );

    return search;
  }


private:
  Point_range const& m_points;
  Tree m_tree;
};

} // namespace spatial_searching
} // namespace Gudhi

#endif // GUDHI_SPATIAL_TREE_DS_H_
