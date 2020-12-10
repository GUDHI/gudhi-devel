/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL and Eigen3
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef KD_TREE_SEARCH_H_
#define KD_TREE_SEARCH_H_

#include <gudhi/Debug_utils.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

#include <boost/property_map/property_map.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <cstddef>
#include <vector>

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Kd_tree_search is only available for CGAL >= 4.11
#endif

#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Kd_tree_search is only available for Eigen3 >= 3.1.0 installed with CGAL
#endif

namespace Gudhi {
namespace spatial_searching {

  /**
  * \class Kd_tree_search Kd_tree_search.h gudhi/Kd_tree_search.h
  * \brief Spatial tree data structure to perform (approximate) nearest and furthest neighbor search.
  *
  * \ingroup spatial_searching
  *
  * \details
  * The class Kd_tree_search is a tree-based data structure, based on
  * <a target="_blank" href="http://doc.cgal.org/latest/Spatial_searching/index.html">CGAL dD spatial searching data structures</a>.
  * It provides a simplified API to perform (approximate) nearest and furthest neighbor searches. Contrary to CGAL default behavior, the tree
  * does not store the points themselves, but stores indices.
  *
  * There are two types of queries: the <i>k-nearest or k-furthest neighbor query</i>, where <i>k</i> is fixed and the <i>k</i> nearest 
  * or furthest points are computed right away,
  * and the <i>incremental nearest or furthest neighbor query</i>, where no number of neighbors is provided during the call, as the
  * neighbors will be computed incrementally when the iterator on the range is incremented.
  *
  * \tparam Search_traits must be a model of the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Spatial_searching/classSearchTraits.html">SearchTraits</a>
  *   concept, such as the <a target="_blank"
  *   href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class, which
  *   can be static if you know the ambiant dimension at compile-time, or dynamic if you don't.
  * \tparam Point_range is the type of the range that provides the points.
  *   It must be a range whose iterator type is a `RandomAccessIterator`.
  */
template <typename Search_traits, typename Point_range>
class Kd_tree_search {
  typedef boost::iterator_property_map<
    typename Point_range::const_iterator,
    CGAL::Identity_property_map<std::ptrdiff_t> >           Point_property_map;

 public:
  /// The Traits.
  typedef Search_traits                                     Traits;
  /// Number type used for distances.
  typedef typename Traits::FT                               FT;
  /// The point type.
  typedef typename Point_range::value_type                  Point;

  typedef CGAL::Search_traits<
    FT, Point,
    typename Traits::Cartesian_const_iterator_d,
    typename Traits::Construct_cartesian_const_iterator_d,
    typename Traits::Dimension>  Traits_base;

  typedef CGAL::Search_traits_adapter<
    std::ptrdiff_t,
    Point_property_map,
    Traits_base>                                            STraits;
  typedef CGAL::Distance_adapter<
    std::ptrdiff_t,
    Point_property_map,
    CGAL::Euclidean_distance<Traits_base> >                 Orthogonal_distance;

  typedef CGAL::Orthogonal_k_neighbor_search<STraits>       K_neighbor_search;
  typedef typename K_neighbor_search::Tree                  Tree;
  typedef typename K_neighbor_search::Distance              Distance;
  /// \brief The range returned by a k-nearest or k-furthest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef K_neighbor_search                                 KNS_range;

  typedef CGAL::Orthogonal_incremental_neighbor_search<
    STraits, Distance, CGAL::Sliding_midpoint<STraits>, Tree>
                                                   Incremental_neighbor_search;
  /// \brief The range returned by an incremental nearest or furthest neighbor search.
  /// Its value type is `std::pair<std::size_t, FT>` where `first` is the index
  /// of a point P and `second` is the squared distance between P and the query point.
  typedef Incremental_neighbor_search                       INS_range;

  // Because CGAL::Fuzzy_sphere takes the radius and not its square
  struct Sphere_for_kdtree_search
  {
    typedef typename Traits::Point_d Point_d;
    typedef typename Traits::FT FT;
    typedef typename Traits::Dimension D;
    typedef D Dimension;

    private:
    STraits traits;
    Point_d c;
    FT sqradmin, sqradmax;
    bool use_max;

    public:
    // `prefer_max` means that we prefer outputting more points at squared distance between r2min and r2max,
    // while `!prefer_max` means we prefer fewer.
    Sphere_for_kdtree_search(Point_d const& c_, FT const& r2min, FT const& r2max, bool prefer_max=true, STraits const& traits_ = {})
      : traits(traits_), c(c_), sqradmin(r2min), sqradmax(r2max), use_max(prefer_max)
      { GUDHI_CHECK(r2min >= 0 && r2max >= r2min, "0 <= r2min <= r2max"); }

    bool contains(std::ptrdiff_t i) const {
      const Point_d& p = get(traits.point_property_map(), i);
      auto ccci = traits.construct_cartesian_const_iterator_d_object();
      return contains_point_given_as_coordinates(ccci(p), ccci(p, 0));
    }

    template <typename Coord_iterator>
      bool contains_point_given_as_coordinates(Coord_iterator pi, Coord_iterator CGAL_UNUSED) const {
        FT distance = 0;
        auto ccci = traits.construct_cartesian_const_iterator_d_object();
        auto ci = ccci(c);
        auto ce = ccci(c, 0);
        FT const& limit = use_max ? sqradmax : sqradmin;
        while (ci != ce) {
          distance += CGAL::square(*pi++ - *ci++);
          // I think Clément advised to check the distance at every step instead of
          // just at the end, especially when the dimension becomes large. Distance
          // isn't part of the concept anyway.
          if (distance > limit) return false;
        }
        return true;
      }

    bool inner_range_intersects(CGAL::Kd_tree_rectangle<FT, D> const& rect) const {
      auto ccci = traits.construct_cartesian_const_iterator_d_object();
      FT distance = 0;
      auto ci = ccci(c);
      auto ce = ccci(c, 0);
      for (int i = 0; ci != ce; ++i, ++ci) {
        distance += CGAL::square(CGAL::max<FT>(CGAL::max<FT>(*ci - rect.max_coord(i), rect.min_coord(i) - *ci), 0 ));
        if (distance > sqradmin) return false;
      }
      return true;
    }


    bool outer_range_contains(CGAL::Kd_tree_rectangle<FT, D> const& rect) const {
      auto ccci = traits.construct_cartesian_const_iterator_d_object();
      FT distance = 0;
      auto ci = ccci(c);
      auto ce = ccci(c, 0);
      for (int i = 0; ci != ce; ++i, ++ci) {
        distance += CGAL::square(CGAL::max<FT>(*ci - rect.min_coord(i), rect.max_coord(i) - *ci));
        if (distance > sqradmax) return false;
      }
      return true;
    }
  };

  /// \brief Constructor
  /// @param[in] points Const reference to the point range. This range
  /// is not copied, so it should not be destroyed or modified afterwards.
  Kd_tree_search(Point_range const& points)
  : m_points(points),
    m_tree(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(points.size()),
           typename Tree::Splitter(),
           STraits(std::begin(points))) {
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
        STraits(std::begin(points))) {
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
      STraits(std::begin(points))) {
    // Build the tree now (we don't want to wait for the first query)
    m_tree.build();
  }

  // Be careful, this function invalidates the tree,
  // which will be recomputed at the next query
  void insert(std::ptrdiff_t point_idx) {
    m_tree.insert(point_idx);
  }

  /// \brief Search for the k-nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] k Number of nearest points to search.
  /// @param[in] sorted Indicates if the computed sequence of k-nearest neighbors needs to be sorted.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`) containing the k-nearest neighbors.
  KNS_range k_nearest_neighbors(
    Point const& p,
    unsigned int k,
    bool sorted = true,
    FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      p,
      k,
      eps,
      true,
      Orthogonal_distance(std::begin(m_points)), sorted);

    return search;
  }

  /// \brief Search incrementally for the nearest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`) containing the
  /// neighbors sorted by their distance to p.
  /// All the neighbors are not computed by this function, but they will be
  /// computed incrementally when the iterator on the range is incremented.
  INS_range incremental_nearest_neighbors(Point const& p, FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      p,
      eps,
      true,
      Orthogonal_distance(std::begin(m_points)) );

    return search;
  }

  /// \brief Search for the k-furthest points from a query point.
  /// @param[in] p The query point.
  /// @param[in] k Number of furthest points to search.
  /// @param[in] sorted Indicates if the computed sequence of k-furthest neighbors needs to be sorted.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`) containing the k-furthest neighbors.
  KNS_range k_furthest_neighbors(
    Point const& p,
    unsigned int k,
    bool sorted = true,
    FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    K_neighbor_search search(
      m_tree,
      p,
      k,
      eps,
      false,
      Orthogonal_distance(std::begin(m_points)), sorted);

    return search;
  }

  /// \brief Search incrementally for the furthest neighbors from a query point.
  /// @param[in] p The query point.
  /// @param[in] eps Approximation factor.
  /// @return A range (whose `value_type` is `std::size_t`)
  /// containing the neighbors sorted by their distance to p.
  /// All the neighbors are not computed by this function, but they will be
  /// computed incrementally when the iterator on the range is incremented.
  INS_range incremental_furthest_neighbors(Point const& p, FT eps = FT(0)) const {
    // Initialize the search structure, and search all N points
    // Note that we need to pass the Distance explicitly since it needs to
    // know the property map
    Incremental_neighbor_search search(
      m_tree,
      p,
      eps,
      false,
      Orthogonal_distance(std::begin(m_points)) );

    return search;
  }

  /// \brief Search for all the neighbors in a ball.
  /// @param[in] p The query point.
  /// @param[in] radius The search radius
  /// @param[out] it The points that lie inside the sphere of center `p` and radius `radius`.
  ///                Note: `it` is used this way: `*it++ = each_point`.
  /// @param[in] eps Approximation factor.
  template <typename OutputIterator>
  void all_near_neighbors(Point const& p,
                   FT const& radius,
                   OutputIterator it,
                   FT eps = FT(0)) const {
    all_near_neighbors2(p, CGAL::square(radius - eps), CGAL::square(radius + eps), it);
  }

  /// \brief Search for all the neighbors in a ball. This is similar to `all_near_neighbors` but takes directly
  /// the square of the minimum distance below which points must be considered neighbors and square of the
  /// maximum distance above which they cannot be.
  /// @param[in] p The query point.
  /// @param[in] sq_radius_min The square of the minimum search radius
  /// @param[in] sq_radius_max The square of the maximum search radius
  /// @param[out] it The points that lie inside the sphere of center `p` and squared radius `sq_radius`.
  ///                Note: `it` is used this way: `*it++ = each_point`.
  template <typename OutputIterator>
  void all_near_neighbors2(Point const& p,
                   FT const& sq_radius_min,
                   FT const& sq_radius_max,
                   OutputIterator it) const {
    m_tree.search(it, Sphere_for_kdtree_search(p, sq_radius_min, sq_radius_max, true, m_tree.traits()));
  }

  int tree_depth() const {
    return m_tree.root()->depth();
  }

 private:
  Point_range const& m_points;
  Tree m_tree;
};

}  // namespace spatial_searching
}  // namespace Gudhi

#endif  // KD_TREE_SEARCH_H_
