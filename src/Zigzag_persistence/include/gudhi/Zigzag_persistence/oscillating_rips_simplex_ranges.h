/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file oscillating_rips_simplex_ranges.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Oscillating_rips_edge_order_policy enum,
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_simplex_iterator_base class,
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_simplex_iterator_range class and
 * @ref Gudhi::zigzag_persistence::Oscillating_rips_simplex_vector_range_constructor class.
 */

#ifndef ZIGZAG_OSCILLATING_RIPS_SIMPLEX_RANGES_H_
#define ZIGZAG_OSCILLATING_RIPS_SIMPLEX_RANGES_H_

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>

#ifdef GUDHI_USE_TBB
#include <tbb/tbb.h>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Zigzag_persistence/Zigzag_edge.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @private
 * @brief
 *
 * @tparam Filtration_value
 */
class Oscillating_rips_edge_expander
{
 public:
  /**
   * @brief Reverse lexicographical order for the simplex handles.
   */
  template <class StableFilteredComplex>
  struct reverse_lexicographic_order {
    using Simplex_handle = typename StableFilteredComplex::Simplex_handle;

    explicit reverse_lexicographic_order(StableFilteredComplex* st) : st_(st) {}

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const
    {
      auto rg1 = st_->simplex_vertex_range(sh1);
      auto rg2 = st_->simplex_vertex_range(sh2);
      auto it1 = rg1.begin();
      auto it2 = rg2.begin();
      while (it1 != rg1.end() && it2 != rg2.end()) {
        if (*it1 == *it2) {
          ++it1;
          ++it2;
        } else {
          return *it1 < *it2;
        }
      }
      return ((it1 == rg1.end()) && (it2 != rg2.end()));
    }

    StableFilteredComplex* st_;
  };

  template <class StableFilteredComplex, typename EdgeRangeIterator>
  static void expand_edges_of_same_filtration_value(
      EdgeRangeIterator& currentEdgeIt,
      const EdgeRangeIterator& edgeEndIterator,
      int maxDimension,
      StableFilteredComplex& complex,
      std::vector<typename StableFilteredComplex::Simplex_handle>& currentSimplices)
  {
    if (currentEdgeIt == edgeEndIterator) return;

    auto fil = currentEdgeIt->get_filtration_value();
    while (currentEdgeIt != edgeEndIterator && currentEdgeIt->get_direction() &&
           currentEdgeIt->get_filtration_value() == fil) {
      complex.insert_edge_as_flag(currentEdgeIt->get_smallest_vertex(),
                                  currentEdgeIt->get_biggest_vertex(),
                                  fil,
                                  maxDimension,
                                  currentSimplices);
      ++currentEdgeIt;
    }
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(currentSimplices.begin(), currentSimplices.end(), reverse_lexicographic_order(&complex));
#else
    std::sort(currentSimplices.begin(), currentSimplices.end(), reverse_lexicographic_order(&complex));
#endif
  }

  /**
   * @brief Updates @p currentSimplices_ for removals.
   *
   * @param fil Current filtration value.
   */
  template <class StableFilteredComplex, typename EdgeRangeIterator>
  static void collapse_edges_of_same_filtration_value(
      EdgeRangeIterator& currentEdgeIt,
      const EdgeRangeIterator& edgeEndIterator,
      StableFilteredComplex& complex,
      std::vector<typename StableFilteredComplex::Simplex_handle>& currentSimplices)
  {
    using Filtration_value = typename StableFilteredComplex::Filtration_value;
    using Simplex_handle = typename StableFilteredComplex::Simplex_handle;

    if (currentEdgeIt == edgeEndIterator) return;

    unsigned int count = 0;
    Filtration_value fil = currentEdgeIt->get_filtration_value();
    while (currentEdgeIt != edgeEndIterator && !currentEdgeIt->get_direction() &&
           currentEdgeIt->get_filtration_value() == fil) {
      Simplex_handle sh = complex.find({currentEdgeIt->get_smallest_vertex()});
      if (currentEdgeIt->get_smallest_vertex() != currentEdgeIt->get_biggest_vertex()) {
        sh = sh->second.children()->members().find(currentEdgeIt->get_biggest_vertex());
      }
      auto toRemove = complex.star_simplex_range(sh);
      currentSimplices.insert(currentSimplices.end(), toRemove.begin(), toRemove.end());
      ++currentEdgeIt;
      ++count;
    }
#ifdef GUDHI_USE_TBB
    tbb::parallel_sort(
        currentSimplices.begin(), currentSimplices.end(), [&](Simplex_handle sh1, Simplex_handle sh2) -> bool {
          return complex.key(sh1) > complex.key(sh2);
        });
#else
    std::sort(currentSimplices.begin(), currentSimplices.end(), [&](Simplex_handle sh1, Simplex_handle sh2) -> bool {
      return complex.key(sh1) > complex.key(sh2);
    });
#endif
    if (count > 1) {  // more than 1 edge inserted: cofaces can be duplicated
      auto last = std::unique(
          currentSimplices.begin(), currentSimplices.end(), [&](Simplex_handle sh1, Simplex_handle sh2) -> bool {
            return complex.key(sh1) == complex.key(sh2);
          });                                                // equal simplex handles means equal key
      currentSimplices.erase(last, currentSimplices.end());  // remove duplicated cofaces
    }
  }
};

/**
 * @class Oscillating_rips_simplex_iterator oscillating_rips_simplex_ranges.h \
 * gudhi/Zigzag_persistence/oscillating_rips_simplex_ranges.h
 * @brief Custom iterator over the simplices of an oscillating rips filtration. Is a forward iterator only.
 *
 * It inherits from boost::iterator_facade.
 *
 * For each positive edge given by the edge range iterator, the given complex computes and adds all possible
 * cofaces and they respective simplex handle are stored. For each negative edge, the given complex computes
 * the simplex handles of all cofaces and they are stored.
 * The simplex handles induced by an edge are stored only the time needed and are removed when reaching the next edge.
 * That is, we start by storing the possible cofaces of an edge, then at each increment,
 * the next element within the stored simplex handles is retrieved and once the end is reached,
 * the simplex handles are erased and replaced by the possible cofaces of the next edge and so on...
 * If the edge is negative, then a simplex is removed from the complex at the next increment.
 * So, the maximum size of the complex corresponds to the maximum size of a complex in the zigzag filtration.
 *
 * The simplices are returned by the iterator as tuples of three elements:
 * the first is the simplex handle of the simplex in the given complex,
 * the second is the filtration value of the corresponding arrow,
 * the third is the direction of the arrow, i.e., indicates if the simplex is inserted or removed.
 *
 * @warning As the iterator stores possibly large ranges, avoid copying it.
 */
template <class StableFilteredComplex, typename EdgeRangeIterator>
class Oscillating_rips_simplex_iterator_base
{
 public:
 public:
  using Filtration_value = typename StableFilteredComplex::Filtration_value; /**< Filtration value type. */
  using Simplex_handle = typename StableFilteredComplex::Simplex_handle;     /**< Simplex handle type. */
  using Simplex_key = typename StableFilteredComplex::Simplex_key;           /**< Key type. */

  /**
   * @brief Constructor. The edge iterators and the complex are not copied, so do not modify or invalidate
   * them outside as long as this iterator is different from the end iterator.
   *
   * @param edgeStartIterator Begin iterator of the oscillating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param edgeEndIterator End iterator of the oscillating Rips edge range.
   * Is moved, so the original iterator will probably be invalidated.
   * @param complex Empty complex which will be used to store current simplices.
   * Only the address of the complex is stored, so the state of the complex can be
   * consulted between each increment.
   * @param maxDimension Maximal dimension to which to expand the complex. If set to -1, there is no limit.
   * Default value: -1.
   */
  Oscillating_rips_simplex_iterator_base(const EdgeRangeIterator& edgeStartIterator,
                                         const EdgeRangeIterator& edgeEndIterator,
                                         StableFilteredComplex* complex,
                                         int maxDimension = -1)
      : complex_(complex),
        currentSimplexIndex_(0),
        currentEdgeIt_(edgeStartIterator),
        endEdgeIt_(edgeEndIterator),
        currentDirection_(true),
        maxDimension_(maxDimension),
        currentArrowNumber_(0)
  {
    if (currentEdgeIt_ == endEdgeIt_) {
      _set_end();
      return;
    }

    // first simplex is the vertex (0,0) which is the only one with its filtration value.
    complex_->insert_edge_as_flag(currentEdgeIt_->get_smallest_vertex(),
                                  currentEdgeIt_->get_biggest_vertex(),
                                  currentEdgeIt_->get_filtration_value(),
                                  maxDimension_,
                                  currentSimplices_);
    ++currentEdgeIt_;

    std::get<0>(currentArrow_) = currentSimplices_[currentSimplexIndex_];
    std::get<1>(currentArrow_) = complex_->filtration(currentSimplices_[currentSimplexIndex_]);
    std::get<2>(currentArrow_) = currentDirection_;

    complex_->assign_key(currentSimplices_[currentSimplexIndex_], 0);
  }

  /**
   * @brief Default constructor. Corresponds to the end iterator.
   */
  Oscillating_rips_simplex_iterator_base()
      : complex_(nullptr), currentSimplexIndex_(0), currentDirection_(true), maxDimension_(0), currentArrowNumber_(0)
  {
  }

  /**
   * @brief Mandatory for the boost::iterator_facade inheritance. Indicates if two iterators are equal.
   *
   * @param other Iterator to compare.
   * @return True, if the two iterators point to the same position.
   * @return False, otherwise.
   */
  bool equal(Oscillating_rips_simplex_iterator_base const& other) const
  {
    if (complex_ == nullptr) return other.complex_ == nullptr;

    return complex_ == other.complex_ && currentEdgeIt_ == other.currentEdgeIt_ &&
           currentSimplexIndex_ == other.currentSimplexIndex_;
  }

  /**
   * @brief Mandatory for the boost::iterator_facade inheritance. Dereference the iterator.
   *
   * @return Value of the current return element.
   */
  const std::tuple<Simplex_handle, Filtration_value, bool>& dereference() const { return currentArrow_; }

  /**
   * @brief Mandatory for the boost::iterator_facade inheritance. Increments the iterator.
   */
  void increment()
  {
    if (!currentDirection_) complex_->remove_maximal_simplex(currentSimplices_[currentSimplexIndex_]);

    ++currentSimplexIndex_;
    ++currentArrowNumber_;

    if (currentSimplexIndex_ == currentSimplices_.size()) {
      if (currentEdgeIt_ == endEdgeIt_) {
        _set_end();
        return;
      }

      currentSimplices_.clear();

      auto fil = currentEdgeIt_->get_filtration_value();
      currentDirection_ = currentEdgeIt_->get_direction();

      std::get<1>(currentArrow_) = fil;
      std::get<2>(currentArrow_) = currentDirection_;

      if (currentDirection_) {
        Oscillating_rips_edge_expander::expand_edges_of_same_filtration_value(
            currentEdgeIt_, endEdgeIt_, maxDimension_, *complex_, currentSimplices_);
      } else {
        Oscillating_rips_edge_expander::collapse_edges_of_same_filtration_value(
            currentEdgeIt_, endEdgeIt_, *complex_, currentSimplices_);
      }
      currentSimplexIndex_ = 0;
    }

    std::get<0>(currentArrow_) = currentSimplices_[currentSimplexIndex_];
    if (currentDirection_) complex_->assign_key(currentSimplices_[currentSimplexIndex_], currentArrowNumber_);
  }

 private:
  /**
   * @brief Reverse lexicographical order for the simplex handles.
   */
  struct reverse_lexicographic_order {
    explicit reverse_lexicographic_order(StableFilteredComplex* st) : st_(st) {}

    bool operator()(const Simplex_handle sh1, const Simplex_handle sh2) const
    {
      auto rg1 = st_->simplex_vertex_range(sh1);
      auto rg2 = st_->simplex_vertex_range(sh2);
      auto it1 = rg1.begin();
      auto it2 = rg2.begin();
      while (it1 != rg1.end() && it2 != rg2.end()) {
        if (*it1 == *it2) {
          ++it1;
          ++it2;
        } else {
          return *it1 < *it2;
        }
      }
      return ((it1 == rg1.end()) && (it2 != rg2.end()));
    }

    StableFilteredComplex* st_;
  };

  std::vector<Simplex_handle> currentSimplices_; /**< Stores current simplex handles. */
  StableFilteredComplex* complex_;               /**< Pointer to the complex. */
  std::size_t currentSimplexIndex_;              /**< Index to current position in currentSimplices_. */
  EdgeRangeIterator currentEdgeIt_;              /**< Iterator pointing to the next edge. */
  EdgeRangeIterator endEdgeIt_;                  /**< End edge iterator. */
  bool currentDirection_;                        /**< Current direction. */
  const int maxDimension_;                       /**< Maximal dimension of expansion. */
  std::tuple<Simplex_handle, Filtration_value, bool> currentArrow_; /**< Current return element. */
  Simplex_key currentArrowNumber_;                                  /**< Number of increments. */

  /**
   * @brief Sets the iterator as the end iterator.
   */
  void _set_end() { complex_ = nullptr; }
};

/**
 * @class Oscillating_rips_simplex_iterator_range oscillating_rips_simplex_ranges.h
 * gudhi/Zigzag_persistence/oscillating_rips_simplex_ranges.h
 * @brief Gives access to a range of the ordered simplices in an oscillating Rips filtration.
 *
 * @ingroup zigzag_persistence
 *
 * @details The simplices are returned as tuples of three elements:
 * the first is the simplex handle of the simplex in the given complex,
 * the second is the filtration value of the corresponding arrow,
 * the third is the direction of the arrow, i.e., indicates if the simplex is inserted or removed.
 * For more information, see @ref Oscillating_rips_iterator.
 *
 * @tparam StableFilteredComplex Filtered complex structure that has stable simplex handles,
 * that is, they do not invalidates after an insertion or a removal (except for the removed simplices).
 * @tparam EdgeRangeIterator Type of the edge range iterator. Each element of the range is assumed to have
 * type @ref Zigzag_edge or at least to be a class with the same public methods than @ref Zigzag_edge.
 *
 * @warning As the custom iterator stores a possibly large range, avoid copying it.
 * Use @ref get_iterator_range, @ref begin and @ref end wisely.
 */
template <class StableFilteredComplex, typename EdgeRangeIterator>
class Oscillating_rips_simplex_iterator_range
{
 public:
  /**
   * @class Oscillating_rips_simplex_iterator oscillating_rips_simplex_ranges.h \
   * gudhi/Zigzag_persistence/oscillating_rips_simplex_ranges.h
   * @brief Custom iterator over the edges of an oscillating rips filtration.
   *
   * It inherits from boost::iterator_facade.
   *
   * The iterator stores the distance matrix and epsilon values to compute the current edge on the fly
   * at each increment.
   *
   * @warning As the iterator stores possibly large ranges, avoid copying it.
   */
  class Oscillating_rips_simplex_iterator
      : public boost::iterator_facade<Oscillating_rips_simplex_iterator,
                                      const std::tuple<typename StableFilteredComplex::Simplex_handle,
                                                       typename StableFilteredComplex::Filtration_value,
                                                       bool>&,
                                      boost::forward_traversal_tag>
  {
   public:
    using Filtration_value = typename StableFilteredComplex::Filtration_value; /**< Filtration value type. */
    using Simplex_handle = typename StableFilteredComplex::Simplex_handle;     /**< Simplex handle type. */

    /**
     * @brief Constructor. Takes already computed epsilon values as input, but assumes that the points are
     * already ordered accordingly. Assumes also that the last epsilon value is 0.
     *
     * @tparam PointRange Point range type.
     * @tparam DistanceFunction Type of the distance function.
     * @param nu Lower multiplier.
     * @param mu Upper multiplier.
     * @param orderedPoints Point cloud as an ordered range. The format of a point has to correspond to the input
     * format of the distance function. The order of the points should be in correspondence with the order of
     * the epsilon values.
     * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
     * and return the distance between those points.
     * @param epsilonValues Epsilon values for the oscillating rips filtration. Should be in decreasing order.
     * And the last value should be 0.
     */
    Oscillating_rips_simplex_iterator(
        Oscillating_rips_simplex_iterator_base<StableFilteredComplex, EdgeRangeIterator>* base)
        : base_iterator_(base)
    {
    }

   private:
    // mandatory for the boost::iterator_facade inheritance.
    friend class boost::iterator_core_access;

    /**
     * @brief Pointer to heavy iterator, to avoid copying it.
     */
    std::shared_ptr<Oscillating_rips_simplex_iterator_base<StableFilteredComplex, EdgeRangeIterator> > base_iterator_;

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Indicates if to iterators are equal.
     *
     * @param other Iterator to compare.
     * @return True, the iterators are pointing to the same position.
     * @return False, otherwise.
     */
    bool equal(Oscillating_rips_simplex_iterator const& other) const
    {
      return base_iterator_->equal(*other.base_iterator_);
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Returns the value of the dereferenced iterator.
     *
     * @return Current edge.
     */
    const std::tuple<Simplex_handle, Filtration_value, bool>& dereference() const
    {
      return base_iterator_->dereference();
    }

    /**
     * @brief Mandatory for the boost::iterator_facade inheritance. Increments the iterator.
     */
    void increment() { base_iterator_->increment(); }
  };

  Oscillating_rips_simplex_iterator_range(const EdgeRangeIterator& edgeStartIterator,
                                          const EdgeRangeIterator& edgeEndIterator,
                                          StableFilteredComplex& complex,
                                          int maxDimension = -1)
      : edgeStartIterator_(edgeStartIterator),
        edgeEndIterator_(edgeEndIterator),
        complex_(&complex),
        maxDimension_(maxDimension)
  {
  }

  /**
   * @brief Returns the begin iterator of a the range of edges based on @ref Oscillating_rips_simplex_iterator.
   * See the @ref zigzagrips "introduction page" for more details about the arguments.
   *
   * @tparam PointRange Point range type.
   * @tparam DistanceFunction Type of the distance function.
   * @param nu Lower multiplier.
   * @param mu Upper multiplier.
   * @param points Point cloud as a range.The format of a point has to correspond to the input format of the
   * distance function.
   * @param distance Distance function. Has to take two points as it from the range @p points as input parameters
   * and return the distance between those points.
   * @param orderPolicy Order policy for the points. Can be either
   * @ref Oscillating_rips_edge_order_policy::FARTHEST_POINT_ORDERING,
   * @ref Oscillating_rips_edge_order_policy::ALREADY_ORDERED or
   * @ref Oscillating_rips_edge_order_policy::RANDOM_POINT_ORDERING.
   *
   * @return Instantiation of @ref Oscillating_rips_simplex_iterator.
   *
   * @warning Avoid copying the iterator as it is heavier than usual iterators.
   */
  Oscillating_rips_simplex_iterator begin()
  {
    // shared pointer on the other side will take ownership
    // enables begin() to be called several times without invalidating precedent iterators
    // still have the inconvenience that all copies of a same iterator (sharing the same base) increment simultaneously
    return Oscillating_rips_simplex_iterator(
        new Oscillating_rips_simplex_iterator_base<StableFilteredComplex, EdgeRangeIterator>(
            edgeStartIterator_, edgeEndIterator_, complex_, maxDimension_));
  }

  /**
   * @brief Returns the end iterator of a the range of edges based on @ref Oscillating_rips_simplex_iterator.
   *
   * @return Default instantiation of @ref Oscillating_rips_simplex_iterator.
   */
  Oscillating_rips_simplex_iterator end() { return endIt_; }

 private:
  EdgeRangeIterator edgeStartIterator_;
  EdgeRangeIterator edgeEndIterator_;
  StableFilteredComplex* complex_;
  int maxDimension_;
  inline static const Oscillating_rips_simplex_iterator endIt_ = Oscillating_rips_simplex_iterator(
      new Oscillating_rips_simplex_iterator_base<StableFilteredComplex, EdgeRangeIterator>());
};

template <class StableFilteredComplex>
class Oscillating_rips_simplex_vector_range_constructor
{
 public:
  using Vertex_handle = typename StableFilteredComplex::Vertex_handle;
  using Simplex_handle = typename StableFilteredComplex::Simplex_handle;
  using Filtration_value = typename StableFilteredComplex::Filtration_value; /**< Filtration value type. */

  struct Simplex {
    std::vector<Vertex_handle> vertices;
    Filtration_value filtration_value;
    bool direction;
  };

  template <typename EdgeRangeIterator>
  static std::vector<Simplex> make_range(EdgeRangeIterator edgeStartIterator,
                                         EdgeRangeIterator edgeEndIterator,
                                         int maxDimension = -1)
  {
    StableFilteredComplex complex;
    std::vector<Simplex> simplices;
    std::vector<Simplex_handle> currentSimplices;
    unsigned int currentArrowNumber = 0;

    while (edgeStartIterator != edgeEndIterator) {
      const auto& edge = *edgeStartIterator;
      currentSimplices.clear();
      Filtration_value fil = edge.get_filtration_value();
      if (edge.get_direction()) {
        Oscillating_rips_edge_expander::expand_edges_of_same_filtration_value(
            edgeStartIterator, edgeEndIterator, maxDimension, complex, currentSimplices);
        unsigned int i = simplices.size();
        simplices.resize(i + currentSimplices.size());
        for (Simplex_handle sh : currentSimplices) {
          Simplex& s = simplices[i++];
          for (auto v : complex.simplex_vertex_range(sh)) s.vertices.push_back(v);
          std::reverse(s.vertices.begin(), s.vertices.end());
          s.filtration_value = fil;
          s.direction = true;
          complex.assign_key(sh, currentArrowNumber++);
        }
      } else {
        Oscillating_rips_edge_expander::collapse_edges_of_same_filtration_value(
            edgeStartIterator, edgeEndIterator, complex, currentSimplices);
        unsigned int i = simplices.size();
        simplices.resize(i + currentSimplices.size());
        for (Simplex_handle sh : currentSimplices) {
          Simplex& s = simplices[i++];
          for (auto v : complex.simplex_vertex_range(sh)) s.vertices.push_back(v);
          std::reverse(s.vertices.begin(), s.vertices.end());
          s.filtration_value = fil;
          s.direction = false;
          complex.remove_maximal_simplex(sh);
        }
      }
    }

    return simplices;
  }
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_OSCILLATING_RIPS_SIMPLEX_RANGES_H_
