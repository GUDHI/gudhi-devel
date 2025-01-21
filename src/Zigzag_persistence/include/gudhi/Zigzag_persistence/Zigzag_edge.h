/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria and Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Zigzag_edge.h
 * @author Clément Maria, Hannah Schreiber
 * @brief Contains the implementation of the @ref Gudhi::zigzag_persistence::Zigzag_edge class.
 */

#ifndef ZIGZAG_ZZ_EDGE_H_
#define ZIGZAG_ZZ_EDGE_H_

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <ostream>
// #include <iomanip>

// #include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace zigzag_persistence {

/**
 * @class Zigzag_edge Zigzag_edge.h gudhi/Zigzag_persistence/Zigzag_edge.h
 * @brief Edge structure for the oscillating Rips filtration.
 *
 * @ingroup zigzag_persistence
 *
 * @details The edges of the filtration are computed first and the remaining simplices are deduced from them.
 * Is also used to represent the vertices for technical reasons, by giving both vertices the same value.
 *
 * @tparam Filtration_value Type of the filtration value. Is recommended to be easy to copy, like an arithmetic type.
 */
template <typename Filtration_value>
class Zigzag_edge
{
 public:
  /**
   * @brief Constructor.
   *
   * @param u First boundary vertex ID
   * @param v Second boundary vertex ID. If same than @p u, the edge is considered a vertex.
   * @param fil Filtration value
   * @param direction If true, forwards. If false, backwards.
   */
  Zigzag_edge(int u, int v, Filtration_value fil, bool direction) : u_(u), v_(v), fil_(fil), direction_(direction)
  {
    if (u > v) std::swap(u_, v_);
  }

  /**
   * @brief Default constructor. Initialize everything to zero and the direction to true.
   */
  Zigzag_edge() : u_(0), v_(0), fil_(0), direction_(true) {}

  /**
   * @brief Returns vertex with smaller ID.
   *
   * @return Smallest ID among the boundary vertices.
   */
  int get_smallest_vertex() const { return u_; }

  /**
   * @brief Returns vertex with bigger ID.
   *
   * @return Biggest ID among the boundary vertices.
   */
  int get_biggest_vertex() const { return v_; }

  /**
   * @brief Returns the filtration value of the edge.
   *
   * @return Filtration value of the edge.
   */
  Filtration_value get_filtration_value() const { return fil_; }

  /**
   * @brief Gives the direction of the arrow corresponding to the edge.
   *
   * @return True, if forward, i.e., an insertion.
   * @return False, if backward, i.e., a removal.
   */
  bool get_direction() const { return direction_; }

  void set(int u, int v, Filtration_value fil, bool direction)
  {
    u_ = u;
    v_ = v;
    fil_ = fil;
    direction_ = direction;
  }

  /**
   * @brief Equality test
   *
   * @param e Edge to compare.
   * @return True, if both edges are equal.
   * @return False, if both edges are not equal.
   */
  bool operator==(const Zigzag_edge& e) const
  {
    return ((e.u_ == u_) && (e.v_ == v_) && (e.fil_ == fil_) && (e.direction_ == direction_));
  }

  friend std::ostream& operator<<(std::ostream& stream, const Zigzag_edge& ze)
  {
    // stream << std::setprecision(6);
    // stream << std::setprecision(std::numeric_limits<Filtration_value>::digits);
    stream << "(" << ze.u_ << ", " << ze.v_ << ") ";
    if (ze.direction_) {
      stream << "-- " << ze.fil_ << " -->";
    } else {
      stream << "<-- " << ze.fil_ << " --";
    }
    return stream;
  }

 private:
  int u_;                /**< Smaller vertex. */
  int v_;                /**< Bigger vertex. */
  Filtration_value fil_; /**< Filtration value. */
  bool direction_;       /**< Direction. True = forward, false = backward. */
};

}  // namespace zigzag_persistence
}  // namespace Gudhi

#endif  // ZIGZAG_ZZ_EDGE_H_
