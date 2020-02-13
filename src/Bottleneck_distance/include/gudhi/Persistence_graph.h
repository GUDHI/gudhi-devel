/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_GRAPH_H_
#define PERSISTENCE_GRAPH_H_

#include <gudhi/Internal_point.h>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

#include <vector>
#include <algorithm>
#include <limits>  // for numeric_limits

namespace Gudhi {

namespace persistence_diagram {

/** \internal \brief Structure representing a Euclidean bipartite graph containing
 *  the points from the two persistence diagrams (including the projections).
 *
 * \ingroup bottleneck_distance
 */
class Persistence_graph {
 public:
  /** \internal \brief Constructor taking 2 PersistenceDiagrams (concept) as parameters. */
  template<typename Persistence_diagram1, typename Persistence_diagram2>
  Persistence_graph(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e);
  /** \internal \brief Is the given point from U the projection of a point in V ? */
  bool on_the_u_diagonal(int u_point_index) const;
  /** \internal \brief Is the given point from V the projection of a point in U ? */
  bool on_the_v_diagonal(int v_point_index) const;
  /** \internal \brief Given a point from V, returns the corresponding (projection or projector) point in U. */
  int corresponding_point_in_u(int v_point_index) const;
  /** \internal \brief Given a point from U, returns the corresponding (projection or projector) point in V. */
  int corresponding_point_in_v(int u_point_index) const;
  /** \internal \brief Given a point from U and a point from V, returns the distance between those points. */
  double distance(int u_point_index, int v_point_index) const;
  /** \internal \brief Returns size = |U| = |V|. */
  int size() const;
  /** \internal \brief Is there as many infinite points (alive components) in both diagrams ? */
  double bottleneck_alive() const;
  /** \internal \brief Returns the O(n^2) sorted distances between the points. */
  std::vector<double> sorted_distances() const;
  /** \internal \brief Returns an upper bound for the diameter of the convex hull of all non infinite points */
  double diameter_bound() const;
  /** \internal \brief Returns the corresponding internal point */
  Internal_point get_u_point(int u_point_index) const;
  /** \internal \brief Returns the corresponding internal point */
  Internal_point get_v_point(int v_point_index) const;

 private:
  std::vector<Internal_point> u;
  std::vector<Internal_point> v;
  double b_alive;
};

template<typename Persistence_diagram1, typename Persistence_diagram2>
Persistence_graph::Persistence_graph(const Persistence_diagram1 &diag1,
                                     const Persistence_diagram2 &diag2, double e)
    : u(), v(), b_alive(0.) {
  std::vector<double> u_alive;
  std::vector<double> v_alive;
  for (auto it = std::begin(diag1); it != std::end(diag1); ++it) {
    if (std::get<1>(*it) == std::numeric_limits<double>::infinity())
      u_alive.push_back(std::get<0>(*it));
    else if (std::get<1>(*it) - std::get<0>(*it) > e)
      u.push_back(Internal_point(std::get<0>(*it), std::get<1>(*it), u.size()));
  }
  for (auto it = std::begin(diag2); it != std::end(diag2); ++it) {
    if (std::get<1>(*it) == std::numeric_limits<double>::infinity())
      v_alive.push_back(std::get<0>(*it));
    else if (std::get<1>(*it) - std::get<0>(*it) > e)
      v.push_back(Internal_point(std::get<0>(*it), std::get<1>(*it), v.size()));
  }
  if (u.size() < v.size())
    swap(u, v);
  std::sort(u_alive.begin(), u_alive.end());
  std::sort(v_alive.begin(), v_alive.end());
  if (u_alive.size() != v_alive.size()) {
    b_alive = std::numeric_limits<double>::infinity();
  } else {
    for (auto it_u = u_alive.cbegin(), it_v = v_alive.cbegin(); it_u != u_alive.cend(); ++it_u, ++it_v)
      b_alive = (std::max)(b_alive, std::fabs(*it_u - *it_v));
  }
}

inline bool Persistence_graph::on_the_u_diagonal(int u_point_index) const {
  return u_point_index >= static_cast<int> (u.size());
}

inline bool Persistence_graph::on_the_v_diagonal(int v_point_index) const {
  return v_point_index >= static_cast<int> (v.size());
}

inline int Persistence_graph::corresponding_point_in_u(int v_point_index) const {
  return on_the_v_diagonal(v_point_index) ?
      v_point_index - static_cast<int> (v.size()) : v_point_index + static_cast<int> (u.size());
}

inline int Persistence_graph::corresponding_point_in_v(int u_point_index) const {
  return on_the_u_diagonal(u_point_index) ?
      u_point_index - static_cast<int> (u.size()) : u_point_index + static_cast<int> (v.size());
}

inline double Persistence_graph::distance(int u_point_index, int v_point_index) const {
  if (on_the_u_diagonal(u_point_index) && on_the_v_diagonal(v_point_index))
    return 0.;
  Internal_point p_u = get_u_point(u_point_index);
  Internal_point p_v = get_v_point(v_point_index);
  return (std::max)(std::fabs(p_u.x() - p_v.x()), std::fabs(p_u.y() - p_v.y()));
}

inline int Persistence_graph::size() const {
  return static_cast<int> (u.size() + v.size());
}

inline double Persistence_graph::bottleneck_alive() const {
  return b_alive;
}

inline std::vector<double> Persistence_graph::sorted_distances() const {
  std::vector<double> distances;
  distances.push_back(0.);  // for empty diagrams
  for (int u_point_index = 0; u_point_index < size(); ++u_point_index) {
    distances.push_back(distance(u_point_index, corresponding_point_in_v(u_point_index)));
    for (int v_point_index = 0; v_point_index < size(); ++v_point_index)
      distances.push_back(distance(u_point_index, v_point_index));
  }
#ifdef GUDHI_USE_TBB
  tbb::parallel_sort(distances.begin(), distances.end());
#else
  std::sort(distances.begin(), distances.end());
#endif
  return distances;
}

inline Internal_point Persistence_graph::get_u_point(int u_point_index) const {
  if (!on_the_u_diagonal(u_point_index))
    return u.at(u_point_index);
  Internal_point projector = v.at(corresponding_point_in_v(u_point_index));
  double m = (projector.x() + projector.y()) / 2.;
  return Internal_point(m, m, u_point_index);
}

inline Internal_point Persistence_graph::get_v_point(int v_point_index) const {
  if (!on_the_v_diagonal(v_point_index))
    return v.at(v_point_index);
  Internal_point projector = u.at(corresponding_point_in_u(v_point_index));
  double m = (projector.x() + projector.y()) / 2.;
  return Internal_point(m, m, v_point_index);
}

inline double Persistence_graph::diameter_bound() const {
  double max = 0.;
  for (auto it = u.cbegin(); it != u.cend(); it++)
    max = (std::max)(max, it->y());
  for (auto it = v.cbegin(); it != v.cend(); it++)
    max = (std::max)(max, it->y());
  return max;
}

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // PERSISTENCE_GRAPH_H_
