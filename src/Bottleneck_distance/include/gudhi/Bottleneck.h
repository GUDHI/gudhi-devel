/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - 2019/06 Vincent Rouvreau : Fix doxygen warning.
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef BOTTLENECK_H_
#define BOTTLENECK_H_

#include <gudhi/Graph_matching.h>

#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <vector>
#include <algorithm>  // for max
#include <limits>  // for numeric_limits

#include <cmath>
#include <cfloat>  // FLT_EVAL_METHOD

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error bottleneck_distance is only available for CGAL >= 4.11
#endif

namespace Gudhi {

namespace persistence_diagram {

inline double bottleneck_distance_approx(Persistence_graph& g, double e) {
  double b_lower_bound = 0.;
  double b_upper_bound = g.max_dist_to_diagonal();
  int siz = g.size();
  if (siz <= 1)
    // The value of alpha would be wrong in this case
    return b_upper_bound;
  const double alpha = std::pow(siz, 1. / 5.);
  Graph_matching m(g);
  Graph_matching biggest_unperfect(g);
  while (b_upper_bound - b_lower_bound > 2 * e) {
    double step = b_lower_bound + (b_upper_bound - b_lower_bound) / alpha;
#if !defined FLT_EVAL_METHOD || FLT_EVAL_METHOD < 0 || FLT_EVAL_METHOD > 1
    // On platforms where double computation is done with excess precision,
    // we force it to its true precision so the following test is reliable.
    volatile double drop_excess_precision = step;
    step = drop_excess_precision;
    // Alternative: step = CGAL::IA_force_to_double(step);
#endif
    if (step <= b_lower_bound || step >= b_upper_bound)  // Avoid precision problem
      break;
    m.set_r(step);
    while (m.multi_augment()) {}  // compute a maximum matching (in the graph corresponding to the current r)
    if (m.perfect()) {
      m = biggest_unperfect;
      b_upper_bound = step;
    } else {
      biggest_unperfect = m;
      b_lower_bound = step;
    }
  }
  return (b_lower_bound + b_upper_bound) / 2.;
}

inline double bottleneck_distance_exact(Persistence_graph& g) {
  std::vector<double> sd = g.sorted_distances();
  long lower_bound_i = 0;
  long upper_bound_i = sd.size() - 1;
  const double alpha = std::pow(g.size(), 1. / 5.);
  Graph_matching m(g);
  Graph_matching biggest_unperfect(g);
  while (lower_bound_i != upper_bound_i) {
    long step = lower_bound_i + static_cast<long> ((upper_bound_i - lower_bound_i - 1) / alpha);
    m.set_r(sd.at(step));
    while (m.multi_augment()) {}  // compute a maximum matching (in the graph corresponding to the current r)
    if (m.perfect()) {
      m = biggest_unperfect;
      upper_bound_i = step;
    } else {
      biggest_unperfect = m;
      lower_bound_i = step + 1;
    }
  }
  return sd.at(lower_bound_i);
}

/** \brief Function to compute the Bottleneck distance between two persistence diagrams.
 *
 * \tparam Persistence_diagram1,Persistence_diagram2
 * models of the concept `PersistenceDiagram`.
 *
 * \param[in] diag1 The first persistence diagram.
 * \param[in] diag2 The second persistence diagram.
 *
 * \param[in] e
 * \parblock
 * If `e` is 0, this uses an expensive algorithm to compute the exact distance.
 *
 * If `e` is not 0, it asks for an additive `e`-approximation, and currently
 * also allows a small multiplicative error (the last 2 or 3 bits of the
 * mantissa may be wrong). This version of the algorithm takes advantage of the
 * limited precision of `double` and is usually a lot faster to compute,
 * whatever the value of `e`.
 *
 * Thus, by default, `e` is the smallest positive double.
 * \endparblock
 *
 * \ingroup bottleneck_distance
 */
template<typename Persistence_diagram1, typename Persistence_diagram2>
double bottleneck_distance(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2,
                           double e = (std::numeric_limits<double>::min)()) {
  Persistence_graph g(diag1, diag2, e);
  if (g.bottleneck_alive() == std::numeric_limits<double>::infinity())
    return std::numeric_limits<double>::infinity();
  return (std::max)(g.bottleneck_alive(), e == 0. ? bottleneck_distance_exact(g) : bottleneck_distance_approx(g, e));
}

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // BOTTLENECK_H_
