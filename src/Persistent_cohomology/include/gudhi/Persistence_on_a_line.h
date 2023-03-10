/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s): Marc Glisse
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENCE_ON_A_LINE_H_
#define PERSISTENCE_ON_A_LINE_H_

#include <vector>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <gudhi/Debug_utils.h>

namespace Gudhi::persistent_cohomology {
/**
 * \brief Computes the persistent homology of the sublevelsets of a PL function defined on \f$\mathbb{R}\f$ in linear time.
 *
 * \ingroup persistent_cohomology
 *
 * Not all input values appear in the output (they would be part of pairs of length 0
 * with a simplicial / cubical complex).
 *
 * @param[in] input Range of filtration values.
 * @param[out] out Functor that is called as `out(birth, death)` for each persistence pair.
 * By convention, it is also called one last time on the minimum and `std::numeric_limits<Filtration>::infinity()`.
 * @param[in] lt Functor that compares 2 filtration values.
 *
 * \author Marc Glisse
 */
template<class FiltrationRange, class OutputFunctor, class Compare = std::less<>>
void compute_persistence_of_function_on_line(FiltrationRange const& input, OutputFunctor&& out, Compare&& lt = {}) {
  // We process the elements of input 1 by 1, simplifying and outputting as much as possible.
  // The simplified sequence with the elements that are still active is stored in data.
  // Invariant: data contains a sequence of type 1 9 2 8 3 7 ...
  using std::begin;
  using std::end;
  auto it = begin(input);
  auto stop = end(input);
  if (it == stop) return;
  typedef std::decay_t<decltype(*it)> Filtration;
  std::vector<Filtration> data;
  Filtration v;
  auto le = [&lt](auto& x, auto& y){ return !lt(y, x); };
  auto ge = [&lt](auto& x, auto& y){ return !lt(x, y); };
  auto gt = [&lt](auto& x, auto& y){ return  lt(y, x); };
  data.push_back(*it++);
state1: // data contains a single element
  if (it == stop) goto infinite;
  v = *it++;
  if (le(v, data[0])) {
state1down:
    data[0] = v;
    goto state1;
  }
  data.push_back(v);
  goto state12;
state12: // data contains only 2 elements, necessarily data[0] < data[1]
  if (it == stop) goto endup;
  v = *it++;
  if (ge(v, data[1])) {
    data[1] = v;
    goto state12;
  }
state12down:
  if (le(v, data[0])) {
    out(data[0], data[1]);
    data.pop_back();
    data[0] = v;
    goto state1;
  }
  data.push_back(v);
  goto state132;
state132: // data[-3] < data[-1] < data[-2]
  if (it == stop) goto enddown;
  v = *it++;
  if (le(v, data.back())) {
    if (gt(v,data.end()[-3])) {
      data.back() = v; goto state132;
    } else {
      out(data.end()[-3], data.end()[-2]);
      data.erase(data.end()-3, data.end());
      if (data.empty()) { data.push_back(v); goto state1; }
      goto down;
    }
  } else {
state132up:
    if (ge(v, data.end()[-2])) {
      out(data.end()[-1], data.end()[-2]);
      data.erase(data.end()-2, data.end());
      goto up;
    } else {
      data.push_back(v);
      goto state312;
    }
  }
state312: // data[-2] < data[-1] < data[-3]
  if (it == stop) goto endup;
  v = *it++;
  if (ge(v, data.back())) {
    if (lt(v,data.end()[-3])) {
      data.back() = v; goto state312;
    } else {
      out(data.end()[-2], data.end()[-3]);
      data.erase(data.end()-3, data.end());
      GUDHI_CHECK (!data.empty(), std::logic_error("Bug in Gudhi"));
      goto up;
    }
  } else {
state312down:
    if (le(v, data.end()[-2])) {
      out(data.end()[-2], data.end()[-1]);
      data.erase(data.end()-2, data.end());
      goto down;
    } else {
      data.push_back(v);
      goto state132;
    }
  }
up: // data[-1] < v after a simplification
  if (data.size() == 1) { data.push_back(v); goto state12; }
  goto state132up;
down: // v < data[-1] after a simplification
  switch (data.size()) {
    case 1:
      goto state1down;
    case 2:
      goto state12down;
    default:
      goto state312down;
  }
  // From here on, we have finished reading input
endup: // data[-2] < data[-1]
  data.pop_back();
  goto enddown;
enddown: // data[-1] < data[-2]
  if (data.size() > 1) {
    out(data.end()[-1], data.end()[-2]);
    data.erase(data.end()-2, data.end());
    goto enddown;
  }
  goto infinite;
infinite: // data only contains the global minimum
  out(data[0], std::numeric_limits<Filtration>::infinity());
}
} // namespace Gudhi::persistent_cohomology
#endif
