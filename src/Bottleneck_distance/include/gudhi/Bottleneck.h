/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       Francois Godi
 *
 *    Copyright (C) 2015  INRIA (France)
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

#ifndef BOTTLENECK_H_
#define BOTTLENECK_H_

#include <gudhi/Graph_matching.h>
#include <cmath>

namespace Gudhi {

namespace persistence_diagram {

/** \brief Function to use in order to compute the Bottleneck distance between two persistence diagrams (see Concepts).
 * If the last parameter e is not 0 (default value if not explicited), you get an additive e-approximation.
 *
 * \ingroup bottleneck_distance
 */
template<typename Persistence_diagram1, typename Persistence_diagram2>
double bottleneck_distance(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2, double e=0.) {
    Persistence_graph g(diag1, diag2, e);
    double b = g.bottleneck_alive();
    if(b == std::numeric_limits<double>::infinity())
        return std::numeric_limits<double>::infinity();
    std::vector<double> sd;
    if(e == 0.)
        sd = g.sorted_distances();
    long idmin = 0;
    long idmax = e==0. ? sd.size() - 1 : g.diameter_bound()/e + 1;
    // alpha can change the complexity
    double alpha = std::pow(idmax, 0.25);
    Graph_matching m(g);
    Graph_matching biggest_unperfect(g);
    while (idmin != idmax) {
        long step = static_cast<long>((idmax - idmin) / alpha);
        m.set_r(e == 0. ?  sd.at(idmin + step) : e*(idmin + step));
        while (m.multi_augment());
        //The above while compute a maximum matching (according to the r setted before)
        if (m.perfect()) {
            idmax = idmin + step;
            m = biggest_unperfect;
        } else {
            biggest_unperfect = m;
            idmin = idmin + step + 1;
        }
    }
    b = std::max(b, e == 0. ? sd.at(idmin) : e*(idmin));
    return b;
}



}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // BOTTLENECK_H_
