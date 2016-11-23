/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Francois Godi
 *
 *    Copyright (C) 2015  INRIA Sophia-Antipolis (France)
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

namespace Gudhi {

namespace bottleneck_distance {

/** \brief Function to use in order to compute the Bottleneck distance between two persistence diagrams. You get an additive e-approximation.
 *
 *
 * \ingroup bottleneck_distance
 */
template<typename Persistence_diagram1, typename Persistence_diagram2>
double compute(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2, double e = 0.);

template<typename Persistence_diagram1, typename Persistence_diagram2>
double compute_exactly(const Persistence_diagram1& diag1, const Persistence_diagram2& diag2);

template<typename Persistence_diagram1, typename Persistence_diagram2>
double compute_exactly(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2) {
    G::initialize(diag1, diag2, 0.);
    std::shared_ptr< std::vector<double> > sd(G::sorted_distances());
    int idmin = 0;
    int idmax = sd->size() - 1;
    // alpha can be modified, this will change the complexity
    double alpha = pow(sd->size(), 0.25);
    Graph_matching m;
    Graph_matching biggest_unperfect;
    while (idmin != idmax) {
        int step = static_cast<int>((idmax - idmin) / alpha);
        m.set_r(sd->at(idmin + step));
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
    return sd->at(idmin);
}

template<typename Persistence_diagram1, typename Persistence_diagram2>
double compute(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2, double e) {
    if(e< std::numeric_limits<double>::min())
        return compute_exactly(diag1, diag2);
    G::initialize(diag1, diag2, e);
    int in = G::diameter()/e;
    int idmin = 0;
    int idmax = in;
    // alpha can be modified, this will change the complexity
    double alpha = pow(in, 0.25);
    Graph_matching m;
    Graph_matching biggest_unperfect;
    while (idmin != idmax) {
        int step = static_cast<int>((idmax - idmin) / alpha);
        m.set_r(e*(idmin + step));
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
    return e*(idmin);
}



}  // namespace bottleneck_distance

}  // namespace Gudhi

#endif  // BOTTLENECK_H_

/* Dichotomic version
template<typename Persistence_diagram1, typename Persistence_diagram2>
double compute(const Persistence_diagram1 &diag1, const Persistence_diagram2 &diag2, double e) {
    if(e< std::numeric_limits<double>::min())
        return compute_exactly(diag1, diag2);
    G::initialize(diag1, diag2, e);
    double d = 0.;
    double f = G::diameter();
    while (f-d > e){
        Graph_matching m;
        m.set_r((d+f)/2.);
        while (m.multi_augment());
        //The above while compute a maximum matching (according to the r setted before)
        if (m.perfect())
            f = (d+f)/2.;
         else
            d= (d+f)/2.;
    }
    return d;
} */

