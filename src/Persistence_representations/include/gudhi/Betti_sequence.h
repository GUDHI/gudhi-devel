/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carriere
 *
 *    Copyright (C) 2018  INRIA (France)
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

#ifndef BETTI_SEQUENCE_H_
#define BETTI_SEQUENCE_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Debug_utils.h>

// standard include
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>
#include <functional>

namespace Gudhi {
namespace Persistence_representations {

/**
 * \class Betti_sequence gudhi/Betti_sequence.h
 * \brief A class implementing Betti sequences
 *
 * \ingroup Persistence_representations
 *
 * \details
**/

class Betti_sequence {

  protected:
   Persistence_diagram diagram;
   int res_x, nb_cv;
   double min_x, max_x;

  public:

   /** \brief Betti_sequence constructor.
    * \ingroup Betti_sequence
    *
    * @param[in] _diagram    persistence diagram.
    * @param[in] _min_x      minimum value of samples.
    * @param[in] _max_x      maximum value of samples.
    * @param[in] _res_x      number of samples.
    *
    */
   Betti_sequence(const Persistence_diagram & _diagram, double _min_x = 0.0, double _max_x = 1.0, int _res_x = 10){diagram = _diagram; min_x = _min_x; max_x = _max_x; res_x = _res_x;}

   /** \brief Computes the Betti sequences of a diagram.
    * \ingroup Betti_sequence
    *
    */
   std::vector<int> vectorize() const {
     int num_pts = diagram.size(); double step = (max_x - min_x)/(res_x - 1);
     std::vector<int> bs(res_x); for(int i = 0; i < res_x; i++)  bs[i] = 0;
     for(int j = 0; j < num_pts; j++){
       double px = diagram[j].first; double py = diagram[j].second;
       int first = std::ceil((px-min_x)/step); int last = std::ceil((py-min_x)/step);
       for(int i = first; i < last; i++)  bs[i] += 1;
     }

     return bs;
   }

}; // class Betti_sequence
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // BETTI_SEQUENCE_H_
