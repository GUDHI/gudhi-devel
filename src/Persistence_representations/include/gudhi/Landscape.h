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

#ifndef LANDSCAPE_H_
#define LANDSCAPE_H_

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
 * \class Landscape gudhi/Landscape.h
 * \brief A class implementing landscapes.
 *
 * \ingroup Persistence_representations
 *
 * \details
 *
 * The landscape is a way to turn a persistence diagram into \f$L^2\f$ functions. Roughly, the idea is to see the boundaries of the rank functions as scalar functions taking values on the diagonal.
 * See \cite bubenik_landscapes_2015 for more details. Here we provide a way to approximate such functions by computing their values on a set of samples.
 *
**/

class Landscape {

  protected:
   Persistence_diagram diagram;
   int res_x, nb_ls;
   double min_x, max_x;

  public:

   /** \brief Landscape constructor.
    * \ingroup Landscape
    *
    * @param[in] _diagram    persistence diagram.
    * @param[in] _nb_ls      number of landscape functions.
    * @param[in] _min_x      minimum value of samples.
    * @param[in] _max_x      maximum value of samples.
    * @param[in] _res_x      number of samples.
    *
    */
   Landscape(const Persistence_diagram & _diagram, int _nb_ls = 5, double _min_x = 0.0, double _max_x = 1.0, int _res_x = 10){diagram = _diagram; nb_ls = _nb_ls; min_x = _min_x; max_x = _max_x; res_x = _res_x;}

   /** \brief Computes the landscape of a diagram.
    * \ingroup Landscape
    *
    */
   std::vector<std::vector<double> > vectorize() const {
     std::vector<std::vector<double> >  ls; for(int i = 0; i < nb_ls; i++)  ls.emplace_back();
     int num_pts = diagram.size(); double step = (max_x - min_x)/res_x;

     for(int i = 0; i < res_x; i++){
       double x = min_x + i*step; double t = x / std::sqrt(2); std::vector<double> events;
       for(int j = 0; j < num_pts; j++){
         double px = diagram[j].first; double py = diagram[j].second;
         if(t >= px && t <= py){  if(t >= (px+py)/2)  events.push_back(std::sqrt(2)*(py-t)); else  events.push_back(std::sqrt(2)*(t-px));  }
       }

       std::sort(events.begin(), events.end(), [](const double & a, const double & b){return a > b;}); int nb_events = events.size();
       for (int j = 0; j < nb_ls; j++){  if(j < nb_events)  ls[j].push_back(events[j]);  else  ls[j].push_back(0);  }
     }
     return ls;
   }




}; // class Landscape
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // LANDSCAPE_H_
