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
 * \class Persistence_landscape_on_grid_exact gudhi/Persistence_landscape_on_grid_exact.h
 * \brief A class implementing exact persistence landscapes by approximating them on a collection of grid points
 *
 * \ingroup Persistence_representations
 *
 * \details
 * In this class, we propose a way to approximate landscapes by sampling the x-axis of the persistence diagram and evaluating the (exact) landscape functions on the sample projections onto the diagonal. Note that this is a different approximation scheme 
 * from the one proposed in Persistence_landscape_on_grid, which puts a grid on the diagonal, maps the persistence intervals on this grid and computes the (approximate) landscape functions on the samples. 
 * Hence, the difference is that we do not modify the diagram in this implementation, but the code can be slower to run.
**/

class Persistence_landscape_on_grid_exact {

  protected:
   Persistence_diagram diagram;
   int res_x, nb_ls;
   double min_x, max_x;

  public:

   /** \brief Persistence_landscape_on_grid_exact constructor.
    * \ingroup Persistence_landscape_on_grid_exact
    *
    * @param[in] _diagram    persistence diagram.
    * @param[in] _nb_ls      number of landscape functions.
    * @param[in] _min_x      minimum value of samples.
    * @param[in] _max_x      maximum value of samples.
    * @param[in] _res_x      number of samples.
    *
    */
   Persistence_landscape_on_grid_exact(const Persistence_diagram & _diagram, int _nb_ls = 5, double _min_x = 0.0, double _max_x = 1.0, int _res_x = 10){diagram = _diagram; nb_ls = _nb_ls; min_x = _min_x; max_x = _max_x; res_x = _res_x;}

   /** \brief Computes the landscape approximation of a diagram.
    * \ingroup Persistence_landscape_on_grid_exact
    *
    */
   std::vector<std::vector<double> > vectorize() const {
     std::vector<std::vector<double> >  ls; for(int i = 0; i < nb_ls; i++)  ls.emplace_back();
     int num_pts = diagram.size(); double step = (max_x - min_x)/(res_x - 1);

     std::vector<std::vector<double> > ls_t; for(int i = 0; i < res_x; i++)  ls_t.emplace_back();
     for(int j = 0; j < num_pts; j++){
       double px = diagram[j].first; double py = diagram[j].second; double mid = (px+py)/2;
       int first = std::ceil((px-min_x)/step); int middle = std::ceil((mid-min_x)/step); int last = std::ceil((py-min_x)/step); double x = min_x + first*step;
       for(int i = first; i < middle; i++){  double value = std::sqrt(2)*(x-px); ls_t[i].push_back(value); x += step;  }
       for(int i = middle; i < last; i++){   double value = std::sqrt(2)*(py-x); ls_t[i].push_back(value); x += step;  }
     }

     for(int i = 0; i < res_x; i++){
       std::sort(ls_t[i].begin(), ls_t[i].end(), [](const double & a, const double & b){return a > b;});
       int nb_events_i = ls_t[i].size();
       for (int j = 0; j < nb_ls; j++){  if(j < nb_events_i)  ls[j].push_back(ls_t[i][j]);  else  ls[j].push_back(0);  }
     }

     return ls;
   }

}; // class Persistence_landscape_on_grid_exact
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // LANDSCAPE_H_
