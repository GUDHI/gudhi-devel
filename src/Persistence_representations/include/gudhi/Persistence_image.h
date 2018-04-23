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

#ifndef PERSISTENCE_IMAGE_H_
#define PERSISTENCE_IMAGE_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>
#include <gudhi/Weight_functions.h>
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
 * \class Persistence_image gudhi/Persistence_image.h
 * \brief A class implementing the persistence images.
 *
 * \ingroup Persistence_representations
 *
 * \details
 *
 * Persistence images are a way to build images from persistence diagrams. Roughly, the idea is to center Gaussians on each diagram point, with a weight that usually depends on
 * the distance to the diagonal, so that the diagram is turned into a function, and then to discretize the plane into pixels, and integrate this function on each pixel.
 * See \cite Persistence_Images_2017 for more details.
 *
**/

class Persistence_image {

 protected:
  Persistence_diagram diagram;
  int res_x, res_y;
  double min_x, max_x, min_y, max_y;
  Weight weight;
  double sigma;

 public:

  /** \brief Persistence Image constructor.
   * \ingroup Persistence_image
   *
   * @param[in] _diagram   persistence diagram.
   * @param[in] _min_x     minimum value of pixel abscissa.
   * @param[in] _max_x     maximum value of pixel abscissa.
   * @param[in] _res_x     number of pixels for the x-direction.
   * @param[in] _min_y     minimum value of pixel ordinate.
   * @param[in] _max_y     maximum value of pixel ordinate.
   * @param[in] _res_y     number of pixels for the y-direction.
   * @param[in] _weight    weight function for the Gaussians.
   * @param[in] _sigma     bandwidth parameter for the Gaussians.
   *
   */
  Persistence_image(const Persistence_diagram & _diagram, double _min_x = 0.0, double _max_x = 1.0, int _res_x = 10, double _min_y = 0.0, double _max_y = 1.0, int _res_y = 10, const Weight & _weight = arctan_weight(1,1), double _sigma = 1.0){
      diagram = _diagram; min_x = _min_x; max_x = _max_x; res_x = _res_x; min_y = _min_y; max_y = _max_y; res_y = _res_y, weight = _weight; sigma = _sigma;
  }

  /** \brief Computes the persistence image of a diagram.
   * \ingroup Persistence_image
   *
   */
  std::vector<std::vector<double> > vectorize() const {
    std::vector<std::vector<double> > im; for(int i = 0; i < res_y; i++)  im.emplace_back();
    double step_x = (max_x - min_x)/res_x; double step_y = (max_y - min_y)/res_y;

    int num_pts = diagram.size();

    for(int i = 0; i < res_y; i++){
      double y = min_y + i*step_y;
      for(int j = 0; j < res_x; j++){
        double x = min_x + j*step_x;

        double pixel_value = 0;
        for(int k = 0; k < num_pts; k++){
          double px = diagram[k].first; double py = diagram[k].second;
          pixel_value += weight(std::pair<double,double>(px,py)) * std::exp(   -((x-px)*(x-px) + (y-(py-px))*(y-(py-px))) / (2*sigma*sigma)   ) / (sigma*std::sqrt(2*pi));
        }
        im[i].push_back(pixel_value);

      }
    }

    return im;

  }




}; // class Persistence_image
}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_IMAGE_H_
