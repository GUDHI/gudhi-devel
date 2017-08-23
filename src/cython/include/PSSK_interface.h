/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University
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

#ifndef PSSK_INTERFACE_H_
#define PSSK_INTERFACE__H_


#include <gudhi/PSSK.h>


namespace Gudhi {
namespace Persistence_representations {

/**
* This is a version of a representation presented in https://arxiv.org/abs/1412.6821
* In that paper the authors are using the representation just to compute kernel. Over here, we extend the usability by
*far.
* Note that the version presented here is not exact, since we are discretizing the kernel.
* The only difference with respect to the original class is the method of creation. We have full (square) image, and for
*every point (p,q), we add a kernel at (p,q) and the negative kernel
* at (q,p)
**/

class PSSK_interface : public PSSK {
 public:
  PSSK_interface(){}

  PSSK_interface(const std::vector<std::pair<double, double> >& interval,
       std::vector<std::vector<double> > filter = create_Gaussian_filter(5, 1), size_t number_of_pixels = 1000,
       double min_ = -1, double max_ = -1)
      :
      PSSK(interval,filter,number_of_pixels,min_,max_){}

  PSSK_interface(const char* filename, std::vector<std::vector<double> > filter = create_Gaussian_filter(5, 1),
       size_t number_of_pixels = 1000, double min_ = -1, double max_ = -1,
       unsigned dimension = std::numeric_limits<unsigned>::max())
      :PSSK(filename,filter,number_of_pixels,min_,max_,dimension){}
 
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PSSK_INTERFACE__H_
