/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/06 Hannah Schreiber: Various small bug fixes (missing `inline`s, `DEBUG_TRACES`s etc.)
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PSSK_H_
#define PSSK_H_

// standard include
#ifdef DEBUG_TRACES
#include <iostream> // std::clog
#endif
#include <cmath>    // std::fabs
#include <limits>   // std::numeric_limits
#include <utility>  // std::pair
#include <vector>

// gudhi include
#include <gudhi/Persistence_heat_maps.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace Persistence_representations {

// TODO: none of the methods are documented and the class description is lacking.

/**
 * This is a version of a representation presented in @cite Reininghaus_Huber_ALL_PSSK.
 * In that paper the authors are using the representation just to compute kernel. Over here, we extend the usability by
 * far.
 * Note that the version presented here is not exact, since we are discretizing the kernel.
 * The only difference with respect to the original class is the method of creation. We have full (square) image, and
 * f or every point (p,q), we add a kernel at (p,q) and the negative kernel at (q,p).
 *
 * @ingroup Persistence_representations
 **/
class PSSK : public Persistence_heat_maps<constant_scaling_function>
{
 public:
  using Base = Persistence_heat_maps<constant_scaling_function>;

  PSSK() : Base() {}

  PSSK(const std::vector<std::pair<double, double> >& interval,
       const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
       std::size_t number_of_pixels = 1000,
       double min = -1,
       double max = -1)
      : Base()
  {
    this->_construct(interval, filter, number_of_pixels, min, max);
  }

  PSSK(const char* filename,
       const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
       std::size_t number_of_pixels = 1000,
       double min = -1,
       double max = -1,
       unsigned int dimension = std::numeric_limits<unsigned int>::max())
      : Base()
  {
    std::vector<std::pair<double, double> > intervals;
    if (dimension == std::numeric_limits<unsigned int>::max()) {
      intervals = read_persistence_intervals_in_one_dimension_from_file(filename);
    } else {
      intervals = read_persistence_intervals_in_one_dimension_from_file(filename, dimension);
    }
    this->_construct(intervals, filter, number_of_pixels, min, max);
  }

 private:
  // if min == max, then the program is requested to set up the values itself based on persistence intervals
  void _construct(const std::vector<std::pair<double, double> >& intervals,
                  const std::vector<std::vector<double> >& filter = create_Gaussian_filter(5, 1),
                  std::size_t number_of_pixels = 1000,
                  double min = -1,
                  double max = -1)
  {
#ifdef DEBUG_TRACES
    std::clog << "Entering construct procedure \n";
#endif

    if (min == max) {
      // in this case, we want the program to set up the min and max values by itself.
      min = std::numeric_limits<double>::max();
      max = -std::numeric_limits<double>::max();

      for (std::size_t i = 0; i != intervals.size(); ++i) {
        if (intervals[i].first < min) min = intervals[i].first;
        if (intervals[i].second > max) max = intervals[i].second;
      }
      // now we have the structure filled in, and moreover we know min and max values of the interval, so we know the
      // range.

      // add some more space:
      auto pad = std::fabs(max - min) / 100;
      min -= pad;
      max += pad;
    }

#ifdef DEBUG_TRACES
    std::clog << "min : " << min << std::endl;
    std::clog << "max : " << max << std::endl;
    std::clog << "number_of_pixels : " << number_of_pixels << std::endl;
#endif

    Base::min_ = min;
    Base::max_ = max;

    // initialization of the structure heat_map
    std::vector<std::vector<double> > heat_map;
    for (std::size_t i = 0; i != number_of_pixels; ++i) {
      std::vector<double> v(number_of_pixels, 0);
      heat_map.push_back(v);
    }
    Base::heat_map_.swap(heat_map);

#ifdef DEBUG_TRACES
    std::clog << "Done creating of the heat map, now we will fill in the structure \n";
#endif

    for (std::size_t pt_nr = 0; pt_nr != intervals.size(); ++pt_nr) {
      // compute the value of intervals_[pt_nr] in the grid:
      int x_grid =
          static_cast<int>((intervals[pt_nr].first - Base::min_) / (Base::max_ - Base::min_) * number_of_pixels);
      int y_grid =
          static_cast<int>((intervals[pt_nr].second - Base::min_) / (Base::max_ - Base::min_) * number_of_pixels);

#ifdef DEBUG_TRACES
      std::clog << "point : " << intervals[pt_nr].first << " , " << intervals[pt_nr].second << std::endl;
      std::clog << "x_grid : " << x_grid << std::endl;
      std::clog << "y_grid : " << y_grid << std::endl;
#endif

      // x_grid and y_grid gives a center of the kernel. We want to have its lower left corner. To get this, we need to
      // shift x_grid and y_grid by a grid diameter.
      x_grid -= filter.size() / 2;
      y_grid -= filter.size() / 2;
      // note that the numbers x_grid and y_grid may be negative.

#ifdef DEBUG_TRACES
      std::clog << "After shift : \n";
      std::clog << "x_grid : " << x_grid << std::endl;
      std::clog << "y_grid : " << y_grid << std::endl;
      std::clog << "filter.size() : " << filter.size() << std::endl;
#endif

      for (std::size_t i = 0; i != filter.size(); ++i) {
        for (std::size_t j = 0; j != filter.size(); ++j) {
          // if the point (x_grid+i,y_grid+j) is the correct point in the grid.
          if (((x_grid + i) >= 0) && (x_grid + i < Base::heat_map_.size()) && ((y_grid + j) >= 0) &&
              (y_grid + j < Base::heat_map_.size())) {
#ifdef DEBUG_TRACES
            std::clog << y_grid + j << " " << x_grid + i << std::endl;
#endif
            Base::heat_map_[y_grid + j][x_grid + i] += filter[i][j];
            Base::heat_map_[x_grid + i][y_grid + j] += -filter[i][j];
          }
        }
      }
    }
  }  // construct
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PSSK_H_
