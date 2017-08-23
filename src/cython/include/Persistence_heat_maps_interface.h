/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016  INRIA (France)
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

#ifndef PERSISTENCE_HEAT_MAPS_INTERFACE_H_
#define PERSISTENCE_HEAT_MAPS_INTERFACE_H_

// gudhi include
#include <gudhi/read_persistence_from_file.h>
#include <gudhi/common_persistence_representations.h>

// standard include
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>
#include <string>
#include <functional>

namespace Gudhi {
namespace Persistence_representations {

class Persistence_heat_maps_interface : public Persistence_heat_maps {
 public:

  Persistence_heat_maps_interface():Persistence_heat_maps(){}

 
  Persistence_heat_maps_interface(const std::vector<std::pair<double, double> >& interval,
                        std::vector<std::vector<double> > filter = create_Gaussian_filter(5, 1),
                        bool erase_below_diagonal = false, size_t number_of_pixels = 1000,
                        double min_ = std::numeric_limits<double>::max(),
                        double max_ = std::numeric_limits<double>::max()):
   Persistence_heat_maps(interval,filter,erase_below_diagonal,number_of_pixels,min_max_){}                     

  
  Persistence_heat_maps_interface(const char* filename, std::vector<std::vector<double> > filter = create_Gaussian_filter(5, 1),
                        bool erase_below_diagonal = false, size_t number_of_pixels = 1000,
                        double min_ = std::numeric_limits<double>::max(),
                        double max_ = std::numeric_limits<double>::max(),
                        unsigned dimension = std::numeric_limits<unsigned>::max()):
  Persistence_heat_maps(filename,filter,erase_below_diagonal,number_of_pixels,min_,max_,dimension){}

  void compute_mean_interface(const std::vector<Persistence_heat_maps*>& maps)
  {
	  this->compute_mean(maps);
  }

  void compute_median_interface(const std::vector<Persistence_heat_maps*>& maps)
  {
	  this->compute_median(maps);
  }
  
  void compute_percentage_of_active_interface(const std::vector<Persistence_heat_maps*>& maps, size_t cutoff = 1)
  {
	  this->compute_percentage_of_active(maps,cutoff);
  }
  
  void print_to_file_interface(const char* filename) const
  {
	  this->print_to_file(filename);
  }  

  void load_from_file_interface(const char* filename)
  {
	  this->load_from_file( filename );
  }  

  inline bool check_if_the_same_interface(const Persistence_heat_maps& second) const 
  {
	  return this->check_if_the_same( second );
  }  

  inline double get_min_interface() const
  {
	  return this->get_min();
  }  

  inline double get_max_interface() const
  {
	  return this->get_max();
  }  

  std::vector<double> vectorize_interface(int number_of_function) const
  {
	  return this->vectorize(number_of_function);
  }  
  
  size_t number_of_vectorize_functions_interface() const
  {
	  return this->number_of_vectorize_functions();
  }  
  
  double project_to_R_interface(int number_of_function) const
  {
	  return this->project_to_R( number_of_function );
  }  
  
  size_t number_of_projections_to_R_interface() const
  {
	  return this->number_of_projections_to_R();
  }  

  double distance_interface(const Persistence_heat_maps& second_, double power = 1) const
  {
	  return this->distance( second, power );
  }  

  void compute_average_interface(const std::vector<Persistence_heat_maps*>& to_average)
  {
	  this->compute_average( to_average );
  }  
  
  double compute_scalar_product_interface(const Persistence_heat_maps& second_) const
  {
	  return this->compute_scalar_product( second_ );
  }  
  
  std::pair<double, double> get_x_range_interface() const
  {
	  return this->get_x_range();
  }  

  std::pair<double, double> get_y_range_interface() const
  {
	  return this->get_y_range();
  }  

};


}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_HEAT_MAPS_INTERFACE_H_
