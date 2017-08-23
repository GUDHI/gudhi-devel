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

#ifndef PERSISTENCE_VECTORS_INTERFACE_H_
#define PERSISTENCE_VECTORS_INTERFACE_H_

// gudhi include
#include <gudhi/persistence_vectors.h>
s
namespace Gudhi {
namespace Persistence_representations {



template <typename F>
class Vector_distances_in_diagram_interface : Vector_distances_in_diagram<Euclidean_distance> {
 public:
  Vector_distances_in_diagram_interface():Vector_distances_in_diagram(){}

  Vector_distances_in_diagram_interface(const std::vector<std::pair<double, double> >& intervals, size_t where_to_cut):
  Vector_distances_in_diagram(intervals,where_to_cut){}

  Vector_distances_in_diagram_interface(const char* filename, size_t where_to_cut,
                              unsigned dimension = std::numeric_limits<unsigned>::max()):
  Vector_distances_in_diagram_interface(filename,where_to_cut,dimension){}                             
                              
  inline double vector_in_position_interface(size_t position) const 
  {
	  return this->vector_in_position(position):
  }
    
  inline size_t size_interface() const 
  {
	  return this->size();
  }

  void write_to_file_interface(const char* filename) const
  {
	  this->write_to_file( filename );
  }
  
  void print_to_file_interface(const char* filename) const
  {
	  this->print_to_file(filename);
  }

  void load_from_file_interface(const char* filename)
  {
	  this->load_from_file(filename);
  }

  double project_to_R_interface(int number_of_function) const
  {
	  return this->project_to_R(number_of_function);
  }
  
  size_t number_of_projections_to_R_interface() const
  {
	  return this->number_of_projections_to_R();
  }

  std::vector<double> vectorize_interface(int number_of_function) const
  {
	  return this->vectorize(number_of_function);
  }
  
  size_t number_of_vectorize_functions_interface() const 
  {
	  return this->number_of_vectorize_functions();
  }

  void compute_average_interface(const std::vector<Vector_distances_in_diagram*>& to_average)
  {
	  this->compute_average(to_average);
  }

  double distance_interface(const Vector_distances_in_diagram& second, double power = 1) const
  {
	  return this->distance(second,power);
  }

  double compute_scalar_product_interface(const Vector_distances_in_diagram& second) const
  {
	  return this->compute_scalar_product(second);
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

#endif  // PERSISTENCE_VECTORS_INTERFACE_H_
