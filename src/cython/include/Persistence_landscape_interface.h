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

#ifndef PERSISTENCE_LANDSCAPE_INTERFACE_H_
#define PERSISTENCE_LANDSCAPE_INTERFACE_H_

#include <gudhi/Persistence_landscape.h>

namespace Gudhi {
namespace Persistence_representations {


class Persistence_landscape_interface : public Persistence_landscape 
{
 public:
  Persistence_landscape_interface():Persistence_landscape(){}

  Persistence_landscape_interface(const std::vector<std::pair<double, double> >& p, size_t number_of_levels = std::numeric_limits<size_t>::max() ):Persistence_landscape(p,number_of_levels){}

  Persistence_landscape_interface(const char* filename, size_t dimension = std::numeric_limits<unsigned>::max() , size_t number_of_levels = std::numeric_limits<size_t>::max() ):Persistence_landscape(filename,dimension,number_of_levels){}
  
  //****************
  static Persistence_landscape_interface* construct_from_file( const char* filename, size_t dimension = std::numeric_limits<unsigned>::max() , size_t number_of_levels = std::numeric_limits<size_t>::max() )
  {
    Persistence_landscape_interface* result = new Persistence_landscape_interface(filename,dimension,number_of_levels);
	  return result;  
  }
  static Persistence_landscape_interface* construct_from_vector_of_pairs( const std::vector<std::pair<double, double> >& p, size_t number_of_levels = std::numeric_limits<size_t>::max() )
  {
    Persistence_landscape_interface* result = new Persistence_landscape_interface(p,number_of_levels);
	  return result;  
  }
    
  //****************

  
  Persistence_landscape_interface* new_abs_interface()
  {
	   return (Persistence_landscape_interface*)this->new_abs();
  }
  
  void new_compute_average(const std::vector<Persistence_landscape_interface*>& to_average) 
  {
	  std::vector<Persistence_landscape*> to_average_new;
	  to_average_new.reserve( to_average.size() );
	  for ( size_t i = 0 ; i != to_average.size() ; ++i )
	  {
		  to_average_new.push_back( (Persistence_landscape*)to_average[i] );
	  }
	  this->compute_average(to_average_new);
  }

/*
  void load_landscape_from_file_interface(const char* filename)
  {
	  this->load_landscape_from_file(filename);
  }


  void print_to_file_interface(const char* filename) const
  {
	  this->print_to_file(filename);
  }


  double compute_integral_of_landscape_interface() const
  {
	  return this->compute_integral_of_landscape();
  }


  double compute_integral_of_a_level_of_a_landscape_interface(size_t level) const
  {
	  return this->compute_integral_of_a_level_of_a_landscape(level);
  }


  double compute_integral_of_landscape_interface(double p) const
  {
	  return this->compute_integral_of_landscape(p);
  }


  double compute_value_at_a_given_point_interface(unsigned level, double x) const
  {
	  return this->compute_value_at_a_given_point(level,x);
  }


  double compute_maximum_interface() const
  {
	  return this->compute_maximum();
  }
  
  
  double compute_minimum_interface() const 
  {
	  return this->compute_minimum();
  }
  
  
  
  double compute_norm_of_landscape_interface(double i) 
  {
	  return this->compute_norm_of_landscape(i);
  }
  
  
  Persistence_landscape abs_interface()
  {
	  return this->abs();
  }

  size_t size_interface() const
  {
	  return this->size();
  }

  double find_max_interface(unsigned lambda) const
  {
	  return this->find_max();
  }

  friend double compute_inner_product_interface(const Persistence_landscape& l1, const Persistence_landscape& l2)
  {
	  return this->compute_inner_product(l1,l2);
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
	  return this->vectorize( number_of_function );
  }
  
  
  size_t number_of_vectorize_function_interface() const 
  {
	  return this->number_of_vectorize_function();
  }
  
  
  void compute_average_interface(const std::vector<Persistence_landscape*>& to_average) 
  {
	  return this->compute_average(to_average);
  }
  
  
  double distance_interface(const Persistence_landscape& second, double power = 1) 
  {
	  return this->distance( second, power );
  }
 
 
  double compute_scalar_product_interface(const Persistence_landscape& second) const 
  {
	  return this->compute_scalar_product( second );
  }


  std::pair<double, double> get_y_range_interface(size_t level = 0) const 
  {
	  return this->get_y_range( level );
  }
  */
};

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // PERSISTENCE_LANDSCAPE_INTERFACE_H_
