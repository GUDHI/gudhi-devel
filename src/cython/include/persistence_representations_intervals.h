/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University
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

#ifndef INCLUDE_PERSISTENCE_REPRESENTATIONS_INTERVALS_
#define INCLUDE_PERSISTENCE_REPRESENTATIONS_INTERVALS_

#include <gudhi/Persistence_intervals.h>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {
	
//But if we want to have the same names of classes in C++ and cyton side we ned this interface, because othervise we will have a name conflict. And we want to have the same names on the 
//C++ and python side for various reasonc (clarity, documentantions etc.).
//If the C++ class we inherid from are template class, we are inherid from concretization, for instance Persistence_intervals<double>.
//Also in this class, we create an interface functions that will be used in the python side. That will allow to have the same name of the functions in the C++ and python side. 

namespace Persistence_representations {

class Persistence_intervals_interface : public Persistence_intervals 
{
 public: 
	Persistence_intervals_interface(const char* filename, unsigned dimension = std::numeric_limits<unsigned>::max())
	: Persistence_intervals(filename, dimension) {
		std::cout << "CONTSTRUCTOR \n";
	}

	Persistence_intervals_interface(const std::vector<std::pair<double, double> >& intervals)
	: Persistence_intervals(intervals) {
	}

	std::pair<double, double> get_x_range_interface() const 
	{
		std::pair<double, double> range = 	this->get_x_range();
		std::cout << std::endl << std::endl << range.first << " " << range.second << std::endl << std::endl <<std::endl;
		return range;
	}

	std::pair<double, double> get_y_range_interface() const 	
	{
		return this->get_y_range();
	}

	std::vector<double> length_of_dominant_intervals_interface(size_t where_to_cut = 100) const	
	{
		return this->length_of_dominant_intervals(where_to_cut);
	}

	std::vector<std::pair<double, double> > dominant_intervals_interface(size_t where_to_cut = 100) const
	{
		return this->dominant_intervals(where_to_cut);
	}	

	std::vector<size_t> histogram_of_lengths_interface(size_t number_of_bins = 10) const
	{
		return this->histogram_of_lengths(number_of_bins);
	}

	std::vector<size_t> cumulative_histogram_of_lengths_interface(size_t number_of_bins = 10) const
	{
		return this->cumulative_histogram_of_lengths(number_of_bins);
	}

	std::vector<double> characteristic_function_of_diagram_interface(double x_min, double x_max, size_t number_of_bins = 10) const
	{
		return this->characteristic_function_of_diagram(x_min,x_max,number_of_bins);
	}	

	std::vector<double> cumulative_characteristic_function_of_diagram_interface(double x_min, double x_max, size_t number_of_bins = 10) const
	{
		return this->cumulative_characteristic_function_of_diagram(x_min,x_max,number_of_bins);
	}
                                                                    
	std::vector<std::pair<double, size_t> > compute_persistent_betti_numbers_interface() const
	{
		return this->compute_persistent_betti_numbers();
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
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // INCLUDE_PERSISTENCE_REPRESENTATIONS_DIAGRAMS_

