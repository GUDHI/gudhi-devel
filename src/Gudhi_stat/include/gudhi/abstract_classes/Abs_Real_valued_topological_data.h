/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

#ifndef Abs_Real_valued_topological_data_H_
#define Abs_Real_valued_topological_data_H_

#include <gudhi/abstract_classes/Abs_Topological_data.h>

namespace Gudhi 
{
namespace Gudhi_stat 
{



/**
* This is specialization of a topological data class allows computing various real-valued characterizations of topological data.
**/

class Abs_Real_valued_topological_data : public virtual Abs_Topological_data
{
public:
	 Abs_Real_valued_topological_data():number_of_functions_for_projections_to_reals(0){} 
	 Abs_Real_valued_topological_data( size_t number_of_functions_ ):number_of_functions_for_projections_to_reals(number_of_functions_){} 
	 size_t number_of_projections_to_R(){return this->number_of_functions_for_projections_to_reals;};
     virtual double project_to_R( int number_of_function ) = 0;
     virtual ~Abs_Real_valued_topological_data(){}
protected:
	size_t number_of_functions_for_projections_to_reals;
};

}//namespace Gudhi_stat
}//namespace Gudhi 

#endif
