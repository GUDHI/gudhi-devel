/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA (France)
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

#ifndef Abs_Vectorized_topological_data_H_
#define Abs_Vectorized_topological_data_H_

#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include "Abs_Topological_data.h"
#include <vector>

namespace Gudhi 
{
namespace Gudhi_stat 
{
using namespace std;

/**
* This is specialization of a topological data class that allows to create vectors based on topological data.
**/

class Abs_Vectorized_topological_data  :  public virtual Abs_Topological_data 	
{
public:
    Abs_Vectorized_topological_data():where_to_cut(10){}
    Abs_Vectorized_topological_data( size_t where_to_cut ):where_to_cut(where_to_cut){}
    size_t number_of_vectorize_functions(){return number_of_functions_for_vectorization;}
    virtual std::vector<double> vectorize( int number_of_function ) = 0;
    virtual ~Abs_Vectorized_topological_data(){}
protected:
    size_t number_of_functions_for_vectorization;
    size_t where_to_cut;
};

}//namespace Gudhi_stat
}//namespace Gudhi 

#endif
