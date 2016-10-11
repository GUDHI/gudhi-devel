
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

#ifndef Abs_Topological_data_with_averages_H_
#define Abs_Topological_data_with_averages_H_
#include <gudhi/abstract_classes/Abs_Topological_data.h>


namespace Gudhi 
{
namespace Gudhi_stat 
{

/**
* This is an abstract container to store topological information. Most typically, this information will be some representation of persistent homology.
**/

class Abs_Topological_data_with_averages : public virtual Abs_Topological_data
{
public:
    Abs_Topological_data_with_averages(){};
    
    //options: 
    //1) Function returns void, take std::vector< Abs_Topological_data_with_averages* > and change the object underneeth to the average of this vector.
    //2) Function returns  Abs_Topological_data*, and compute the average of all the objects in the vector, plus (*this). 
    //virtual Abs_Topological_data* compute_average( std::vector< Abs_Topological_data_with_averages* > to_average ) = 0;
    //At the moment I will try to implement option (1). 
    virtual void compute_average( const std::vector< Abs_Topological_data_with_averages* >& to_average ) = 0;
    virtual ~Abs_Topological_data_with_averages(){}
protected:
};


}//namespace Gudhi_stat
}//namespace Gudhi 


 


#endif
