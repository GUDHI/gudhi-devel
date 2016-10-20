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

#ifndef Abs_Topological_data_with_scalar_product_H_
#define Abs_Topological_data_with_scalar_product_H_

#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <vector>

namespace Gudhi 
{
namespace Gudhi_stat 
{

/**
* This is an abstract container to store topological information. Most typically, this information will be some representation of persistent homology.
**/

class Abs_Topological_data_with_scalar_product : public virtual Abs_Topological_data
{
public:
    Abs_Topological_data_with_scalar_product(){};
    virtual double compute_scalar_product( const Abs_Topological_data_with_scalar_product* second ) = 0;
    virtual ~Abs_Topological_data_with_scalar_product(){}
protected:
};

}//namespace Gudhi_stat
}//namespace Gudhi 

#endif
