/*    This file is part of the Gudhi hiLibrary. The Gudhi library
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

#ifndef Persistence_intervals_WITH_DISTANCES_H_
#define Persistence_intervals_WITH_DISTANCES_H_


#include <gudhi/persistence_representations/Persistence_intervals.h>
#include <gudhi/Bottleneck.h>

namespace Gudhi 
{
namespace Gudhi_stat 
{

class Persistence_intervals_with_distances : public Persistence_intervals
{
public:
	/**
	 * This is a constructor of a class Persistence_intervals_with_distances from a text file. Each line of the input file is supposed to contain two numbers of a type doube (or convertable to double)
	 * representing the birth and the death of the persistence interval. If the pairs are not sorted so that birth <= death, then the constructor will sort then that way.
	 * The second parameter of a constructor is a dimension of intervals to be read from a file. If your file contains only birt-death pairs, use the default value. 
	**/ 
    Persistence_intervals_with_distances( const char* filename , unsigned dimension = std::numeric_limits<unsigned>::max() ):Persistence_intervals( filename , dimension ){}
    
    /**
     * This is a constructor of a class Persistence_intervals_with_distances from a vector of pairs. Each pair is assumed to represent a persistence interval. We assume that the first elemnets of pairs 
     * are smaller or equal the second elements of pairs.
    **/ 
    Persistence_intervals_with_distances( const std::vector< std::pair< double,double > >& intervals ):Persistence_intervals( intervals ){}

    /**
     *Computations of distance from the current persistnce diagram to the persistence diagram given as a parameter of this function.
     *The last but one parameter, power, is here in case we would like to compute p=th Wasserstein distance. At the moment, this method only implement Bottleneck distance,
     * which is infinity Wasserstein distance. Therefore any power which is not the default std::numeric_limits< double >::max() will be ignored and an
     * exception will be thrown.
     * The last parameter, tolerance, it is an additiv error of the approimation, set by default to zero.
    **/
     double distance( const Persistence_intervals_with_distances& second , double power = std::numeric_limits< double >::max() , double tolerance = 0) const
    {
		if ( power >= std::numeric_limits< double >::max() )
		{
			return Gudhi::persistence_diagram::bottleneck_distance(this->intervals, second.intervals, tolerance);			
		}
		else
		{
			std::cerr << "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance." << std::endl;
			throw "At the moment Gudhi do not support Wasserstein distances. We only support Bottleneck distance.";
		}
	}
};


}//namespace gudhi stat
}//namespace gudhi

#endif
