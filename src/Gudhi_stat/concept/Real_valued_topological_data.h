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

/** \brief The concept Real_valued_topological_data describes the requirements 
  * for a type to implement a container that allows computations of its projections to R. 
  */
class Real_valued_topological_data
{
public:
	/**
	 * Typically there are various ways data can be projected to R. This function give us the number of functions for vectorization provided by a given class. 
	**/ 
	 size_t number_of_projections_to_R();
	   /**
     * This is a function to compute the projection from this container to reals. The parameter of a function have to be between 0 and the value returned by number_of_projections_to_R().
    **/ 
     double project_to_R( size_t number_of_projection );
};

