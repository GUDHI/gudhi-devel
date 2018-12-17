/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
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

#ifndef CONCEPT_VECTORIZED_TOPOLOGICAL_DATA_H_
#define CONCEPT_VECTORIZED_TOPOLOGICAL_DATA_H_

namespace Gudhi {

namespace Persistence_representations {

/** \brief The concept Vectorized_topological_data describes the requirements
  * for a type to implement a container that allows vectorization.
  */
class Vectorized_topological_data {
 public:
  /**
   * There are various ways data can be vectorized. This function give us the number of functions for vectorization
   *provided by a given class.
  **/
  size_t number_of_vectorize_functions();
  /**
   * This is a function to vectorize given container. The parameter of a function have to be between 0 and the value
   *returned by number_of_vectorize_functions().
  **/
  std::vector<double> vectorize(int number_of_function);
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // CONCEPT_VECTORIZED_TOPOLOGICAL_DATA_H_
