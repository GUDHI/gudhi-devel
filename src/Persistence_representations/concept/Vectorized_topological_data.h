/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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
