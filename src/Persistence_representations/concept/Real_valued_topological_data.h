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

#ifndef CONCEPT_REAL_VALUED_TOPOLOGICAL_DATA_H_
#define CONCEPT_REAL_VALUED_TOPOLOGICAL_DATA_H_

namespace Gudhi {

namespace Persistence_representations {

/** \brief The concept Real_valued_topological_data describes the requirements
  * for a type to implement a container that allows computations of its projections to R.
  */
class Real_valued_topological_data {
 public:
  /**
* Typically there are various ways data can be projected to R. This function gives us the number of functions for
* vectorization provided by a given class.
  **/
  size_t number_of_projections_to_R();
  /**
* This is a function to compute the projection from this container to reals. The parameter of a function have to
* be between 0 and the value returned by number_of_projections_to_R().
**/
  double project_to_R(size_t number_of_projection);
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // CONCEPT_REAL_VALUED_TOPOLOGICAL_DATA_H_
