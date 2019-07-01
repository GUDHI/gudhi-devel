
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

#ifndef CONCEPT_TOPOLOGICAL_DATA_WITH_AVERAGES_H_
#define CONCEPT_TOPOLOGICAL_DATA_WITH_AVERAGES_H_

namespace Gudhi {

namespace Persistence_representations {

/** \brief The concept Topological_data_with_averages describes the requirements
  * for a type to implement a container that allows computations of averages.
  * Note that the average object after being computed is stored in *this.
  */
class Topological_data_with_averages {
 public:
  void compute_average(const std::vector<Topological_data_with_averages*>& to_average);
};

}  // namespace Persistence_representations

}  // namespace Gudhi

#endif  // CONCEPT_TOPOLOGICAL_DATA_WITH_AVERAGES_H_
