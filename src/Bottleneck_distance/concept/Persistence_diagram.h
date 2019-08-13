/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author:       François Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_
#define CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_

namespace Gudhi {

namespace persistence_diagram {

/** \brief Concept of point in a persistence diagram. std::get<0>(point) must return the birth of the corresponding component and std::get<1>(point) its death.
 * Both should be convertible to `double`.
 * A valid implementation of this concept is std::pair<double,double>.
 * Death should be larger than birth, death can be std::numeric_limits<double>::infinity() for components which stay alive.
 *
 * \ingroup bottleneck_distance
 */
struct DiagramPoint{};

/** \brief Concept of persistence diagram. It is a range of `DiagramPoint`.
 * std::begin(diagram) and std::end(diagram) must return corresponding iterators.
 *
 * \ingroup bottleneck_distance
 */
struct PersistenceDiagram{};

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // CONCEPT_BOTTLENECK_DISTANCE_PERSISTENCE_DIAGRAM_H_
