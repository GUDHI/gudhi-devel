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

#ifndef DOC_BOTTLENECK_DISTANCE_INTRO_BOTTLENECK_DISTANCE_H_
#define DOC_BOTTLENECK_DISTANCE_INTRO_BOTTLENECK_DISTANCE_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace persistence_diagram {

/**  \defgroup bottleneck_distance Bottleneck distance
 * 
 * \author    Fran&ccedil;ois Godi
 * @{
 * 
 * \section bottleneckdefinition Definition
 * 
 * The bottleneck distance measures the similarity between two persistence diagrams. It is the shortest distance b for
 * which there exists a perfect matching between the points of the two diagrams (completed with all the points on the
 * diagonal in order to ignore cardinality mismatchs) such that any couple of matched points are at distance at most b.
 *
 * \image html perturb_pd.png On this picture, the red edges represent the matching. The bottleneck distance is the length of the longest edge.
 *
 * This implementation is based on ideas from "Geometry Helps in Bottleneck Matching and Related Problems"
 * \cite DBLP:journals/algorithmica/EfratIK01. Another relevant publication, although it was not used is
 * "Geometry Helps to Compare Persistence Diagrams" \cite Kerber:2017:GHC:3047249.3064175.
 */
/** @} */  // end defgroup bottleneck_distance

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // DOC_BOTTLENECK_DISTANCE_INTRO_BOTTLENECK_DISTANCE_H_
