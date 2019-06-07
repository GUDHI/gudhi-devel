/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author:       François Godi
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modifications:
 *      - 2019/06 Vincent Rouvreau : Fix #11 - Distance computation shall be better documented.
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
 * diagonal in order to ignore cardinality mismatchs) such that any couple of matched points are at distance at most b,
 * where the distance between points is the sup norm in \f$\mathbb{R}^2\f$ (not the Euclidean distance).
 *
 * \image html perturb_pd.png On this picture, the red edges represent the matching. The bottleneck distance is the length of the longest edge.
 *
 * This implementation is based on ideas from "Geometry Helps in Bottleneck Matching and Related Problems"
 * \cite DBLP:journals/algorithmica/EfratIK01. Another relevant publication, although it was not used is
 * "Geometry Helps to Compare Persistence Diagrams" \cite Kerber:2017:GHC:3047249.3064175.
 * 
 * \section bottleneckdistanceprecision Distance computation
 *
 * The following example explains how the distance is computed:
 *
 * \code{.cpp}
#include <gudhi/Bottleneck.h>

#include <iostream>
#include <vector>
#include <utility>  // for pair

int main() {
  std::vector< std::pair<double, double> > diag1, diag2;
  diag1.emplace_back(0., 0.);
  diag2.emplace_back(0., 13.);

  double b = Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2);
  std::cout << "Bottleneck distance = " << b << std::endl;
}
 * \endcode
 *
 * \code Bottleneck distance = 6.5
 * \endcode
 *
 * \image html bottleneck_distance_example.png The point (0, 13) is at distance 6.5 from the diagonal and more specifically from the point (6.5, 6.5)
 *
 * \section bottleneckbasicexample Basic example
 *
 * This other example computes the bottleneck distance from 2 persistence diagrams:
 * \include Bottleneck_distance/bottleneck_basic_example.cpp
 *
 * \code
    Bottleneck distance = 0.75
    Approx bottleneck distance = 0.808176
 * \endcode

 */
/** @} */  // end defgroup bottleneck_distance

}  // namespace persistence_diagram

}  // namespace Gudhi

#endif  // DOC_BOTTLENECK_DISTANCE_INTRO_BOTTLENECK_DISTANCE_H_
