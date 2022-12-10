/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_SUBSAMPLING_INTRO_SUBSAMPLING_H_
#define DOC_SUBSAMPLING_INTRO_SUBSAMPLING_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace subsampling {

/**  \defgroup subsampling Subsampling
 * 
 * \author Cl&eacute;ment Jamin, Siargey Kachanovich, Marc Glisse
 * 
 * @{
 * 
 * \section subsamplingintroduction Introduction
 * 
 * This Gudhi component offers methods to subsample a set of points.
 * 
 * \section sparsifyexamples Example: sparsify_point_set
 * 
 * This example outputs a subset of the input points so that the 
 * squared distance between any two points
 * is greater than or equal to 0.4.
 * 
 * \include example_sparsify_point_set.cpp
 * 
 * \section farthestpointexamples Example: choose_n_farthest_points
 *
 * This example outputs a subset of 100 points obtained by Gonz&aacute;lez algorithm,
 * starting with a random point, and greedily adding the point farthest from those already picked.
 *
 * \include example_choose_n_farthest_points.cpp
 * 
 * There is an alternative version `choose_n_farthest_points_metric()` with the same syntax,
 * which can be faster in many cases.
 *
 * \section randompointexamples Example: pick_n_random_points
 *
 * This example outputs a subset of 100 points picked randomly.
 *
 * \include example_pick_n_random_points.cpp
 */
/** @} */  // end defgroup subsampling

}  // namespace subsampling

}  // namespace Gudhi

#endif  // DOC_SUBSAMPLING_INTRO_SUBSAMPLING_H_
