/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_
#define DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_

namespace Gudhi {

namespace cech_complex {

/**  \defgroup cech_complex Čech complex
 * 
 * \author    Vincent Rouvreau, Hind montassif, Marc Glisse
 * 
 * @{
 * 
 * \section cechdefinition Čech complex definition
 * 
 * Čech complex
 * <a target="_blank" href="https://en.wikipedia.org/wiki/%C4%8Cech_complex">(Wikipedia)</a> is a
 * <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a> constructed
 *  from a proximity graph. The set of all simplices is filtered by the radius of their minimal enclosing ball.
 *
 * The input shall be a range of points where a point is defined as <a target="_blank" href="https://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Point__d.html">CGAL kernel Point_d.</a>
 * 
 * \remark For people only interested in the topology of the \ref cech_complex (for instance persistence),
 * \ref alpha_complex is equivalent to the \ref cech_complex and much smaller if you do not bound the radii.
 * \ref cech_complex can still make sense in higher dimension precisely because you can bound the radii.
 *
 * \subsection cechalgorithm Algorithm
 *
 * Cech_complex first builds a proximity graph from a point cloud.
 * The filtration value of each edge of the `Gudhi::Proximity_graph` is computed using CGAL kernel functions.
 * 
 * All edges that have a filtration value strictly greater than a user given maximal radius value, \f$max\_radius\f$,
 * are not inserted into the complex.
 * 
 * Vertex name correspond to the index of the point in the given range (aka. the point cloud).
 * 
 * \image html "cech_one_skeleton.png" "Čech complex proximity graph representation"
 * 
 * When creating a simplicial complex from this proximity graph, Cech_complex inserts the proximity graph into the
 * simplicial complex data structure, and then expands the simplicial complex when required.
 *
 * On this example, as edges \f$(x,y)\f$, \f$(y,z)\f$ and \f$(z,y)\f$ are in the complex, the minimal ball radius
 * containing the points \f$(x,y,z)\f$ is computed.
 *
 * \f$(x,y,z)\f$ is inserted to the simplicial complex with the filtration value set with
 * \f$mini\_ball\_radius(x,y,z))\f$ iff \f$mini\_ball\_radius(x,y,z)) \leq max\_radius\f$.
 *
 * And so on for higher dimensions.
 *
 * \image html "cech_complex_representation.png" "Čech complex expansion"
 *
 * This radius computation is the reason why the Cech_complex is taking much more time to be computed than the
 * \ref rips_complex but it offers more topological guarantees.
 *
 * If you already have a simplicial complex, it is possible to assign to each simplex a filtration value corresponding
 * to the squared radius of its minimal enclosing ball using `assign_MEB_filtration()`. This can provide an alternate
 * way of computing a Čech filtration, or it can be used on a Delaunay triangulation to compute a Delaunay-Čech
 * filtration.
 *
 * \subsection cechpointscloudexample Example from a point cloud
 * 
 * This example builds the proximity graph from the given points, and maximal radius values.
 * Then it creates a `Simplex_tree` with it.
 * 
 * Then, it is asked to display information about the simplicial complex.
 * 
 * \include cech_complex_example_from_points.cpp
 * 
 * When launching (maximal enclosing ball radius is 1., is expanded until dimension 2):
 * 
 * \code $> ./Cech_complex_example_from_points
 * \endcode
 *
 * the program output is:
 * 
 * \include cech_complex_example_from_points_for_doc.txt
 * 
 */
/** @} */  // end defgroup cech_complex

}  // namespace cech_complex

}  // namespace Gudhi

#endif  // DOC_CECH_COMPLEX_INTRO_CECH_COMPLEX_H_
