/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_
#define DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
// needs namespace for Doxygen to link on classes
namespace alpha_complex {

/**  \defgroup alpha_complex Alpha complex
 *
 * \author    Vincent Rouvreau
 *
 * @{
 *
<div class="toc">
Table of Contents
<ul>
<li class="level1"><a href="#definition">Definition</a></li>
<li class="level1"><a href="#pointsexample">Example from points</a></li>
<li class="level1"><a href="#createcomplexalgorithm">Create complex algorithm</a></li>
<li class="level1"><a href="#offexample">Example from OFF file</a></li>
<li class="level1"><a href="#weighted3dexample">3d specific version</a></li>
</ul>
</div>

 * \section definition Definition
 *
 * Alpha_complex is a <a target="_blank" href="https://en.wikipedia.org/wiki/Simplicial_complex">simplicial complex</a>
 * constructed from the finite cells of a Delaunay Triangulation.
 *
 * The filtration value of each simplex is computed as the square of the circumradius of the simplex if the
 * circumsphere is empty (the simplex is then said to be Gabriel), and as the minimum of the filtration
 * values of the codimension 1 cofaces that make it not Gabriel otherwise.
 *
 * All simplices that have a filtration value \f$ > \alpha^2 \f$ are removed from the Delaunay complex
 * when creating the simplicial complex if it is specified.
 *
 * \image html "alpha_complex_representation.png" "Alpha-complex representation"
 *
 * Alpha_complex is constructing a <a target="_blank"
 * href="http://doc.cgal.org/latest/Triangulation/index.html#Chapter_Triangulations">Delaunay Triangulation</a>
 * \cite cgal:hdj-t-19b from <a target="_blank" href="http://www.cgal.org/">CGAL</a> (the Computational Geometry
 * Algorithms Library \cite cgal:eb-19b) and is able to create a `SimplicialComplexForAlpha`.
 *
 * The complex is a template class requiring an Epick_d <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/index.html#Chapter_dD_Geometry_Kernel">dD Geometry Kernel</a>
 * \cite cgal:s-gkd-19b from CGAL as template parameter.
 *
 * \remark
 * - When the simplicial complex is constructed with an infinite value of \f$ \alpha^2 \f$, the complex is a Delaunay
 * complex with special filtration values. The Delaunay complex without filtration values is also available by passing
 * `default_filtration_value=true` to `Alpha_complex::create_complex`.
 * - For people only interested in the topology of the \ref alpha_complex (for instance persistence),
 * \ref alpha_complex is equivalent to the \ref cech_complex and much smaller if you do not bound the radii.
 * \ref cech_complex can still make sense in higher dimension precisely because you can bound the radii.
 * - Using the default `CGAL::Epeck_d` makes the construction safe. If you pass `exact=true` to
 * `Alpha_complex::create_complex`, the filtration values are the exact ones converted to the filtration value type of
 * the simplicial complex. This can be very slow. If you pass `exact=false` (the default), the filtration values are
 * only guaranteed to have a small multiplicative error compared to the exact value, see <code>
 * <a class="el" target="_blank" href="https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html">
 * CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double</a></code> for details. A drawback, when computing
 * persistence, is that an empty exact interval [10^12,10^12] may become a non-empty approximate interval
 * [10^12,10^12+10^6]. Using `CGAL::Epick_d` makes the computations slightly faster, and the combinatorics are still
 * exact, but the computation of filtration values can exceptionally be arbitrarily bad. In all cases, we still
 * guarantee that the output is a valid filtration (faces have a filtration value no larger than their cofaces).
 * - For performances reasons, it is advised to use `Alpha_complex` with \ref cgal &ge; 5.0.0.
 *
 * \section pointsexample Example from points
 *
 * This example builds the Delaunay triangulation from the given points in a 2D static kernel, and creates a
 * `Simplex_tree` with it.
 *
 * Then, it is asked to display information about the simplicial complex.
 *
 * \include Alpha_complex/Alpha_complex_from_points.cpp
 *
 * When launching:
 *
 * \code $> ./Alpha_complex_example_from_points
 * \endcode
 *
 * the program output is:
 *
 * \include Alpha_complex/alphaoffreader_for_doc_60.txt
 *
 * \section createcomplexalgorithm Create complex algorithm
 *
 * \subsection datastructure Data structure
 *
 * In order to create the simplicial complex, first, it is built from the cells of the Delaunay Triangulation.
 * The filtration values are set to NaN, which stands for unknown value.
 *
 * In example, :
 * \image html "alpha_complex_doc.png" "Simplicial complex structure construction example"
 *
 * \subsection filtrationcomputation Filtration value computation algorithm
 * <br>
 * \f$
 * \textbf{for } \text{i : dimension } \rightarrow 0 \textbf{ do}\\
 * \quad \textbf{for all } \sigma \text{ of dimension i}\\
 * \quad\quad \textbf{if } \text{filtration(} \sigma ) \text{ is NaN} \textbf{ then}\\
 * \quad\quad\quad \text{filtration(} \sigma ) = \alpha^2( \sigma )\\
 * \quad\quad \textbf{end if}\\
 * \quad\quad \textbf{for all } \tau \text{ face of } \sigma \textbf{ do}\quad\quad
 * \textit{// propagate alpha filtration value}\\
 * \quad\quad\quad \textbf{if } \text{filtration(} \tau ) \text{ is not NaN} \textbf{ then}\\
 * \quad\quad\quad\quad \text{filtration(} \tau \text{) = min( filtration(} \tau \text{), filtration(} \sigma
 * \text{) )}\\
 * \quad\quad\quad \textbf{else}\\
 * \quad\quad\quad\quad \textbf{if } \tau \text{ is not Gabriel for } \sigma \textbf{ then}\\
 * \quad\quad\quad\quad\quad \text{filtration(} \tau \text{) = filtration(} \sigma \text{)}\\
 * \quad\quad\quad\quad \textbf{end if}\\
 * \quad\quad\quad \textbf{end if}\\
 * \quad\quad \textbf{end for}\\
 * \quad \textbf{end for}\\
 * \textbf{end for}\\
 * \text{make_filtration_non_decreasing()}\\
 * \text{prune_above_filtration()}\\
 * \f$
 *
 * \subsubsection dimension2 Dimension 2
 *
 * From the example above, it means the algorithm looks into each triangle ([0,1,2], [0,2,4], [1,2,3], ...),
 * computes the filtration value of the triangle, and then propagates the filtration value as described
 * here :
 * \image html "alpha_complex_doc_420.png" "Filtration value propagation example"
 *
 * \subsubsection dimension1 Dimension 1
 *
 * Then, the algorithm looks into each edge ([0,1], [0,2], [1,2], ...),
 * computes the filtration value of the edge (in this case, propagation will have no effect).
 *
 * \subsubsection dimension0 Dimension 0
 *
 * Finally, the algorithm looks into each vertex ([0], [1], [2], [3], [4], [5] and [6]) and
 * sets the filtration value (0 in case of a vertex - propagation will have no effect).
 *
 * \subsubsection nondecreasing Non decreasing filtration values
 *
 * As the squared radii computed by CGAL are an approximation, it might happen that these \f$ \alpha^2 \f$ values do
 * not quite define a proper filtration (i.e. non-decreasing with respect to inclusion).
 * We fix that up by calling `SimplicialComplexForAlpha::make_filtration_non_decreasing()`.
 *
 * \subsubsection pruneabove Prune above given filtration value
 *
 * The simplex tree is pruned from the given maximum \f$ \alpha^2 \f$ value (cf.
 * `SimplicialComplexForAlpha::prune_above_filtration()`).
 * In the following example, the value is given by the user as argument of the program.
 *
 * \section weightedversion Weighted specific version
 * <b>Requires:</b> \ref eigen &ge; 3.1.0 and \ref cgal &ge; 5.1.0.
 * For performances reasons, it is advised to use \ref eigen &ge; 3.5.8 and \ref cgal &ge; 5.2.0
 * A weighted version for Alpha complex is available (cf. Alpha_complex). 
 *
 * This example builds the CGAL weighted alpha shapes from a small molecule, and initializes the alpha complex with
 * it. This example is taken from <a href="https://doc.cgal.org/latest/Alpha_shapes_3/index.html#title13">CGAL 3d
 * weighted alpha shapes</a>.
 *
 * Then, it is asked to display information about the alpha complex.
 *
 * \include Alpha_complex/Weighted_alpha_complex_from_points.cpp
 *
 * When launching:
 *
 * \code $> ./Weighted_alpha_complex_example_from_points
 * \endcode
 *
 * the program output is:
 *
 * \include Alpha_complex/weightedalpha3dfrompoints_for_doc.txt
 *
 *
 * \section offexample Example from OFF file
 *
 * This example builds the Delaunay triangulation in a dynamic kernel, and initializes the alpha complex with it.
 *
 *
 * Then, it is asked to display information about the alpha complex.
 *
 * \include Alpha_complex/Alpha_complex_from_off.cpp
 *
 * When launching:
 *
 * \code $> ./Alpha_complex_example_from_off ../../data/points/alphacomplexdoc.off 32.0
 * \endcode
 *
 * the program output is:
 *
 * \include Alpha_complex/alphaoffreader_for_doc_32.txt
 *
 *
 * \section weighted3dexample 3d specific version
 *
 * A specific module for Alpha complex is available in 3d (cf. Alpha_complex_3d) and allows to construct standard,
 * weighted, periodic or weighted and periodic versions of alpha complexes. Alpha values computation can be
 * Gudhi::alpha_complex::complexity::FAST, Gudhi::alpha_complex::complexity::SAFE (default value) or
 * Gudhi::alpha_complex::complexity::EXACT.
 *
 * This example builds the CGAL 3d weighted alpha shapes from a small molecule, and initializes the alpha complex with
 * it. This example is taken from <a href="https://doc.cgal.org/latest/Alpha_shapes_3/index.html#title13">CGAL 3d
 * weighted alpha shapes</a>.
 *
 * Then, it is asked to display information about the alpha complex.
 *
 * \include Alpha_complex/Weighted_alpha_complex_3d_from_points.cpp
 *
 * When launching:
 *
 * \code $> ./Alpha_complex_example_weighted_3d_from_points
 * \endcode
 *
 * the program output is:
 *
 * \include Alpha_complex/weightedalpha3dfrompoints_for_doc.txt
 *
 */
/** @} */  // end defgroup alpha_complex

}  // namespace alpha_complex

namespace alphacomplex = alpha_complex;

}  // namespace Gudhi

#endif  // DOC_ALPHA_COMPLEX_INTRO_ALPHA_COMPLEX_H_
