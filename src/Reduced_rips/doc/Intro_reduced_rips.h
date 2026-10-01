/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_
#define DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_

namespace Gudhi {

namespace reduced_rips {

/**  @defgroup reduced_rips Reduced Vietoris-Rips degree-1 persistent homology
 *
 * @author    Thomas Burnett, Musashi Koyama
 *
 * @{
 *
 * @section reducedripsdefinition Definition
 *
 * This module computes the degree-1 (i.e. @f$H_1@f$) Vietoris-Rips persistent homology, with coefficients in
 * @f$\mathbb{Z}/2\mathbb{Z}@f$, of a point cloud or a distance matrix, without ever building the full
 * Vietoris-Rips complex. It implements the <em>Reduced Vietoris-Rips filtration</em> of Koyama, M&eacute;moli,
 * Robins and Turner @cite koyama2026computationdegree1persistenthomology.
 *
 * The @f$H_1@f$ barcode of the Vietoris-Rips filtration is determined by a small part of the complex, and the
 * reduction builds only that part:
 * - <b>Edges (1-simplices).</b> The only edges that can create or fill an @f$H_1@f$ class are those of the
 *   <em>relative neighborhood graph</em> (RNG). The <em>lune</em> of an edge @f$(a,b)@f$ of length
 *   @f$r = d(a,b)@f$ is the set of points strictly closer than @f$r@f$ to <em>both</em> @f$a@f$ and @f$b@f$;
 *   the edge belongs to the RNG exactly when its lune is empty. An edge whose lune contains a point is
 *   "dominated" and never needed for @f$H_1@f$.
 * - <b>Triangles (2-simplices).</b> Only the triangles certified by the same lune construction are added to
 *   the boundary matrix and reduced, instead of every triangle of the complex.
 *
 * Because so few simplices are ever materialised, the reduction scales to far larger inputs than building the
 * full complex would, while returning the <em>exact</em> degree-1 barcode. This holds for any symmetric matrix
 * of non-negative dissimilarities.
 *
 *
 * @image html "reduced_rips_lune.png" "The lune of an edge (a,b), which contains two connected components"
 *
 * Each connected component of a non-empty lune corresponds to one 2-simplex boundary column for the reduction, so
 * the edge pictured above contributes two columns. An edge with an empty lune belongs to the RNG and can only
 * give birth to a cycle.
 *
 * @section reducedripsinput Input
 *
 * Two inputs are accepted:
 * - a <b>Euclidean point cloud</b> (@ref Gudhi::reduced_rips::Reduced_rips::from_points), for which the
 *   relative neighborhood graph is built from a Delaunay triangulation in dimension 2 and 3 and a direct
 *   @f$O(n^2)@f$ construction otherwise.
 * - an <b>arbitrary symmetric distance matrix</b> (@ref Gudhi::reduced_rips::Reduced_rips::from_distance_matrix),
 *   for which the same reduction is driven purely by the supplied distances.
 *
 * Both factories take an optional initial per-point neighbor budget. The default of @f$\sqrt{n}@f$ (the
 * reference paper's choice) grows automatically whenever a point needs more neighbors, so it only tunes the
 * starting allocation.
 *
 * The class is templated on its `Filtration_value` type: `from_points` requires a floating-point type
 * (`double` by default), while `from_distance_matrix` accepts any arithmetic type, including exact integer
 * distances, and then computes the barcode exactly in that type.
 *
 *
 * @section reducedripsscope Scope and limitations
 *
 * - Only homological dimension 1 is computed (no @f$H_0@f$ or higher degrees).
 * - Coefficients are fixed to the field @f$\mathbb{Z}/2\mathbb{Z}@f$.
 * - The relative-neighborhood-graph construction costs @f$O(n^2)@f$ time in ambient dimension @f$\geq 4@f$
 *   and for distance matrices (in dimension 2 and 3 it is Delaunay-based and much cheaper). The reduction
 *   itself is output-sensitive, driven by the number of edges processed before the last bar closes.
 *
 * @section reducedripsexamples Examples
 *
 * @subsection reducedripsminimal Minimalistic examples
 *
 * \li \gudhi_example_link{Reduced_rips,example_reduced_rips_from_points.cpp} - Degree-1 persistence of a
 * hand-typed point cloud (the four corners of a unit square).
 * <details>
 *   @dontinclude example_reduced_rips_from_points.cpp
 *   @skip #include
 *   @until return 0;
 *   @skipline }
 * </details>
 * \li \gudhi_example_link{Reduced_rips,example_reduced_rips_from_distance_matrix.cpp} - The same square given
 * as a hand-typed (lower-triangular) distance matrix, with no coordinates.
 * <details>
 *   @dontinclude example_reduced_rips_from_distance_matrix.cpp
 *   @skip #include
 *   @until return 0;
 *   @skipline }
 * </details>
 *
 * @subsection reducedripsfromfile Reading a point cloud from a file
 *
 * \li \gudhi_example_link{Reduced_rips,example_reduced_rips_from_off.cpp} - Degree-1 persistence of a point
 * cloud read from an OFF file.
 *
 * @}
 */

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // DOC_REDUCED_RIPS_INTRO_REDUCED_RIPS_H_
