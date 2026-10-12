/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber, David Loiseaux
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_DOC_MULTI_FILTRATION_INTRO_H_
#define GUDHI_DOC_MULTI_FILTRATION_INTRO_H_

// needs namespace for Doxygen to link on classes
namespace Gudhi {
namespace multi_filtration {

/** \defgroup multi_filtration Multi-parameter Filtration Values
 * @{
 * \author    Hannah Schreiber, David Loiseaux
 *
 * \section multifiltrationintro Multi-parameter Filtration Values
 *
 * A multi-parameter filtration is a family of spaces or complexes indexed by vectors in \f$ \mathbb R^n \f$.
 * The index set is ordered coordinatewise: for \f$ u = (u_1,\ldots,u_n) \f$ and \f$ v = (v_1,\ldots,v_n) \f$,
 * we write \f$ u \preceq v \f$ when \f$ u_i \leq v_i \f$ for every \f$ i\f$. A filtration
 * \f$ (K_u)_{u\in\mathbb R^n} \f$ satisfies \f$ K_u \subseteq K_v \f$ whenever \f$ u \preceq v \f$.
 *
 * For a cell (for example, a simplex) \f$ \sigma \f$, the set of parameter values at which it exists is therefore
 * an upward-closed subset of \f$ \mathbb R^n \f$. This module represents such sets when they are finitely generated:
 * a finite set of vectors \f$ G = \{g_1,\ldots,g_k\} \subseteq \mathbb R^n \f$ represents the union of positive cones
 * \f[
 *   L(G)=\bigcup_{g\in G}\{u\in\mathbb R^n\mid g\preceq u\}.
 * \f]
 * In other words, the cell exists at \f$ u \f$ if at least one generator \f$ g \in G \f$ satisfies
 * \f$ g_i \leq u_i \f$ for every parameter \f$ i \f$. Each generator is an **apparition/birth time** (also called a
 * **corner**) for the cell. It is a vector of \f$ n \f$ values, one for each parameter. Note that it is not a
 * generator of a homology module in the algebraic sense.
 *
 * We say that the set of generators in the filtration value of a cell is minimal if no generators is covered by
 * the cone of another generator in the set. I.e., all generators are incomparable. If the represented set has
 * \f$ k \f$ minimal generators, its value is called \f$ k \f$-critical.
 *
 * A filtration value can thus be viewed either as its set of corners or, equivalently, as the upward-closed set
 * they generate. The order on filtration values is defined by inclusion of these sets: `a <= b` means that the set
 * represented by `b` is contained in the set represented by `a`. This convention agrees with the usual order for
 * 1-parameter filtration values: a cell appearing earlier has a smaller value. To define a valid filtration, the set
 * of parameter values at which a coface is present must be contained in the corresponding set for each of its faces.
 *
 * With more than one parameter, values are in general **not totally ordered**. For instance, the values with the
 * single corners \f$(1,3)\f$ and \f$(3,1)\f$ are incomparable: neither of their generated quadrants contains the other.
 * The relational operators `operator<`, `operator<=`, and their opposites use this cone-inclusion order. When an
 * application needs a total order merely to arrange values, the explicit lexicographical comparison functions
 * (such as @ref
 * Multi_parameter_filtration_value::is_less_or_equal_than_lexicographically "is_less_or_equal_than_lexicographically")
 * should be used instead. The lexicographical order is not the filtration order.
 *
 * \subsection multifiltrationinterface Implementation
 *
 * The main class is @ref Multi_parameter_filtration_value. It stores the corners of a finitely generated
 * multi-parameter filtration value of a cell and provides operations to inspect and modify them, compare values by
 * cone inclusion, simplify a set of corners, perform arithmetic operations, convert between storage representations,
 * and handle special values such as infinities and NaNs. It is a value representation utility: it does not, by itself,
 * compute persistent homology or a multi-parameter persistence module. It can instead be used by data structures and
 * algorithms that need to store multi-parameter filtration values. It is usable as filtration value in a
 * @ref Gudhi::Simplex_tree at it realizes the @ref FiltrationValue concept.
 *
 * The first template argument of @ref Multi_parameter_filtration_value is a storage
 * policy. The module currently provides:
 * - @ref Flat_array_filtration, which stores the coordinates of all corners in one flat
 *   array. Its representation is assumed to be minimal and lexicographically ordered.
 * - @ref Nested_array_filtration, which stores the corners in a vector of vectors and
 *   likewise assumes a minimal, ordered representation.
 * - @ref Degree_bifiltration, a specialized compact representation for a two-parameter
 *   filtration in which one coordinate follows a regular `shift + index * step` mapping.
 *
 * `Co` is a compile-time option that reverses the order and uses negative rather than positive cones.
 * `Ensure1Criticality` is another compile-time option that constrains a value to remain 1-critical.
 * The storage policy determines how the corners are stored and which representation-specific operations are available.
 *
 * \subsection multifiltrationdegreerips Implicit axis
 *
 * The storage policy @ref Degree_bifiltration only stores the first parameter explicitly and not the second parameter.
 * This specialization is useful when a bifiltration has a predictable sequence of corners along the implicit parameter.
 * For example, a Degree-Rips bifiltration combines the Rips scale with a degree-based parameter. For a value with
 * corners \f$ [(r_0,s), (r_1,s+\Delta), \ldots, (r_{k-1}, s+(k-1)\Delta)] \f$, the second coordinate is determined by
 * the corner index, so @ref Degree_bifiltration stores only the first coordinates \f$ r_i \f$ and reconstructs the
 * second coordinates from `shift` and `step` (with \f$ s = \texttt{shift} \f$ and \f$ \Delta = \texttt{step} \f$).
 * This avoids storing both coordinates for every corner and can save memory when such values are attached to many
 * simplices. An applied example is available in the
 * [Degree Rips bifiltration notebook](https://davidlapous.github.io/multipers/notebooks/degree_rips_interface.html).
 *
 * \section multifiltrationexamples Examples
 *
 * \subsection mfminusage Minimalistic example
 *
 * The following example creates a value with two corners, adds a second value with one corner, and illustrates both
 * comparability and incomparability. The two corners of \f$ f \f$ generate upper-right quadrants rooted at
 * \f$ (1,3) \f$ and \f$(3,1)\f$. The quadrant rooted at \f$ (2,2) \f$ is not contained in their union, and neither
 * is their union contained in that quadrant.
 *
 * \code{.cpp}
 * #include <cassert>
 * #include <iostream>
 * #include <vector>
 *
 * #include <gudhi/Multi_parameter_filtration_value.h>
 * #include <gudhi/Multi_filtration/Flat_array_filtration.h>
 *
 * using Storage_policy = Gudhi::multi_filtration::Flat_array_filtration<double>;
 * using Filtration_value = Gudhi::multi_filtration::Multi_parameter_filtration_value<Storage_policy>;
 *
 * int main() {
 *   // One generator (corner) with two parameters.
 *   Filtration_value f({1.0, 3.0});
 *
 *   // Add another corner. The flat-array policy keeps the minimal set ordered.
 *   const std::vector<double> second_corner{3.0, 1.0};
 *   f.add_generator(second_corner);
 *
 *   Filtration_value g({2.0, 2.0});
 *   assert(f.num_parameters() == 2);
 *   assert(f.num_generators() == 2);
 *   assert(!(f <= g));
 *   assert(!(g <= f));  // The two values are incomparable.
 *
 *   Filtration_value h({4.0, 4.0});
 *   assert(g <= h);     // The quadrant rooted at (4, 4) is contained in that rooted at (2, 2).
 *
 *   std::cout << "f = " << f << '\n';
 *   std::cout << "f has " << f.num_generators() << " corners and "
 *             << f.num_parameters() << " parameters.\n";
 * }
 * \endcode
 *
 * The call to @ref Multi_parameter_filtration_value::add_generator preserves a sorted representation and, for
 * storage policies that support a minimal-set representation, discards redundant corners. If corners are inserted
 * through lower-level operations that do not check these properties, call
 * @ref Multi_parameter_filtration_value::simplify before relying on methods that require a minimal, ordered set.
 *
 * \subsection mfother Other examples
 *
 * \li \gudhi_example_link{Multi_filtration,example_multi_simplex_tree.cpp} - A simple example
 * to showcase how to use the @ref Multi_parameter_filtration_value class with a @ref Gudhi::Simplex_tree.
 *
 * @}
 */

}  // namespace multi_filtration
}  // namespace Gudhi

#endif  // GUDHI_DOC_MULTI_FILTRATION_INTRO_H_
