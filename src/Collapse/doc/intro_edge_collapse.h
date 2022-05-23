/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_EDGE_COLLAPSE_INTRO_EDGE_COLLAPSE_H_
#define DOC_EDGE_COLLAPSE_INTRO_EDGE_COLLAPSE_H_

namespace Gudhi {

namespace collapse {

/**  \defgroup edge_collapse Edge collapse
 * 
 * \author    Siddharth Pritam and Marc Glisse
 * 
 * @{
 * 
 * This module implements edge collapse of a filtered flag complex as described in \cite edgecollapsearxiv, in
 * particular it reduces a filtration of Vietoris-Rips complex represented by a graph to a smaller flag filtration with
 * the same persistent homology.
 * 
 * \section edge_collapse_definition Edge collapse definition
 * 
 * An edge \f$e\f$ in a simplicial complex \f$K\f$ is called a <b>dominated edge</b> if the link of \f$e\f$ in
 * \f$K\f$, \f$lk_K(e)\f$ is a simplicial cone, that is, there exists a vertex \f$v^{\prime} \notin e\f$ and a
 * subcomplex \f$L\f$ in \f$K\f$, such that \f$lk_K(e) = v^{\prime}L\f$. We say that the vertex \f$v^{\prime}\f$
 * \e dominates \f$e\f$ and \f$e\f$ is \e dominated by \f$v^{\prime}\f$.
 * An <b> elementary edge collapse </b> is the removal of a dominated edge \f$e\f$ from \f$K\f$ (the cofaces of \f$e\f$
 * are implicitly removed as well).
 * Domination is used as a simple sufficient condition that ensures that this removal is a homotopy preserving
 * operation.
 *
 * The dominated edges can be easily characterized as follows:
 * 
 * -- For a general simplicial complex: an edge \f$e \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$,
 * if and only if all the maximal simplices of \f$K\f$ that contain \f$e\f$ also contain \f$v^{\prime}\f$.
 * 
 * -- For a flag complex: an edge \f$e \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$, if and only
 * if all the vertices in \f$K\f$ that have an edge with both vertices of \f$e\f$ also have an edge with
 * \f$v^{\prime}\f$. Notice that this only depends on the graph.
 * 
 * In the context of a filtration, an edge collapse may translate into an increase of the filtration value of an edge,
 * or its removal if it already had the largest filtration value.
 * The algorithm to compute the smaller induced filtration is described in \cite edgecollapsearxiv.
 * Edge collapse can be successfully employed to reduce any input filtration of flag complexes to a smaller induced
 * filtration which preserves the persistent homology of the original filtration and is a flag complex as well.
 * 
 * The algorithm implemented here does not produce a minimal filtration. Taking its output and applying the algorithm a
 * second time may further simplify the filtration.
 *
 * \subsection edgecollapseexample Basic edge collapse
 * 
 * This example calls `Gudhi::collapse::flag_complex_collapse_edges()` from a proximity graph represented as a list of
 * `Filtered_edge`.
 * Then it collapses edges and displays a new list of `Filtered_edge` (with fewer edges)
 * that will preserve the persistence homology computation.
 * 
 * \include edge_collapse_basic_example.cpp
 * 
 * When launching the example:
 * 
 * \code $> ./Edge_collapse_example_basic
 * \endcode
 *
 * the program output could be:
 * 
 * \include edge_collapse_example_basic.txt
 */
/** @} */  // end defgroup strong_collapse

}  // namespace collapse

}  // namespace Gudhi

#endif  // DOC_EDGE_COLLAPSE_INTRO_EDGE_COLLAPSE_H_
