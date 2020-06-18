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
 * \author    Siddharth Pritam
 * 
 * @{
 * 
 * This module implements edge collapse of a filtered flag complex, in particular it reduces a filtration of
 * Vietoris-Rips complex from its graph to another smaller flag filtration with same persistence.
 * Where a filtration is a sequence of simplicial (here Rips) complexes connected with inclusions.
 * 
 * \section edge_collapse_definition Edge collapse definition
 * 
 * An edge \f$e\f$ in a simplicial complex \f$K\f$ is called a <b>dominated edge</b> if the link of \f$e\f$ in
 * \f$K\f$, \f$lk_K(e)\f$ is a simplicial cone, that is, there exists a vertex \f$v^{\prime} \notin e\f$ and a
 * subcomplex \f$L\f$ in \f$K\f$, such that \f$lk_K(e) = v^{\prime}L\f$. We say that the vertex  \f$v^{\prime}\f$ is
 * {dominating} \f$e\f$ and \f$e\f$ is {dominated} by \f$v^{\prime}\f$. 
 * An <b> elementary egde collapse </b> is the removal of a dominated edge \f$e\f$ from \f$K\f$, 
 * which we denote with \f$K\f$ \f${\searrow\searrow}^1 \f$ \f$K\setminus e\f$. 
 * The symbol \f$\mathbf{K\setminus e}\f$ (deletion of \f$e\f$ from \f$K\f$) refers to the subcomplex of \f$K\f$ which
 * has all simplices of \f$K\f$ except \f$e\f$ and the ones containing \f$e\f$.
 * There is an <b>edge collapse</b> from a simplicial complex \f$K\f$ to its subcomplex \f$L\f$, 
 * if there exists a series of elementary edge collapses from \f$K\f$ to \f$L\f$, denoted as \f$K\f$
 * \f${\searrow\searrow}\f$ \f$L\f$.
 * 
 * An edge collapse is a homotopy preserving operation, and it can be further expressed as sequence of the classical
 * elementary simple collapse. 
 * A complex without any dominated edge is called a \f$1\f$- minimal complex and the core \f$K^1\f$ of simplicial
 * complex is a minimal complex such that \f$K\f$ \f${\searrow\searrow}\f$ \f$K^1\f$.
 * Computation of a core (not unique) involves computation of dominated edges and the dominated edges can be easily
 * characterized as follows:
 * 
 * -- For general simplicial complex: An edge \f$e \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$,
 * <i>if and only if</i> all the maximal simplices of \f$K\f$ that contain \f$e\f$ also contain \f$v^{\prime}\f$
 * 
 * -- For a flag complex: An edge \f$e \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$, <i>if and only
 * if</i> all the vertices in \f$K\f$ that has an edge with both vertices of \f$e\f$  also has an edge with
 * \f$v^{\prime}\f$.
 * 
 * The algorithm to compute the smaller induced filtration is described in Section 5 \cite edgecollapsesocg2020.
 * Edge collapse can be successfully employed to reduce any given filtration of flag complexes to a smaller induced
 * filtration which preserves the persistent homology of the original filtration and is a flag complex as well.
 * 
 * The general idea is that we consider edges in the filtered graph and sort them according to their filtration value
 * giving them a total order.
 * Each edge gets a unique index denoted as \f$i\f$ in this order.  To reduce the filtration, we move forward with
 * increasing filtration value 
 * in the graph and check if the current edge \f$e_i\f$ is dominated in the current graph \f$G_i := \{e_1, .. e_i\} \f$
 * or not. 
 * If the edge \f$e_i\f$ is dominated we remove it from the filtration and move forward to the next edge \f$e_{i+1}\f$.
 * If \f$e_i\f$ is non-dominated then we keep it in the reduced filtration and then go backward in the current graph
 * \f$G_i\f$ to look for new non-dominated edges that was dominated before but might become non-dominated at this
 * point. 
 * If an edge \f$e_j, j < i \f$ during the backward search is found to be non-dominated, we include \f$e_j\f$ in to the
 * reduced filtration and we set its new filtration value to be \f$i\f$ that is the index of \f$e_i\f$.
 * The precise mechanism for this reduction has been described in Section 5 \cite edgecollapsesocg2020. 
 * Here we implement this mechanism for a filtration of Rips complex.
 * After perfoming the reduction the filtration reduces to a flag-filtration with the same persistence as the original
 * filtration. 
 * 
 * \subsection edgecollapseexample Basic edge collapse
 * 
 * This example builds the `Flag_complex_edge_collapser` from a proximity graph represented as a list of
 * `Flag_complex_edge_collapser::Filtered_edge`.
 * Then it collapses edges and displays a new list of `Flag_complex_edge_collapser::Filtered_edge` (with less edges)
 * that will preserve the persistence homology computation.
 * 
 * \include Collapse/edge_collapse_basic_example.cpp
 * 
 * When launching the example:
 * 
 * \code $> ./Edge_collapse_example_basic
 * \endcode
 *
 * the program output is:
 * 
 * \include Collapse/edge_collapse_example_basic.txt
 */
/** @} */  // end defgroup strong_collapse

}  // namespace collapse

}  // namespace Gudhi

#endif  // DOC_EDGE_COLLAPSE_INTRO_EDGE_COLLAPSE_H_
