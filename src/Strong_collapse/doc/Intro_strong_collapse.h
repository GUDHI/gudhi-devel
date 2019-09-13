/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_STRONG_COLLAPSE_INTRO_STRONG_COLLAPSE_H_
#define DOC_STRONG_COLLAPSE_INTRO_STRONG_COLLAPSE_H_

namespace Gudhi {

namespace strong_collapse {

/**  \defgroup strong_collapse Strong collapse
 * 
 * \author    Siddharth Pritam
 * 
 * @{
 * 
 * \section strong_collapse_definition Strong collapse definition
 * 
 * A vertex \f$v\f$ in a simplicial complex \f$K\f$ is called a <b>dominated vertex</b> if the link of \f$v\f$ in
 * \f$K\f$, \f$lk_K(v)\f$ is a simplicial cone, that is, there exists a vertex \f$v^{\prime} \neq v\f$ and a subcomplex
 * \f$L\f$ in \f$K\f$, such that \f$lk_K(v) = v^{\prime}L\f$. We say that the vertex  \f$v^{\prime}\f$ is {dominating}
 * \f$v\f$ and \f$v\f$ is {dominated} by \f$v^{\prime}\f$. 
 * An <b>elementary strong collapse</b> is the deletion of a dominated vertex \f$v\f$ from \f$K\f$, 
 * which we denote with \f$K\f$ \f${\searrow\searrow}\f$ \f$K\setminus v\f$. 
 * The symbol \f$\mathbf{K\setminus v}\f$ (deletion of \f$v\f$ from \f$K\f$) refers to the subcomplex of \f$K\f$ which
 * has all simplices of \f$K\f$ except the ones containing \f$v\f$.
 * There is a <b>strong collapse</b> from a simplicial complex \f$K\f$ to its subcomplex \f$L\f$, 
 * if there exists a series of elementary strong collapses from \f$K\f$ to \f$L\f$, denoted as \f$K\f$
 * \f${\searrow\searrow}\f$ \f$L\f$.
 * 
 * A strong collapse is a homotopy preserving operation,in fact it a special type of the classical simple collapse. 
 * A complex without any dominated vertex is called a minimal complex and the core \f$K^0\f$ of simplicial comlex is a
 * minimal complex such that \f$K\f$ \f${\searrow\searrow}\f$ \f$K^0\f$.
 * Computation of a core involves computation of dominated vertices and the dominated vertices can be easily
 * characterized as follows:
 * 
 * -- For general simplicial complex: A vertex \f$v \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$,
 * <i>if and only if</i> all the maximal simplices of \f$K\f$ that contain $v$ also contain \f$v^{\prime}\f$
 * 
 * -- For a flag complex: A vertex \f$v \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$, <i>if and only
 * if</i> all the vertices in \f$K\f$ that has an edge with $v$ also has an edge with \f$v^{\prime}\f$, given 
 * \f$[v.v^\prime]\f$ is an egde \f$\in K\f$.
 *  
 * This module implements strong collapse of a flag complex, in particular  Vietoris-Rips (VR) complex from its graph
 * (1-skeleton) based on the above domination criteria. The algorithm to compute the core is described in Section 3
 * \cite strongcollapsesocg2019. 
 * \subsection persitence_computation_via_strong_collapse Persitence computation via strong_collapse
 * Strong collapse can be successfully employed to reduce a given sequence of simplicial complexes to a smaller induced
 * sequence which preserves the persistent homology of the original sequence. 
 * The general idea is to compute the cores of each complex in the sequence and there exists an induced simplicial maps
 * among the cores, which preserves the persistent homology. 
 * The precise mechanism for this reduction has been described in Section 4 \cite strongcollapseesa2018. 
 * Here we implement this mechanism for a filtration of Rips complex, where a filtration is a sequence of simplicial
 * (here Rips) complexes connected with inclusions. 
 * After perfoming the reduction the filtration reduces to a flag-tower, which is a sequence of flag complexes
 * connected through more general (inclusion and contractions) simplicial maps. 
 * However to compute the persistent homology of any tower one needs to convert it to an equivalent filtration. Section
 * 4 of \cite strongcollapsesocg2019 descibes an efficient strategy to convert a flag tower to a flag filtration, 
 * using only the 1-skeletons of the complexes, which we implement here.
 * 
 * \subsection strong_collapse_from_points_example Example from a point cloud and a distance function
 * 
 * This example builds the edge graph from the given points, threshold value, and distance function.
 * Then it creates a `Flag_complex_strong_collapse` (exact version) with it.
 * 
 * Then, it is asked to display the distance matrix after the collapse operation.
 * 
 * \include Strong_collapse/strong_collapse_from_points.cpp
 * 
 * \code $> ./strong_collapse_from_points
 * \endcode
 *
 * the program output is:
 * 
 * \include Strong_collapse/strong_collapse_from_points_for_doc.txt
 * 
 * A `Gudhi::rips_complex::Rips_complex` can be built from the distance matrix if you want to compute persistence on
 * top of it.
 * 
 * \subsection strong_collapse_approximate_and_exact_versions Approximate and Exact versions
 * Given a Rips filtration, one can choose to collapse the original complexes after each edge inclusion. 
 * However, we can also choose to strong collapse the complexes less often, i.e. after several edge inclusions rather
 * than just one. 
 * This will result in a faster algorithm but comes with a cost: the computed PD is then only approximate. 
 * The input filtration of VR complexes associated to a set of increasing values of the scale parameter, we call
 * <b>snapshots</b> the values of the scale parameter at which we choose to strong collapse the complex. 
 * The difference between two consecutive snapshots is called a <b>step</b>. We approximate the {filtration value} of a
 * simplex as the value of the snapshot at which it first appears. 
 * This approach will report all persistence pairs that are separated by at least one snapshot. 
 * Hence if all steps are equal to some \f$\epsilon>0\f$, we will compute all the persistence pairs whose lengths are
 * at least \f$\epsilon\f$. 
 * It follows that the <i>\f$l_{\infty}\f$-bottleneck distance</i> between the computed PD and the exact one is at most
 * \f$\epsilon\f$. 
 * If instead  the ratio between any two consecutive steps is taken to a constant \f$\rho>1\f$, the
 * \f$l_{\infty}\f$-bottleneck distance will be at most \f$\log \rho\f$  after reparameterizing the filtrations on a
 * \f$\log\f$-\f$\log\f$-scale~ \cite sheehy13linear. 

 * For more information about our approach of computing strong collapses and persitent homology via strong collapses,
 * we refer the users to \cite strongcollapsesocg2019 and \cite strongcollapseesa2018.
 * 
 */
/** @} */  // end defgroup strong_collapse

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // DOC_STRONG_COLLAPSE_INTRO_STRONG_COLLAPSE_H_
