/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
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
 * A vertex \f$v\f$ in a simplicial complex \f$K\f$ is called a \textbf{dominated vertex} if the link of \f$v\f$ in \f$K\f$, 
 * \f$lk_K(v)\f$ is a simplicial cone, that is, there exists a vertex \f$v^{\prime} \neq v\f$ and a subcomplex \f$L\f$ in \f$K\f$, 
 * such that \f$lk_K(v) = v^{\prime}L\f$. We say that the vertex  \f$v^{\prime}\f$ is {dominating} \f$v\f$ and \f$v\f$ is {dominated} by \f$v^{\prime}\f$. 
 * An \textbf{elementary strong collapse} is the deletion of a dominated vertex \f$v\f$ from \f$K\f$, 
 * which we denote with \f$K\f$ \f${\searrow\searrow}\f$ \f$K\setminus v\f$. 
 * The symbol \textbf{\f$K \setminus v\f$} (deletion of \f$v\f$ from \f$K\f$) refers to the subcomplex of \f$K\f$ which has all simplices of \f$K\f$ except the ones containing \f$v\f$.
 * There is a \textbf{strong collapse} from a simplicial complex \f$K\f$ to its subcomplex \f$L\f$, 
 * if there exists a series of elementary strong collapses from \f$K\f$ to \f$L\f$, denoted as \f$K\f$ \f${\searrow\searrow}\f$ \f$L\f$.
 * 
 * A strong collapse is a homotopy preserving operation,in fact it a special type of the classical simple collapse. 
 * A complex without any dominated vertex is called a minimal complex and the core \f$K^0\f$ of simplicial comlex is a minimal complex such that \f$K\f$ \f${\searrow\searrow}\f$ \f$K^0\f$.
 * Computation of a core involves computation of dominated vertices and the dominated vertices can be easily characterized as follows:
 * -- For general simplicial complex: A vertex \f$v \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$, \textit{if and only if} all the maximal simplices of \f$K\f$ that contain $v$ also contain \f$v^{\prime}\f$
 * -- For a flag complex: A vertex \f$v \in K\f$ is dominated by another vertex \f$v^{\prime} \in K\f$, \textit{if and only if} all the vertices in \f$K\f$ that has an edge with $v$ also has an edge with \f$v^{\prime}\f$, given 
 * \f$[v.v^\prime]\f$ is an egde \f$\in K\f$.
 *  
 * This module implements strong collapse of a flag complex, in particular  Vietoris-Rips (VR) complex from its graph (1-skeleton) based on the above domination criteria. The algorithm to compute the core is described in Section 3 \cite strongcollapsesocg2019. 
 * \subsection persitence_computation_via_strong_collapse Persitence computation via strong_collapse
 * Strong collapse can be successfully employed to reduce a given sequence of simplicial complexes to a smaller induced sequence which preserves the persistent homology of the original sequence. 
 * The general idea is to compute the cores of each complex in the sequence and there exists an induced simplicial maps among the cores, which preserves the persistent homology. 
 * The precise mechanism for this reduction has been described in Section 4 \cite strongcollapseesa2018. 
 * Here we implement this mechanism for a filtration of Rips complex, where a filtration is a sequence of simplicial (here Rips) complexes connected with inclusions. 
 * After perfoming the reduction the filtration reduces to a flag-tower, which is a sequence of flag complexes connected through more general (inclusion and contractions) simplicial maps. 
 * However to compute the persistent homology of any tower one needs to convert it to an equivalent filtration. Section 4 of \cite strongcollapsesocg2019 descibes an efficient strategy to convert a flag tower to a flag filtration, 
 * using only the 1-skeletons of the complexes, which we implement here.
 * 
 * \subsection approximate_and_exact_versions Approximate and Exact versions
 * Given a Rips filtration, one can choose to collapse the original complexes after each edge inclusion. 
 * However, we can also choose to strong collapse the complexes less often, i.e. after several edge inclusions rather than just one. 
 * This will result in a faster algorithm but comes with a cost: the computed PD is then only approximate. 
 * The input filtration of VR complexes associated to a set of increasing values of the scale parameter, we call \textbf{snapshots} the values of the scale parameter at which we choose to strong collapse the complex. 
 * The difference between two consecutive snapshots is called a \textbf{step}. We approximate the {filtration value} of a simplex as the value of the snapshot at which it first appears. 
 * This approach will report all persistence pairs that are separated by at least one snapshot. 
 * Hence if all steps are equal to some \f$\epsilon>0\f$, we will compute all the persistence pairs whose lengths are at least \f$\epsilon\f$. 
 * It follows that the  \textit{\f$l_{\infty}\f$-bottleneck distance} between the computed PD and the exact one is at most \f$\epsilon\f$. 
 * If instead  the ratio between any two consecutive steps is taken to a constant \f$\rho>1\f$, the \f$l_{\infty}\f$-bottleneck distance will be at most \f$\log \rho\f$  after reparameterizing the filtrations on a \f$\log\f$-\f$\log\f$-scale~\cite{sheehy13linear}. 

 * For more information about our approach of computing strong collapses and persitent homology via strong collapses, we refer the users to \cite strongcollapsesocg2019 and \cite strongcollapseesa2018.
 * %\subsection strong_collapse_latex Some LaTeX
 * %Let \f$K\f$ be a simplicial complex.
 *
 * %Or like that :
 * \f[\tau \subseteq \gamma \subseteq \sigma \f]
 *
 * \subsection strong_collapse_doxygen Some Doxygen special commands
 *
 * @warning Please, try to do an "example driven documentation".
 *
 * @note Just a note.
 *
 * <a target="_blank" href="http://www.doxygen.nl/manual/commands.html">Doxygen special commands</a>
 *
 * A link to an image from the doc repository
 * \image html "strong_collapse_representation.png" "Strong collapse representation"
 * 
 * \subsection strong_collapse_citations Link to GUDHI biblio
 *
 * A reference \cite sheehy13linear.
 * Modify biblio/bibliography.bib to add new references.
 *
 * \subsection strong_collapse_code_doc Link to doxygen from code
 *
 * A link to `Flag_complex_sparse_matrix` auto documented class.
 * cf. src/Strong_collapse/include/gudhi/strong_collapse/Flag_complex_sparse_matrix.h
 *
 * \subsection strong_collapse_code Some code
 *
 * Example : We provide some use case examples in the following code.
 *
 * 
 * \include Strong_collapse/strong_collapse_rips_persistence_step_by_step.cpp
 * 
 *
 * 
 */
/** @} */  // end defgroup rips_complex

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // DOC_STRONG_COLLAPSE_INTRO_STRONG_COLLAPSE_H_
