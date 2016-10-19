#ifndef WITNESS_COMPLEX_DOC_H_
#define WITNESS_COMPLEX_DOC_H_

/**
   \defgroup witness_complex Witness complex

   \author Siargey Kachanovich

   \image html "Witness_complex_representation.png" "Witness complex representation"

   \section witnessdefinitions Definitions

   Witness complex is a simplicial complex defined on two sets of points in \f$\mathbb{R}^D\f$:

   \li \f$W\f$ set of **witnesses** and 
   \li \f$L\f$ set of **landmarks**.

   Even though often the set of landmarks \f$L\f$ is a subset of the set of witnesses \f$ W\f$, it is not a requirement for the current implementation.

   Landmarks are the vertices of the simplicial complex
   and witnesses help to decide on which simplices are inserted via a predicate "is witnessed".

   De Silva and Carlsson in their paper \cite de2004topological differentiate **weak witnessing** and **strong witnessing**:

   - *weak*:   \f$ \sigma \subset L \f$ is witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in L \setminus \sigma,\ d(w,l) \leq d(w,l') \f$
   - *strong*: \f$ \sigma \subset L \f$ is witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in L,\ d(w,l) \leq d(w,l') \f$

   where \f$ d(.,.) \f$ is a distance function.

   Both definitions can be relaxed by a real value \f$\alpha\f$:

   - *weak*:   \f$ \sigma \subset L \f$ is \f$\alpha\f$-witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in L \setminus \sigma,\ d(w,l)^2 \leq d(w,l')^2 + \alpha^2 \f$
   - *strong*: \f$ \sigma \subset L \f$ is \f$\alpha\f$-witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in L,\ d(w,l)^2 \leq d(w,l')^2 + \alpha^2 \f$

   which leads to definitions of **weak relaxed witness complex** (or just relaxed witness complex for short) and **strong relaxed witness complex** respectively.

   \section witnessimplementation Implementation
   
   The two complexes described above are implemented in the corresponding classes 
   - Gudhi::witness_complex::Witness_complex
   - Gudhi::witness_complex::Strong_witness_complex

   The construction of both of them follow the same scheme:
   1. Construct a search tree on landmarks (for that Gudhi::spatial_searching::Kd_tree_search is used internally).
   2. Construct lists of nearest landmarks for each witness (special internal structure Gudhi::witness_complex::Active_witness is used internally).
   3. Construct the witness complex for nearest landmark lists.

   The constructors take on the step 1, while the function 'create_complex' executes the steps 2 and 3.

   \section witnessexample1 Example 1: Constructing weak relaxed witness complex from an off file

   Let's start with a simple example, which reads an off point file and computes a weak witness complex.

   \include Witness_complex/example_witness_complex_off.cpp

   \section witnessexample2 Example2: Computing persistence using strong relaxed witness complex

   Here is an example of constructing a strong witness complex filtration and computing persistence on it:

   \include Witness_complex/example_strong_witness_persistence.cpp

   *\image html "bench_Cy8.png" "Running time as function on number of landmarks" width=10cm
   *\image html "bench_sphere.png" "Running time as function on number of witnesses for |L|=300" width=10cm

   \copyright GNU General Public License v3.


 */

#endif  // WITNESS_COMPLEX_DOC_H_
