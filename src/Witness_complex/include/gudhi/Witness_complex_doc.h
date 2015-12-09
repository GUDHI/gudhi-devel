#ifndef WITNESS_COMPLEX_DOC_
#define WITNESS_COMPLEX_DOC_

/**
   \defgroup witness_complex Witness complex

   \author Siargey Kachanovich

   \section Definitions

   Witness complex \f$ Wit(W,L) \f$  is a simplicial complex defined on two sets of points in \f$\mathbb{R}^D\f$:

   \li \f$W\f$ set of **witnesses** and 
   \li \f$L \subseteq W\f$ set of **landmarks**.

   The simplices are based on landmarks
   and a simplex belongs to the witness complex if and only if it is witnessed, that is:

   \f$ \sigma \subset L \f$ is witnessed if there exists a point \f$w \in W\f$ such that 
   w is closer to the vertices of \f$ \sigma \f$ than other points in \f$ L \f$ and all of its faces are witnessed as well. 
      
   \section Implementation

   Two classes are implemented in this module: Gudhi::Witness_complex and Gudhi::Relaxed_witness_complex.

   While Gudhi::Witness_complex represents the classical witness complex, Gudhi::Relaxed_witness_complex takes an additional positive real parameter \f$ \alpha \f$ and constructs simplices \f$ \sigma \f$, for which
   there exists \f$ w \in W \f$, such that \f$ d(p,w) < d(q,w) + \alpha \f$ for all \f$ p \in \sigma, q \in L\setminus \sigma \f$.

   In both cases, the constructors take a {witness}x{closest_landmarks} table,
   which can be constructed by two additional classes Landmark_choice_by_furthest_point and Landmark_choice_by_random_point also included in the module.

   \copyright GNU General Public License v3.


 */

#endif
