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

   The principal class of this module is Gudhi::Witness_complex.

   In both cases, the constructor for this class takes a {witness}x{closest_landmarks} table, where each row represents a witness and consists of landmarks sorted by distance to this witness.
   This table can be constructed by two additional classes Landmark_choice_by_furthest_point and Landmark_choice_by_random_point also included in the module.

   *\image html "bench_Cy8.png" "Running time as function on number of landmarks" width=10cm
   *\image html "bench_sphere.png" "Running time as function on number of witnesses for |L|=300" width=10cm
   

   \copyright GNU General Public License v3.


 */

#endif
