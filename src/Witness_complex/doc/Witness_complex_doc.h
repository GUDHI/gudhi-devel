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

   - *weak*:   \f$ \sigma \subset L \f$ is witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in \mathbf{L \setminus \sigma},\ d(w,l) \leq d(w,l') \f$
   - *strong*: \f$ \sigma \subset L \f$ is witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in \mathbf{L},\ d(w,l) \leq d(w,l') \f$

   where \f$ d(.,.) \f$ is a distance function.

   Both definitions can be relaxed by a real value \f$\alpha\f$:

   - *weak*:   \f$ \sigma \subset L \f$ is \f$\alpha\f$-witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in \mathbf{L \setminus \sigma},\ d(w,l)^2 \leq d(w,l')^2 + \alpha^2 \f$
   - *strong*: \f$ \sigma \subset L \f$ is \f$\alpha\f$-witnessed by \f$ w \in W\f$ if \f$ \forall l \in \sigma,\ \forall l' \in \mathbf{L},\ d(w,l)^2 \leq d(w,l')^2 + \alpha^2 \f$

   which leads to definitions of **weak relaxed witness complex** (or just relaxed witness complex for short) and **strong relaxed witness complex** respectively.

   \image html "swit.svg" "Strongly witnessed simplex"

   In particular case of 0-relaxation, weak complex corresponds to **witness complex** introduced in \cite de2004topological, whereas 0-relaxed strong witness complex consists of just vertices and is not very interesting.
   Hence for small relaxation weak version is preferable.
   However, to capture the homotopy type (for example using Gudhi::persistent_cohomology::Persistent_cohomology) it is often necessary to work with higher filtration values. In this case strong relaxed witness complex is faster to compute and offers similar results.

   \section witnessimplementation Implementation
   
   The two complexes described above are implemented in the corresponding classes 
   - Gudhi::witness_complex::Witness_complex
   - Gudhi::witness_complex::Euclidean_witness_complex
   - Gudhi::witness_complex::Strong_witness_complex
   - Gudhi::witness_complex::Euclidean_strong_witness_complex

   The construction of the Euclidean versions of complexes follow the same scheme:
   1. Construct a search tree on landmarks (for that Gudhi::spatial_searching::Kd_tree_search is used internally).
   2. Construct lists of nearest landmarks for each witness (special structure Gudhi::witness_complex::Active_witness is used internally).
   3. Construct the witness complex for nearest landmark lists.

   In the non-Euclidean classes, the lists of nearest landmarks are supposed to be given as input.

   The constructors take on the steps 1 and 2, while the function 'create_complex' executes the step 3.

   \section witnessexample1 Example 1: Constructing weak relaxed witness complex from an off file

   Let's start with a simple example, which reads an off point file and computes a weak witness complex.

   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   
#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/Points_off_io.h>

#include <CGAL/Epick_d.h>

#include <string>
#include <vector>

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef typename K::Point_d Point_d;
typedef typename Gudhi::witness_complex::Euclidean_witness_complex<K> Witness_complex;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef std::vector< Point_d > Point_vector;

int main(int argc, char * const argv[]) {
  std::string file_name = argv[1];
  int nbL = atoi(argv[2]), lim_dim = atoi(argv[4]);
  double alpha2 = atof(argv[3]);
  Gudhi::Simplex_tree<> simplex_tree;

  // Read the point file
  Point_vector point_vector, landmarks;
  Gudhi::Points_off_reader<Point_d> off_reader(file_name);
  point_vector = Point_vector(off_reader.get_point_cloud());

  // Choose landmarks (one can choose either of the two methods below)
  // Gudhi::subsampling::pick_n_random_points(point_vector, nbL, std::back_inserter(landmarks));
  Gudhi::subsampling::choose_n_farthest_points(K().squared_distance_d_object(), point_vector, nbL,
      Gudhi::subsampling::random_starting_point, std::back_inserter(landmarks));

  // Compute witness complex
  Witness_complex witness_complex(landmarks, point_vector);

  witness_complex.create_complex(simplex_tree, alpha2, lim_dim);
}

   
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   \section witnessexample2 Example2: Computing persistence using strong relaxed witness complex

   Here is an example of constructing a strong witness complex filtration and computing persistence on it:

   \include Witness_complex/strong_witness_persistence.cpp

   \section witnessexample3 Example3: Computing relaxed witness complex persistence from a distance matrix 

   In this example we compute the relaxed witness complex persistence from a given matrix of closest landmarks to each witness.
   Each landmark is given as the couple (index, distance).

   \include Witness_complex/example_nearest_landmark_table.cpp

 */

#endif  // WITNESS_COMPLEX_DOC_H_
