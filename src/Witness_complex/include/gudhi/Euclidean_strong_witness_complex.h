/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL and Eigen3
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef EUCLIDEAN_STRONG_WITNESS_COMPLEX_H_
#define EUCLIDEAN_STRONG_WITNESS_COMPLEX_H_

#include <gudhi/Strong_witness_complex.h>
#include <gudhi/Active_witness/Active_witness.h>
#include <gudhi/Kd_tree_search.h>

#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <Eigen/src/Core/util/Macros.h>  // for EIGEN_VERSION_AT_LEAST

#include <utility>
#include <vector>

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Euclidean_strong_witness_complex is only available for CGAL >= 4.11
#endif

#if !EIGEN_VERSION_AT_LEAST(3,1,0)
# error Euclidean_strong_witness_complex is only available for Eigen3 >= 3.1.0 installed with CGAL
#endif

namespace Gudhi {

namespace witness_complex {

/**
 *  \private
 * \class Euclidean_strong_witness_complex
 * \brief Constructs strong witness complex for given sets of witnesses and landmarks in Euclidean space.
 * \ingroup witness_complex
 *
 * \tparam Kernel_ requires a <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class.
 */
template< class Kernel_ >
class Euclidean_strong_witness_complex
    : public Strong_witness_complex<std::vector<typename Gudhi::spatial_searching::Kd_tree_search<Kernel_,
                                                                                                  std::vector<typename Kernel_::Point_d>>::INS_range>> {
 private:
  typedef Kernel_                                                                      K;
  typedef typename K::Point_d                                                          Point_d;
  typedef std::vector<Point_d>                                                         Point_range;
  typedef Gudhi::spatial_searching::Kd_tree_search<Kernel_, Point_range>               Kd_tree;
  typedef typename Kd_tree::INS_range                                                  Nearest_landmark_range;
  typedef typename std::vector<Nearest_landmark_range>                                 Nearest_landmark_table;

  typedef typename Nearest_landmark_range::Point_with_transformed_distance Id_distance_pair;
  typedef typename Id_distance_pair::first_type Landmark_id;
  typedef Landmark_id Vertex_handle;

 private:
  Point_range                         landmarks_;
  Kd_tree                             landmark_tree_;
  using Strong_witness_complex<Nearest_landmark_table>::nearest_landmark_table_;

 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  /**
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details Records landmarks from the range 'landmarks' into a 
   *           table internally, as well as witnesses from the range 'witnesses'.
   *           Both ranges should have value_type Kernel_::Point_d.
   */
  template< typename LandmarkRange,
            typename WitnessRange >
  Euclidean_strong_witness_complex(const LandmarkRange & landmarks,
                                   const WitnessRange &  witnesses)
    : landmarks_(std::begin(landmarks), std::end(landmarks)), landmark_tree_(landmarks_) {
    nearest_landmark_table_.reserve(boost::size(witnesses));
    for (auto w : witnesses)
      nearest_landmark_table_.push_back(landmark_tree_.incremental_nearest_neighbors(w));
  }

  /** \brief Returns the point corresponding to the given vertex.
   */
  template <typename Vertex_handle>
  Point_d get_point(Vertex_handle vertex) const {
    return landmarks_[vertex];
  }

  //@}
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // EUCLIDEAN_STRONG_WITNESS_COMPLEX_H_
