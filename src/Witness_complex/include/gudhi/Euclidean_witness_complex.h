/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA (France)
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

#ifndef EUCLIDEAN_WITNESS_COMPLEX_H_
#define EUCLIDEAN_WITNESS_COMPLEX_H_

#include <gudhi/Witness_complex.h>
#include <gudhi/Active_witness/Active_witness.h>
#include <gudhi/Kd_tree_search.h>

#include <utility>
#include <vector>
#include <list>
#include <limits>

namespace Gudhi {

namespace witness_complex {

/**
 * \private
 * \class Euclidean_witness_complex
 * \brief Constructs (weak) witness complex for given sets of witnesses and landmarks in Euclidean space.
 * \ingroup witness_complex
 *
 * \tparam Kernel_ requires a <a target="_blank"
 * href="http://doc.cgal.org/latest/Kernel_d/classCGAL_1_1Epick__d.html">CGAL::Epick_d</a> class.
 */
template< class Kernel_ >
class Euclidean_witness_complex
    : public Witness_complex<std::vector<typename Gudhi::spatial_searching::Kd_tree_search<Kernel_,
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
  using Witness_complex<Nearest_landmark_table>::nearest_landmark_table_;

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
  Euclidean_witness_complex(const LandmarkRange & landmarks,
                            const WitnessRange &  witnesses)
    : landmarks_(std::begin(landmarks), std::end(landmarks)), landmark_tree_(landmarks) {
    nearest_landmark_table_.reserve(boost::size(witnesses));
    for (auto w : witnesses)
      nearest_landmark_table_.push_back(landmark_tree_.query_incremental_nearest_neighbors(w));
  }

  /** \brief Returns the point corresponding to the given vertex.
   *  @param[in] vertex Vertex handle of the point to retrieve.
   */
  Point_d get_point(Vertex_handle vertex) const {
    return landmarks_[vertex];
  }

  //@}
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // EUCLIDEAN_WITNESS_COMPLEX_H_
