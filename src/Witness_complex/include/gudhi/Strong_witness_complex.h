/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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

#ifndef STRONG_WITNESS_COMPLEX_H_
#define STRONG_WITNESS_COMPLEX_H_

#include <utility>
#include <vector>
#include <list>
#include <limits>

#include "Active_witness/Active_witness.h"
#include <gudhi/Kd_tree_search.h>

namespace gss = Gudhi::spatial_searching;

namespace Gudhi {
  
namespace witness_complex {

template< class Kernel_ >
class Strong_witness_complex {
private:
  typedef Kernel_                                                 K;
  typedef typename K::Point_d                                     Point_d;
  typedef typename K::FT                                          FT;
  typedef std::vector<Point_d>                                    Point_range;
  typedef gss::Kd_tree_search<Kernel_, Point_range>               Kd_tree;
  typedef typename Kd_tree::INS_range                             Nearest_landmark_range;
  typedef typename std::vector<Nearest_landmark_range>            Nearest_landmark_table;
  typedef typename Nearest_landmark_range::iterator               Nearest_landmark_row_iterator;
  
  typedef std::vector< double > Point_t;
  typedef std::vector< Point_t > Point_Vector;

  typedef FT Filtration_value;

  
  typedef std::ptrdiff_t Witness_id;
  typedef typename Nearest_landmark_range::Point_with_transformed_distance Id_distance_pair;
  typedef typename Id_distance_pair::first_type Landmark_id;
  typedef Active_witness<Id_distance_pair, Nearest_landmark_range> ActiveWitness;
  typedef std::list< ActiveWitness > ActiveWitnessList;
  typedef std::vector< Landmark_id > typeVectorVertex;
  typedef std::pair< typeVectorVertex, Filtration_value> typeSimplex;

 private:
  Point_range                         witnesses_, landmarks_;
  Kd_tree                             landmark_tree_;
  
 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  /*
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details The parameters should satisfy InputIterator C++ concepts.
   */
  template< typename InputIteratorLandmarks,
            typename InputIteratorWitnesses >
  Strong_witness_complex(InputIteratorLandmarks landmarks_first,
                  InputIteratorLandmarks landmarks_last,
                  InputIteratorWitnesses witnesses_first,
                  InputIteratorWitnesses witnesses_last)
    : witnesses_(witnesses_first, witnesses_last), landmarks_(landmarks_first, landmarks_last), landmark_tree_(landmarks_)
  {    
  }

  /** \brief Returns the point corresponding to the given vertex.
   */
  Point_d get_point( std::size_t vertex ) const
  {
    return landmarks_[vertex];
  }
  
  /** \brief Outputs the (strong) witness complex with 
   *         squared relaxation parameter 'max_alpha_square'
   *         to simplicial complex 'complex'.
   *         The parameter 'limit_dimension' represents the maximal dimension of the simplicial complex
   *         (default value = no limit).
   */
  template < typename SimplicialComplexForWitness >
  bool create_complex(SimplicialComplexForWitness& complex,
                      FT  max_alpha_square,
                      Landmark_id limit_dimension = std::numeric_limits<Landmark_id>::max())       
  {
    std::size_t nbL = landmarks_.size();
    Landmark_id complex_dim = 0;
    if (complex.num_vertices() > 0) {
      std::cerr << "Strong witness complex cannot create complex - complex is not empty.\n";
      return false;
    }
    if (max_alpha_square < 0) {
      std::cerr << "Strong witness complex cannot create complex - squared relaxation parameter must be non-negative.\n";
      return false;
    }
    if (limit_dimension < 0) {
      std::cerr << "Strong witness complex cannot create complex - limit dimension must be non-negative.\n";
      return false;
    }
    typeVectorVertex vv;
    for (unsigned i = 0; i != nbL; ++i) {
      // initial fill of 0-dimensional simplices
      vv = {i};
      complex.insert_simplex(vv, Filtration_value(0.0));
    }
    for (auto w: witnesses_) {
      ActiveWitness aw(landmark_tree_.query_incremental_nearest_neighbors(w));
      typeVectorVertex simplex;
      typename ActiveWitness::iterator aw_it = aw.begin();
      float lim_dist2 = aw.begin()->second + max_alpha_square;
      while ((Landmark_id)simplex.size() <= limit_dimension + 1 && aw_it != aw.end() && aw_it->second < lim_dist2) {
        simplex.push_back(aw_it->first);
        complex.insert_simplex_and_subfaces(simplex, aw_it->second - aw.begin()->second);
        aw_it++;
      }
      if ((Landmark_id)simplex.size() - 1 > complex_dim)
        complex_dim = simplex.size() - 1;
    }
    complex.set_dimension(complex_dim);
    return true;
  }

  //@}
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif
