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

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <math.h>
#include <ctime>
#include <iostream>

// Needed for nearest neighbours
#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

// Needed for the adjacency graph in bad link search
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

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

  // Witness_range<Closest_landmark_range<Vertex_handle>>

  /*
   *  \brief Iterative construction of the (weak) witness complex.
   *  \details The witness complex is written in sc_ basing on a matrix knn of
   *  nearest neighbours of the form {witnesses}x{landmarks}.
   *
   *  The type KNearestNeighbors can be seen as 
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges.
   *
   *  Constructor takes into account at most (dim+1) 
   *  first landmarks from each landmark range to construct simplices.
   *
   *  Landmarks are supposed to be in [0,nbL_-1]
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
  
  /** \brief Outputs the (weak) witness complex with 
   *         squared relaxation parameter 'max_alpha_square'
   *         to simplicial complex 'complex'.
   */
  template < typename SimplicialComplexForWitness >
  bool create_complex(SimplicialComplexForWitness& complex,
                      FT  max_alpha_square)       
  {
    unsigned nbL = landmarks_.size();
    if (complex.num_vertices() > 0) {
      std::cerr << "Witness complex cannot create complex - complex is not empty.\n";
      return false;
    }
    if (max_alpha_square < 0) {
      std::cerr << "Witness complex cannot create complex - squared relaxation parameter must be non-negative.\n";
      return false;
    }
    typeVectorVertex vv;
    //ActiveWitnessList active_witnesses;// = new ActiveWitnessList();
    for (unsigned i = 0; i != nbL; ++i) {
      // initial fill of 0-dimensional simplices
      vv = {i};
      complex.insert_simplex(vv, Filtration_value(0.0));
      /* TODO Error if not inserted : normally no need here though*/
    }
    for (auto w: witnesses_) {
      ActiveWitness aw(landmark_tree_.query_incremental_nearest_neighbors(w));
      typeVectorVertex simplex;
      typename ActiveWitness::iterator aw_it = aw.begin();
      float lim_d2 = aw.begin()->second + max_alpha_square;
      while (aw_it != aw.end() && aw_it->second < lim_d2) {
        simplex.push_back(aw_it->first);
        complex.insert_simplex_and_subfaces(simplex, aw_it->second - aw.begin()->second);
        aw_it++;
      }
    }
    return true;
  }

  //@}
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif
