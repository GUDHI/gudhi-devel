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

#ifndef WITNESS_COMPLEX_H_
#define WITNESS_COMPLEX_H_

// Needed for the adjacency graph in bad link search
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/range/size.hpp>

#include "Active_witness/Active_witness.h"
#include <gudhi/distance_functions.h>
#include <gudhi/Kd_tree_search.h>

#include <algorithm>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <ctime>
#include <iostream>

namespace gss = Gudhi::spatial_searching;

namespace Gudhi {
  
namespace witness_complex {

// /*
//  *  \private
//     \class Witness_complex
//     \brief Constructs the witness complex for the given set of witnesses and landmarks.
//     \ingroup witness_complex
//  */
template< class Kernel_ >
class Witness_complex {
private:
  typedef Kernel_                                                 K;
  typedef typename K::Point_d                                     Point_d;
  typedef typename K::FT                                          FT;
  typedef std::vector<Point_d>                                    Point_range;
  typedef gss::Kd_tree_search<Kernel_, Point_range>               Kd_tree;
  typedef typename Kd_tree::INS_range                             Nearest_landmark_range;
  typedef typename std::vector<Nearest_landmark_range>            Nearest_landmark_table;
  typedef typename Nearest_landmark_range::iterator               Nearest_landmark_row_iterator;
  //typedef std::vector<std::pair<unsigned,FT>>                     Nearest_landmarks;
  
  // struct Active_witness {
  //   int witness_id;
  //   int landmark_id;

  //   Active_witness(int witness_id_, int landmark_id_)
  //       : witness_id(witness_id_),
  //       landmark_id(landmark_id_) { }
  // };
  
  typedef std::vector< double > Point_t;
  typedef std::vector< Point_t > Point_Vector;

  typedef FT Filtration_value;

  
  typedef std::ptrdiff_t Witness_id;
  typedef std::ptrdiff_t Landmark_id;
  typedef std::pair<Landmark_id, FT> Id_distance_pair;
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
  Witness_complex(InputIteratorLandmarks landmarks_first,
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
    ActiveWitnessList active_witnesses;// = new ActiveWitnessList();
    for (unsigned i = 0; i != nbL; ++i) {
      // initial fill of 0-dimensional simplices
      // by doing it we don't assume that landmarks are necessarily witnesses themselves anymore
      //counter++;
      vv = {i};
      complex.insert_simplex(vv, Filtration_value(0.0));
      /* TODO Error if not inserted : normally no need here though*/
    }
    unsigned k = 1; /* current dimension in iterative construction */
    for (auto w: witnesses_)
      active_witnesses.push_back(ActiveWitness(landmark_tree_.query_incremental_nearest_neighbors(w)));
    ActiveWitness aw_copy(active_witnesses.front());
    while (!active_witnesses.empty() && k < nbL ) {
      typename ActiveWitnessList::iterator aw_it = active_witnesses.begin();
      while (aw_it != active_witnesses.end()) {
        std::vector<int> simplex;
        //simplex.reserve(k+1);
        bool ok = add_all_faces_of_dimension(k,
                                             max_alpha_square,
                                             std::numeric_limits<double>::infinity(),
                                             aw_it->begin(),
                                             simplex,
                                             complex,
                                             aw_it->end());
        if (!ok)
          //{aw_it++;}
        active_witnesses.erase(aw_it++); //First increase the iterator and then erase the previous element
        else
          aw_it++;
      } 
      std::cout << "Active witnesses after dim=" << k << " is finished: " << active_witnesses.size() << "\n";
      std::cout << complex << "\n";
      k++;
    }
    return true;
  }

  //@}

 private:
  /* \brief Adds recursively all the faces of a certain dimension dim witnessed by the same witness
   * Iterator is needed to know until how far we can take landmarks to form simplexes
   * simplex is the prefix of the simplexes to insert
   * The output value indicates if the witness rests active or not
   */
  template < typename SimplicialComplexForWitness >
  bool add_all_faces_of_dimension(int dim,
                                  double alpha2,
                                  double norelax_dist2,
                                  typename ActiveWitness::iterator curr_l,
                                  std::vector<int>& simplex,
                                  SimplicialComplexForWitness& sc,
                                  typename ActiveWitness::iterator end)
  {
    if (curr_l == end)
      return false;
    bool will_be_active = false;
    typename ActiveWitness::iterator l_it = curr_l;
    if (dim > 0)
      for (; l_it->second - alpha2 <= norelax_dist2 && l_it != end; ++l_it) {
        simplex.push_back(l_it->first);
        if (sc.find(simplex) != sc.null_simplex()) {
          typename ActiveWitness::iterator next_it = l_it;
          will_be_active = add_all_faces_of_dimension(dim-1,
                                                      alpha2,
                                                      norelax_dist2,
                                                      ++next_it,
                                                      simplex,
                                                      sc,
                                                      end) || will_be_active;
        }
        assert(!simplex.empty());
        simplex.pop_back();
        // If norelax_dist is infinity, change to first omitted distance
        if (l_it->second <= norelax_dist2)
          norelax_dist2 = l_it->second;
        typename ActiveWitness::iterator next_it = l_it;
        will_be_active = add_all_faces_of_dimension(dim,
                                                    alpha2,
                                                    norelax_dist2,
                                                    ++next_it,
                                                    simplex,
                                                    sc,
                                                    end) || will_be_active;
      } 
    else if (dim == 0)
      for (; l_it->second - alpha2 <= norelax_dist2 && l_it != end; ++l_it) {
        simplex.push_back(l_it->first);
        double filtration_value = 0;
        // if norelax_dist is infinite, relaxation is 0.
        //std::cout << "landmark_id=" << l_it->first << " distance=" << l_it->second << "\n";
        // std::size_t landmark_id = l_it->first;
        // double distance = l_it->second;
        if (l_it->second > norelax_dist2) 
          filtration_value = l_it->second - norelax_dist2; 
        if (all_faces_in(simplex, &filtration_value, sc)) {
          will_be_active = true;
          sc.insert_simplex(simplex, filtration_value);
        }
        assert(!simplex.empty());
        simplex.pop_back();
        // If norelax_dist is infinity, change to first omitted distance
        if (l_it->second < norelax_dist2)
          norelax_dist2 = l_it->second;
      } 
    return will_be_active;
  }
  
  /** \brief Check if the facets of the k-dimensional simplex witnessed 
   *  by witness witness_id are already in the complex.
   *  inserted_vertex is the handle of the (k+1)-th vertex witnessed by witness_id
   */
  template < typename SimplicialComplexForWitness >
  bool all_faces_in(std::vector<int>& simplex,
                    double* filtration_value,
                    SimplicialComplexForWitness& sc)
  {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;

    std::vector< int > facet;
    for (std::vector<int>::iterator not_it = simplex.begin(); not_it != simplex.end(); ++not_it)
      {
        facet.clear();
        for (std::vector<int>::iterator it = simplex.begin(); it != simplex.end(); ++it)
          if (it != not_it)
            facet.push_back(*it);
        Simplex_handle facet_sh = sc.find(facet);
        if (facet_sh == sc.null_simplex())
          return false;
        else if (sc.filtration(facet_sh) > *filtration_value)
          *filtration_value = sc.filtration(facet_sh);
      }
    return true;
  }

  // bool is_face(Simplex_handle face, Simplex_handle coface)
  // {
  //   // vertex range is sorted in decreasing order
  //   auto fvr = sc.simplex_vertex_range(face);
  //   auto cfvr = sc.simplex_vertex_range(coface);
  //   auto fv_it = fvr.begin();
  //   auto cfv_it = cfvr.begin();
  //   while (fv_it != fvr.end() && cfv_it != cfvr.end()) {
  //     if (*fv_it < *cfv_it)
  //       ++cfv_it;
  //     else if (*fv_it == *cfv_it) {
  //       ++cfv_it;
  //       ++fv_it;
  //     }
  //     else
  //       return false;
      
  //   }
  //   return (fv_it == fvr.end());
  // }

  
 public:
  template < typename SimplicialComplexForWitness >
  void print_complex(SimplicialComplexForWitness& complex)
  {
    std::cout << complex << "\n";
  }

  template < typename Container >
  void print_container(Container& container)
  {
    for (auto l: container)
      std::cout << l << ", ";
    std::cout << "\n";
  }

  // /*
  //  *  \brief Verification if every simplex in the complex is witnessed by witnesses in knn.
  //  *  \param print_output =true will print the witnesses for each simplex
  //  *  \remark Added for debugging purposes.
  //  */
  // template< class KNearestNeighbors >
  // bool is_witness_complex(KNearestNeighbors const & knn, bool print_output) {
  //   for (Simplex_handle sh : sc_.complex_simplex_range()) {
  //     bool is_witnessed = false;
  //     typeVectorVertex simplex;
  //     int nbV = 0;  // number of verticed in the simplex
  //     for (Vertex_handle v : sc_.simplex_vertex_range(sh))
  //       simplex.push_back(v);
  //     nbV = simplex.size();
  //     for (typeVectorVertex w : knn) {
  //       bool has_vertices = true;
  //       for (Vertex_handle v : simplex)
  //         if (std::find(w.begin(), w.begin() + nbV, v) == w.begin() + nbV) {
  //           has_vertices = false;
  //         }
  //       if (has_vertices) {
  //         is_witnessed = true;
  //         if (print_output) {
  //           std::cout << "The simplex ";
  //           print_vector(simplex);
  //           std::cout << " is witnessed by the witness ";
  //           print_vector(w);
  //           std::cout << std::endl;
  //         }
  //         break;
  //       }
  //     }
  //     if (!is_witnessed) {
  //       if (print_output) {
  //         std::cout << "The following simplex is not witnessed ";
  //         print_vector(simplex);
  //         std::cout << std::endl;
  //       }
  //       assert(is_witnessed);
  //       return false;
  //     }
  //   }
  //   return true;
  // }
};

  /**
   *  \ingroup witness_complex
   *  \brief Iterative construction of the witness complex.
   *  \details The witness complex is written in simplicial complex sc_
   *   basing on a matrix knn of
   *  nearest neighbours of the form {witnesses}x{landmarks}.
   *
   *  The type KNearestNeighbors can be seen as
   *  Witness_range<Closest_landmark_range<Vertex_handle>>, where
   *  Witness_range and Closest_landmark_range are random access ranges.
   *
   *  Procedure takes into account at most (dim+1)
   *  first landmarks from each landmark range to construct simplices.
   *
   *  Landmarks are supposed to be in [0,nbL_-1]
   */
  // template <class KNearestNeighbors, class SimplicialComplexForWitness>
  // void witness_complex(KNearestNeighbors const & knn,
  //                      int nbL,
  //                      int dim,
  //                      SimplicialComplexForWitness & sc) {
  //   Witness_complex<SimplicialComplexForWitness>(knn, nbL, dim, sc);
  // }

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // WITNESS_COMPLEX_H_
