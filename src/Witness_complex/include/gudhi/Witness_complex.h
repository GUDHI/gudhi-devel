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

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <gudhi/distance_functions.h>

#include <algorithm>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <ctime>
#include <iostream>

namespace Gudhi {

namespace witness_complex {

/** 
    \class Witness_complex
    \brief Constructs the witness complex for the given set of witnesses and landmarks.
    \ingroup witness_complex
 */
template< class Simplicial_complex>
class Witness_complex {
 private:
  struct Active_witness {
    int witness_id;
    int landmark_id;

    Active_witness(int witness_id_, int landmark_id_)
        : witness_id(witness_id_),
        landmark_id(landmark_id_) { }
  };

 private:
  typedef typename Simplicial_complex::Simplex_handle Simplex_handle;
  typedef typename Simplicial_complex::Vertex_handle Vertex_handle;

  typedef std::vector< double > Point_t;
  typedef std::vector< Point_t > Point_Vector;

  typedef std::vector< Vertex_handle > typeVectorVertex;
  typedef std::pair< typeVectorVertex, Filtration_value> typeSimplex;
  typedef std::pair< Simplex_handle, bool > typePairSimplexBool;

  typedef int Witness_id;
  typedef int Landmark_id;
  typedef std::list< Vertex_handle > ActiveWitnessList;

 private:
  int nbL;  // Number of landmarks
  Simplicial_complex& sc;  // Simplicial complex

 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** @name Constructor
   */

  //@{

  // Witness_range<Closest_landmark_range<Vertex_handle>>

  /**
   *  \brief Iterative construction of the witness complex.
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
  template< typename KNearestNeighbors >
  Witness_complex(KNearestNeighbors const & knn,
                  int nbL_,
                  int dim,
                  Simplicial_complex & sc_) : nbL(nbL_), sc(sc_) {
    // Construction of the active witness list
    int nbW = knn.size();
    typeVectorVertex vv;
    int counter = 0;
    /* The list of still useful witnesses
     * it will diminuish in the course of iterations
     */
    ActiveWitnessList active_w;  // = new ActiveWitnessList();
    for (Vertex_handle i = 0; i != nbL; ++i) {
      // initial fill of 0-dimensional simplices
      // by doing it we don't assume that landmarks are necessarily witnesses themselves anymore
      counter++;
      vv = {i};
      sc.insert_simplex(vv);
      // TODO(SK) Error if not inserted : normally no need here though
    }
    int k = 1; /* current dimension in iterative construction */
    for (int i = 0; i != nbW; ++i)
      active_w.push_back(i);
    while (!active_w.empty() && k < dim) {
      typename ActiveWitnessList::iterator it = active_w.begin();
      while (it != active_w.end()) {
        typeVectorVertex simplex_vector;
        /* THE INSERTION: Checking if all the subfaces are in the simplex tree*/
        bool ok = all_faces_in(knn, *it, k);
        if (ok) {
          for (int i = 0; i != k + 1; ++i)
            simplex_vector.push_back(knn[*it][i]);
          sc.insert_simplex(simplex_vector);
          // TODO(SK) Error if not inserted : normally no need here though
          it++;
        } else {
          active_w.erase(it++);  // First increase the iterator and then erase the previous element
        }
      }
      k++;
    }
  }

  //@}

 private:
  /** \brief Check if the facets of the k-dimensional simplex witnessed 
   *  by witness witness_id are already in the complex.
   *  inserted_vertex is the handle of the (k+1)-th vertex witnessed by witness_id
   */
  template <typename KNearestNeighbors>
  bool all_faces_in(KNearestNeighbors const &knn, int witness_id, int k) {
    std::vector< Vertex_handle > facet;
    // CHECK ALL THE FACETS
    for (int i = 0; i != k + 1; ++i) {
      facet = {};
      for (int j = 0; j != k + 1; ++j) {
        if (j != i) {
          facet.push_back(knn[witness_id][j]);
        }
      }  // endfor
      if (sc.find(facet) == sc.null_simplex())
        return false;
    }  // endfor
    return true;
  }

  template <typename T>
  static void print_vector(const std::vector<T>& v) {
    std::cout << "[";
    if (!v.empty()) {
      std::cout << *(v.begin());
      for (auto it = v.begin() + 1; it != v.end(); ++it) {
        std::cout << ",";
        std::cout << *it;
      }
    }
    std::cout << "]";
  }

 public:
  /**
   *  \brief Verification if every simplex in the complex is witnessed by witnesses in knn.
   *  \param print_output =true will print the witnesses for each simplex
   *  \remark Added for debugging purposes.
   */
  template< class KNearestNeighbors >
  bool is_witness_complex(KNearestNeighbors const & knn, bool print_output) {
    // bool final_result = true;
    for (Simplex_handle sh : sc.complex_simplex_range()) {
      bool is_witnessed = false;
      typeVectorVertex simplex;
      int nbV = 0;  // number of verticed in the simplex
      for (Vertex_handle v : sc.simplex_vertex_range(sh))
        simplex.push_back(v);
      nbV = simplex.size();
      for (typeVectorVertex w : knn) {
        bool has_vertices = true;
        for (Vertex_handle v : simplex)
          if (std::find(w.begin(), w.begin() + nbV, v) == w.begin() + nbV) {
            has_vertices = false;
            // break;
          }
        if (has_vertices) {
          is_witnessed = true;
          if (print_output) {
            std::cout << "The simplex ";
            print_vector(simplex);
            std::cout << " is witnessed by the witness ";
            print_vector(w);
            std::cout << std::endl;
          }
          break;
        }
      }
      if (!is_witnessed) {
        if (print_output) {
          std::cout << "The following simplex is not witnessed ";
          print_vector(simplex);
          std::cout << std::endl;
        }
        assert(is_witnessed);
        return false;
      }
    }
    return true;
  }
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // WITNESS_COMPLEX_H_
