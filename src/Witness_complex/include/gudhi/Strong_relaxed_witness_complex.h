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

#ifndef STRONG_RELAXED_WITNESS_COMPLEX_H_
#define STRONG_RELAXED_WITNESS_COMPLEX_H_

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

namespace Gudhi {

namespace witness_complex {

  /** \addtogroup simplex_tree
   *  Witness complex is a simplicial complex defined on two sets of points in \f$\mathbf{R}^D\f$:
   *  \f$W\f$ set of witnesses and \f$L \subseteq W\f$ set of landmarks. The simplices are based on points in \f$L\f$
   *  and a simplex belongs to the witness complex if and only if it is witnessed (there exists a point \f$w \in W\f$ such that
   *  w is closer to the vertices of this simplex than others) and all of its faces are witnessed as well. 
   */
template< class Simplicial_complex >
class Strong_relaxed_witness_complex {

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
  typedef typename Simplicial_complex::Filtration_value FT;

  typedef std::vector< double > Point_t;
  typedef std::vector< Point_t > Point_Vector;

  // typedef typename Simplicial_complex::Filtration_value Filtration_value;
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
    /** @name Constructor
     */

    //@{
  
    /**
     *  \brief Iterative construction of the relaxed witness complex.
     *  \details The witness complex is written in sc_ basing on a matrix knn
     *  of k nearest neighbours of the form {witnesses}x{landmarks} and 
     *  and a matrix distances of distances to these landmarks from witnesses.
     *  The parameter alpha defines relaxation and
     *  limD defines the 
     *
     *  The line lengths in one given matrix can differ, 
     *  however both matrices have the same corresponding line lengths.
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
    
    template< typename KNearestNeighbours >
    Strong_relaxed_witness_complex(std::vector< std::vector<double> > const & distances,
                            KNearestNeighbours const & knn,
                            Simplicial_complex & sc_, 
                            int nbL_,
                            double alpha2,
                            unsigned limD) : nbL(nbL_), sc(sc_) {
      int nbW = knn.size();
      typeVectorVertex vv;
      //int counter = 0;
      /* The list of still useful witnesses
       * it will diminuish in the course of iterations
       */
      ActiveWitnessList active_w;// = new ActiveWitnessList();
      for (int i = 0; i != nbL; ++i) {
        // initial fill of 0-dimensional simplices
        // by doing it we don't assume that landmarks are necessarily witnesses themselves anymore
        //counter++;
        vv = {i};
        sc.insert_simplex(vv, Filtration_value(0.0));
        /* TODO Error if not inserted : normally no need here though*/
      }
      for (int i=0; i != nbW; ++i) {
        // int i_end = limD+1;
        // if (knn[i].size() < limD+1)
        //   i_end = knn[i].size();
        // double dist_wL = *(distances[i].begin());
        // while (distances[i][i_end] > dist_wL + alpha2)
        //   i_end--;
        // add_all_witnessed_faces(distances[i].begin(),
        //                         knn[i].begin(),
        //                         knn[i].begin() + i_end + 1);
        unsigned j_end = 0;
        while (j_end < distances[i].size() && j_end <= limD && distances[i][j_end] <= distances[i][0] + alpha2) {
          std::vector<int> simplex;
          for (unsigned j = 0; j <= j_end; ++j)
            simplex.push_back(knn[i][j]);
          assert(distances[i][j_end] - distances[i][0] >= 0);
          sc.insert_simplex_and_subfaces(simplex, distances[i][j_end] - distances[i][0]);
          j_end++; 
        }
      }
      sc.set_dimension(limD);
    }

  //@}
  
private:
  /* \brief Adds recursively all the faces of a certain dimension dim witnessed by the same witness
   * Iterator is needed to know until how far we can take landmarks to form simplexes
   * simplex is the prefix of the simplexes to insert
   * The output value indicates if the witness rests active or not
   */
  void add_all_witnessed_faces(std::vector<double>::const_iterator curr_d,
                               std::vector<int>::const_iterator curr_l,
                               std::vector<int>::const_iterator end)
  {
    std::vector<int> simplex;
    std::vector<int>::const_iterator l_end = curr_l;
    for (; l_end != end; ++l_end) {
      std::vector<int>::const_iterator l_it = curr_l;
      std::vector<double>::const_iterator d_it = curr_d;
      simplex = {};
      for (; l_it != l_end; ++l_it, ++d_it)
        simplex.push_back(*l_it);
      sc.insert_simplex_and_subfaces(simplex, *(d_it--) - *curr_d);
    }
  }
  
  /** \brief Check if the facets of the k-dimensional simplex witnessed 
   *  by witness witness_id are already in the complex.
   *  inserted_vertex is the handle of the (k+1)-th vertex witnessed by witness_id
   */
  bool all_faces_in(std::vector<int>& simplex, double* filtration_value)
  {
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

  bool is_face(Simplex_handle face, Simplex_handle coface)
  {
    // vertex range is sorted in decreasing order
    auto fvr = sc.simplex_vertex_range(face);
    auto cfvr = sc.simplex_vertex_range(coface);
    auto fv_it = fvr.begin();
    auto cfv_it = cfvr.begin();
    while (fv_it != fvr.end() && cfv_it != cfvr.end()) {
      if (*fv_it < *cfv_it)
        ++cfv_it;
      else if (*fv_it == *cfv_it) {
        ++cfv_it;
        ++fv_it;
      }
      else
        return false;
      
    }
    return (fv_it == fvr.end());
  }

  // void erase_simplex(Simplex_handle sh)
  // {
  //   auto siblings = sc.self_siblings(sh);
  //   auto oncles = siblings->oncles();
  //   int prev_vertex = siblings->parent();
  //   siblings->members().erase(sh->first);
  //   if (siblings->members().empty()) {
  //     typename typedef Simplicial_complex::Siblings Siblings;
  //     oncles->members().find(prev_vertex)->second.assign_children(new Siblings(oncles, prev_vertex));
  //     assert(!sc.has_children(oncles->members().find(prev_vertex)));
  //     //delete siblings;
  //   }

  // }
  
  void elementary_collapse(Simplex_handle face_sh, Simplex_handle coface_sh)
  {
    erase_simplex(coface_sh);
    erase_simplex(face_sh);
  }

public:
  // void collapse(std::vector<Simplex_handle>& simplices)
  // {
  //   // Get a vector of simplex handles ordered by filtration value
  //   std::cout << sc << std::endl;
  //   //std::vector<Simplex_handle> simplices;
  //   for (Simplex_handle sh: sc.filtration_simplex_range())
  //     simplices.push_back(sh);
  //   // std::sort(simplices.begin(),
  //   //           simplices.end(),
  //   //           [&](Simplex_handle sh1, Simplex_handle sh2)
  //   //           { double f1 = sc.filtration(sh1), f2 = sc.filtration(sh2);
  //   //             return (f1 > f2) || (f1 >= f2 && sc.dimension(sh1) > sc.dimension(sh2)); });
  //   // Double iteration
  //   auto face_it = simplices.rbegin();
  //   while (face_it != simplices.rend() && sc.filtration(*face_it) != 0) {
  //     int coface_count = 0;
  //     auto reduced_coface = simplices.rbegin();
  //     for (auto coface_it = simplices.rbegin(); coface_it != simplices.rend() && sc.filtration(*coface_it) != 0; ++coface_it)
  //       if (face_it != coface_it && is_face(*face_it, *coface_it)) {
  //         coface_count++;
  //         if (coface_count == 1)
  //           reduced_coface = coface_it;
  //         else
  //           break;
  //       }
  //     if (coface_count == 1) {
  //       std::cout << "Erase ( ";
  //       for (auto v: sc.simplex_vertex_range(*(--reduced_coface.base())))
  //         std::cout << v << " ";
        
  //       simplices.erase(--(reduced_coface.base()));
  //       //elementary_collapse(*face_it, *reduced_coface);
  //       std::cout << ") and then ( ";
  //       for (auto v: sc.simplex_vertex_range(*(--face_it.base())))
  //         std::cout << v << " ";
  //       std::cout << ")\n";
  //       simplices.erase(--((face_it++).base()));
  //       //face_it = simplices.rbegin();
  //       //std::cout << "Size of vector: " << simplices.size() << "\n";
  //     } 
  //     else
  //       face_it++;
  //   }
  //   sc.initialize_filtration();
  //   //std::cout << sc << std::endl;
  // }

  template <class Dim_lists>
  void collapse(Dim_lists& dim_lists)
  {
    dim_lists.collapse();
  }
  
private:
  
  /** Collapse recursively boundary faces of the given simplex
   *  if its filtration is bigger than alpha_lim.  
   */
  void rec_boundary_collapse(Simplex_handle sh, FT alpha_lim)
  {
    for (Simplex_handle face_it : sc.boundary_simplex_range()) {
      
    }
      
  }
  
}; //class Strong_relaxed_witness_complex

} // namespace witness_complex
  
} // namespace Gudhi

#endif
