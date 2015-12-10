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

#ifndef GUDHI_WITNESS_COMPLEX_H_
#define GUDHI_WITNESS_COMPLEX_H_

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/distance_functions.h"
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <math.h>
#include <ctime>
#include <iostream>

// Needed for the adjacency graph in bad link search
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace Gudhi {


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
          landmark_id(landmark_id_)
      {}
    };


  private:
  
    typedef typename Simplicial_complex::Simplex_handle Simplex_handle;
    typedef typename Simplicial_complex::Vertex_handle Vertex_handle;
      
    typedef std::vector< double > Point_t;
    typedef std::vector< Point_t > Point_Vector;

    typedef typename Simplicial_complex::Filtration_value Filtration_value;
    typedef std::vector< Vertex_handle > typeVectorVertex;
    typedef std::pair< typeVectorVertex, Filtration_value> typeSimplex;
    typedef std::pair< Simplex_handle, bool > typePairSimplexBool;
    
    typedef int Witness_id;
    typedef int Landmark_id;
    typedef std::list< Vertex_handle > ActiveWitnessList;
    
  private:
    int nbL;                   // Number of landmarks
    double density;            // Desired density
    Simplicial_complex& sc;     // Simplicial complex
    
  public:

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Constructor
     */

    //@{
    
    /**
     *  \brief Iterative construction of the witness complex.
     *  \details The witness complex is written in sc_ basing on a matrix knn of
     *  nearest neighbours of the form {witnesses}x{landmarks}.
     *  Parameter dim serves as the limit for the number of closest landmarks to consider.
     *  Landmarks are supposed to be in [0,nbL_-1]
     */    
    template< typename KNearestNeighbours >
    Witness_complex(KNearestNeighbours & knn,
                    Simplicial_complex & sc_,
                    int nbL_,
                    int dim ): nbL(nbL_), sc(sc_) 
    {
      //Construction of the active witness list
      int nbW = knn.size();
      typeVectorVertex vv;
      typeSimplex simplex;
      typePairSimplexBool returnValue;
      int counter = 0;
      /* The list of still useful witnesses
       * it will diminuish in the course of iterations
       */
      ActiveWitnessList active_w;// = new ActiveWitnessList();
      for (int i=0; i != nbL; ++i) {
        // initial fill of 0-dimensional simplices
        // by doing it we don't assume that landmarks are necessarily witnesses themselves anymore
        counter++;
        vv = {i};
        returnValue = sc.insert_simplex(vv, Filtration_value(0.0));
        /* TODO Error if not inserted : normally no need here though*/
      }
      int k=1; /* current dimension in iterative construction */
      for (int i=0; i != nbW; ++i)
        active_w.push_back(i);
      //std::cout << "Successfully added edges" << std::endl;
      while (!active_w.empty() && k < dim )
        {
          //std::cout << "Started the step k=" << k << std::endl;
          typename ActiveWitnessList::iterator it = active_w.begin();
          while (it != active_w.end())
            {
              typeVectorVertex simplex_vector;
              /* THE INSERTION: Checking if all the subfaces are in the simplex tree*/
              bool ok = all_faces_in(knn, *it, k);
              if (ok)
                {
                  for (int i = 0; i != k+1; ++i)
                    simplex_vector.push_back(knn[*it][i]);
                  returnValue = sc.insert_simplex(simplex_vector,0.0);
                  it++;
                }
              else
                active_w.erase(it++); //First increase the iterator and then erase the previous element
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
    template <typename KNearestNeighbours>
    bool all_faces_in(KNearestNeighbours &knn, int witness_id, int k)
    {
      //std::cout << "All face in with the landmark " << inserted_vertex << std::endl;
      std::vector< Vertex_handle > facet;
      //Vertex_handle curr_vh = curr_sh->first;
      // CHECK ALL THE FACETS
      for (int i = 0; i != k+1; ++i)
        {
          facet = {};
          for (int j = 0; j != k+1; ++j)
            {
              if (j != i)
                {
                  facet.push_back(knn[witness_id][j]);
                }
            }//endfor
          if (sc.find(facet) == sc.null_simplex())
            return false;
          //std::cout << "++++ finished loop safely\n";
        } //endfor
      return true;
    }

   template <typename T>
    void print_vector(std::vector<T> v)
    {
      std::cout << "[";
      if (!v.empty())
        {
          std::cout << *(v.begin());
          for (auto it = v.begin()+1; it != v.end(); ++it)
            {
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
    bool is_witness_complex(KNearestNeighbors & knn, bool print_output)
    {
      //bool final_result = true;
      for (Simplex_handle sh: sc.complex_simplex_range())
        {
          bool is_witnessed = false;
          typeVectorVertex simplex;
          int nbV = 0; //number of verticed in the simplex
          for (int v: sc.simplex_vertex_range(sh))
            simplex.push_back(v);
          nbV = simplex.size();
          for (typeVectorVertex w: knn)
            {
              bool has_vertices = true;  
              for (int v: simplex)
                if (std::find(w.begin(), w.begin()+nbV, v) == w.begin()+nbV)
                  {
                    has_vertices = false;
                    //break;
                  }
              if (has_vertices)
                {
                  is_witnessed = true;
                  if (print_output)
                    {
                      std::cout << "The simplex ";
                      print_vector(simplex);
                      std::cout << " is witnessed by the witness ";
                      print_vector(w);
                      std::cout << std::endl;
                    }
                  break;
                }
            }
          if (!is_witnessed)
            {
              if (print_output)
                {
                  std::cout << "The following simplex is not witnessed ";
                  print_vector(simplex);
                  std::cout << std::endl;
                }
              assert(is_witnessed);
              return false;
            }
        }
      return true; // Arrive here if the not_witnessed check failed all the time
    }

    
}; //class Witness_complex


  
} // namespace Guhdi

#endif
