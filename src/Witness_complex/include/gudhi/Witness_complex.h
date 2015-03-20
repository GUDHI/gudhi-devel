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

#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"
#include <vector>
#include <list>
#include <math.h>

namespace Gudhi {

  /*
template<typename FiltrationValue = double,
         typename SimplexKey      = int,
         typename VertexHandle    = int>
class Simplex_tree;
  */  

    /*
    typedef int Witness_id;
typedef int Landmark_id;
typedef int Simplex_handle; //index in vector complex_
    */

template<typename FiltrationValue = double,
         typename SimplexKey      = int,
         typename VertexHandle    = int>
class Witness_complex: public Simplex_tree<> {
    //class Witness_complex: public Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> {
    //class Witness_complex {
private:

  /*
  template < typename Witness_id = int,
             typename Landmark_id = int,
             typename Simplex_handle = Simplex_handle >
  struct Active_witness {
    Witness_id witness_id;
    Landmark_id landmark_id;
    Simplex_handle simplex_handle;
  };
  */
      
struct Active_witness {
    int witness_id;
    int landmark_id;
    Simplex_handle simplex_handle;

  Active_witness(int witness_id_, int landmark_id_, Simplex_handle simplex_handle_)
    : witness_id(witness_id_),
      landmark_id(landmark_id_),
      simplex_handle(simplex_handle_)
  {}
};
    
      


public:
  
//Simplex_tree<> st;
 
// typedef typename std::vector< Simplex_handle >::iterator    Boundary_simplex_iterator;
// typedef boost::iterator_range<Boundary_simplex_iterator>    Boundary_simplex_range;  
 
// typedef typename std::vector< Simplex_handle >::iterator    Skeleton_simplex_iterator;
// typedef boost::iterator_range< Skeleton_simplex_iterator >  Skeleton_simplex_range;

//  typedef IndexingTag Indexing_tag;
  /** \brief Type for the value of the filtration function.
   *
   * Must be comparable with <. */
//  typedef FiltrationValue Filtration_value;
  /** \brief Key associated to each simplex.
   *
   * Must be a signed integer type. */
//  typedef SimplexKey Simplex_key;

  /** \brief Type for the vertex handle.
   *
   * Must be a signed integer type. It admits a total order <. */
  typedef VertexHandle Vertex_handle;

  /* Type of node in the simplex tree. */ 
  typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;
  /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
  typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;
  typedef typename Dictionary::iterator Simplex_handle;
  
/*
  friend class Simplex_tree_node_explicit_storage< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_siblings< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle>, Dictionary>;
  friend class Simplex_tree_simplex_vertex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_boundary_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_complex_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
  friend class Simplex_tree_skeleton_simplex_iterator< Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> >;
*/

  /* \brief Set of nodes sharing a same parent in the simplex tree. */
  /* \brief Set of nodes sharing a same parent in the simplex tree. */
// typedef Simplex_tree_siblings<Simplex_tree, Dictionary> Siblings;


  typedef std::vector< double > Point_t;
  typedef std::vector< Point_t > Point_Vector;
  
  typedef std::vector< Vertex_handle > typeVectorVertex;
  typedef std::pair< typeVectorVertex, Filtration_value> typeSimplex;
  typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

  typedef int Witness_id;
  typedef int Landmark_id;
  typedef std::list< Active_witness > ActiveWitnessList;

/**
 * /brief Iterative construction of the witness complex basing on a matrix of k nearest neighbours of the form {witnesses}x{landmarks}.
 * Landmarks are supposed to be in [0,nbL-1]
 */

template< typename KNearestNeighbours >
void witness_complex(KNearestNeighbours & knn)
{
    int k=1; /* current dimension in iterative construction */
    //Construction of the active witness list
    int nbW = knn.size();
    int nbL = knn.at(0).size();
    //VertexHandle vh;
    typeVectorVertex vv;
    typeSimplex simplex;
    typePairSimplexBool returnValue;
    /* The list of still useful witnesses
     * it will diminuish in the course of iterations
     */
    ActiveWitnessList active_w;// = new ActiveWitnessList();
    for (int i=0; i != nbL; ++i) {
        // initial fill of 0-dimensional simplices
        // by doing it we don't assume that landmarks are necessarily witnesses themselves anymore
        //vh = (Vertex_handle)i;
        vv = {i};
        /* TODO Filtration */
        simplex = std::make_pair(vv, Filtration_value(0.0));
        returnValue = this->insert_simplex(simplex.first, simplex.second);
        /* TODO Error if not inserted : normally no need here though*/
    }
    int u,v;     // two extremities of an edge
    if (nbL > 1) // if the supposed dimension of the complex is >0
        for (int i=0; i != nbW; ++i) {
            // initial fill of active witnesses list
            u = knn[i][0];
            v = knn[i][1];
            if (u > v) {
                u = v;
                v = knn[i][0];
            }
            //Siblings * curr_sib = &root_;
            //vh = (Vertex_handle)i;
            vv = {u,v};
            returnValue = this->insert_simplex(vv,Filtration_value(0.0));
            active_w.push_back(Active_witness(i,v,returnValue.first));
            //Simplex_handle sh = root_.members_.begin()+u;
            //active_w.push_front(i);
        }
    while (!active_w.empty() && k+1 < nbL ) {
        typename ActiveWitnessList::iterator it = active_w.begin();
        while (it != active_w.end()) {
            typeVectorVertex simplex_vector;
            typeVectorVertex suffix;
            /*
            for (int i=0; i != k; ++i) {
                // create a (k+1) element array for the given active witness
                 simplex_vector.push_back(knn[it->witness_id][i]);
            }
            
            sort(simplex_vector.begin(),simplex_vector.end());
            */
            /* THE INSERTION: Checking if all the subfaces are in the simplex tree is mandatory */
            bool ok = all_faces_in(knn[it->witness_id][k],it->simplex_handle);
            if (ok) 
                returnValue = insert_simplex(simplex_vector,0.0);
            else
                active_w.erase(it++); //First increase the iterator and then erase the previous element
        }
    } 
}

private:

  bool all_faces_in(VertexHandle last, Simplex_handle sh)
  {
    std::list< VertexHandle > suffix;
    Simplex_handle curr_sh = sh;
    Siblings * curr_sibl = self_siblings(sh);
    VertexHandle curr_vh = curr_sh->first;
    while (curr_vh > last) {
      suffix.push_front(curr_vh);
      curr_vh = curr_sibl->parent();
      curr_sibl = curr_sibl->oncles();
      curr_sh = curr_sibl->find(curr_vh);
    }
    suffix.push_front(last);
    typename std::list< VertexHandle >::iterator itVV = suffix.begin();
    Simplex_handle sh_bup = curr_sh; // Back up the pointer
    while (itVV != suffix.end() && curr_sh->second.children()->find(*itVV) != null_simplex()) {
      // If the node doesn't exist then stop, else go down the tree
      curr_sh = curr_sh->second.children()->find(*itVV);
    }
    if (itVV == suffix.end()) {
      // the simplex is already in the tree
      return true;
    }
    else if (itVV != suffix.end()) {
      // the father of the simplex is not in the tree
      return false;
    }
    else {
      // CHECK ALL THE FACETS
      curr_sh = sh_bup;
      while (curr_sibl != root()) {
        suffix.push_front(curr_vh);
        curr_vh = curr_sibl->parent();
        curr_sibl = curr_sibl->oncles();
        curr_sh = curr_sibl->find(curr_vh);
      }
      suffix.push_front(curr_vh);
      sh_bup = curr_sh; // the first vertex lexicographicly
      for (typename std::list< VertexHandle >::iterator itExcl = suffix.begin(); itExcl != suffix.end(); ++itExcl) {
        if (*itExcl != last) {
          itVV = suffix.begin();
          while (itVV != itExcl) {
            if (curr_sibl->find(*itVV) == null_simplex())
              return false;
            curr_sh = curr_sibl->find(*itVV);
            curr_sibl = self_siblings(curr_sh);
            itVV++;
          }
          itVV++;
          while (itVV != suffix.end()) {
            if (curr_sibl->find(*itVV) == null_simplex())
              return false;
            curr_sh = curr_sibl->find(*itVV);
            curr_sibl = self_siblings(curr_sh);
            itVV++;
          } //endwhile
        } //endif
      } //endfor
      return true;
    } //end check all the facets
  }

/**
 * \brief Permutes the vector in such a way that the landmarks appear first
 */

  void furthestPoints(Point_Vector &W, int nbP, std::string file_land, int dim, int nbL, Point_Vector &L)
  {
  //std::cout << "Enter furthestPoints "<< endl;
  //Point_Vector *L = new Point_Vector();
    double density = 5.;
    int current_number_of_landmarks=0;
    double curr_max_dist;
    double curr_dist;
    double mindist = 10005.;
    int curr_max_w=0;
    int curr_w=0;
    srand(354698);
    int rand_int = rand()% nbP;
    //std::cout << rand_int << endl;
    L.push_back(W[rand_int]);// first landmark is random
    current_number_of_landmarks++;
    while (1) {
      curr_w = 0;
      curr_max_dist = -1;
      for(Point_Vector::iterator itW = W.begin(); itW != W.end(); itW++) {
        //compute distance from w and L
        mindist = 100000.;
        for(Point_Vector::iterator itL = L.begin(); itL != L.end(); itL++) {
          //curr_dist = distPoints(*itW,*itL);
          curr_dist = euclidean_distance(*itW,*itL);
          if(curr_dist < mindist) {
            mindist = curr_dist;
          }
        }
        if(mindist > curr_max_dist) {
          curr_max_w = curr_w; //???
          curr_max_dist = mindist;
        }
        curr_w++;
      }
      L.push_back(W[curr_max_w]);
      current_number_of_landmarks++;
      density = sqrt(curr_max_dist);
      //std::cout << "[" << current_number_of_landmarks << ":" << density <<"] ";
      if(L.size() == nbL) break;
    }
    //std::cout << endl;
    return L;
  }
  
}; //class Witness_complex

  
} // namespace Guhdi

#endif
