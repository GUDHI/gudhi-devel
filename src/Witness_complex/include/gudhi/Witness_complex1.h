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
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"
#include <vector>
#include <list>
#include <limits>
#include <math.h>

#include <iostream>

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
  typedef std::list< Vertex_handle > ActiveWitnessList;

/**
 * /brief Iterative construction of the witness complex basing on a matrix of k nearest neighbours of the form {witnesses}x{landmarks}.
 * Landmarks are supposed to be in [0,nbL-1]
 */

  
template< typename KNearestNeighbours >
void witness_complex(KNearestNeighbours & knn)
//void witness_complex(std::vector< std::vector< Vertex_handle > > & knn)
{
  std::cout << "**Start the procedure witness_complex" << std::endl;
    int k=2; /* current dimension in iterative construction */
    //Construction of the active witness list
    int nbW = knn.size();
    int nbL = knn.at(0).size();
    //VertexHandle vh;
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
        //vh = (Vertex_handle)i;
      counter++;
        vv = {i};
        /* TODO Filtration */
        //simplex = std::make_pair(vv, Filtration_value(0.0));
        //returnValue = this->insert_simplex(simplex.first, simplex.second);
        returnValue = insert_simplex(vv, Filtration_value(0.0));
        /* TODO Error if not inserted : normally no need here though*/
    }
    //vv = {0};
    //returnValue = insert_simplex(vv,Filtration_value(0.0));
    std::cout << "Successfully added landmarks" << std::endl;
    // PRINT2
    print_sc(root()); std::cout << std::endl;
    int u,v;     // two extremities of an edge
    if (nbL > 1) // if the supposed dimension of the complex is >0
      {
        for (int i=0; i != nbW; ++i)
          {
            // initial fill of active witnesses list
            u = knn[i][0];
            v = knn[i][1];
            //Siblings * curr_sib = &root_;
            //vh = (Vertex_handle)i;
            vv = {u,v};
            returnValue = this->insert_simplex(vv,Filtration_value(0.0));
            print_sc(root()); std::cout << std::endl;
            //std::cout << "Added edges" << std::endl;
          }
        //print_sc(root());
        for (int i=0; i != nbW; ++i)
          {
            // initial fill of active witnesses list
            u = knn[i][0];
            v = knn[i][1];
            if ( u > v)
              {
                u = v;
                v = knn[i][0];
                knn[i][0] = knn[i][1];
                knn[i][1] = v;
              }
            Simplex_handle sh;
            vv = {u,v};
            sh = (root()->find(u))->second.children()->find(v);
            active_w.push_back(i);
        }
      }
    std::cout << "Successfully added edges" << std::endl;
    while (!active_w.empty() && k+1 < nbL )
      {
        std::cout << "Started the step k=" << k << std::endl;
        typename ActiveWitnessList::iterator it = active_w.begin();
        while (it != active_w.end())
          {
            typeVectorVertex simplex_vector;
            /* THE INSERTION: Checking if all the subfaces are in the simplex tree*/
            // First sort the first k landmarks
            VertexHandle inserted_vertex = knn[*it][k];
            bool ok = all_faces_in(knn, *it, k, inserted_vertex);
            if (ok)
              {
                for (int i = 0; i != k+1; ++i)
                  simplex_vector.push_back(knn[*it][i]);
                returnValue = insert_simplex(simplex_vector,0.0);
                it++;
              }
            else
                active_w.erase(it++); //First increase the iterator and then erase the previous element
            print_sc(root()); std::cout << std::endl;
          }
        k++;
    } 
}

  void witness_complex_from_file(std::string file_name, int nbL)
  {
    //READ THE FILE INTO A POINT VECTOR
    std::vector< std::vector< double > > point_vector;
    read_points(file_name, point_vector);
    std::vector<std::vector< int > > WL;
    furthestPoints(point_vector, point_vector.size(), nbL, WL);
    witness_complex(WL);
  }
  
private:

  void print_sc(Siblings * sibl)
  {
    if (sibl == NULL)
      std::cout << "&";
    else
      print_children(sibl->members_);
  }

  void print_children(Dictionary map)
  {
    std::cout << "(";
    if (!map.empty())
      {
        std::cout << map.begin()->first;
        /*if (map.begin()->second.children() == root())
          std::cout << "Sweet potato"; */
        if (has_children(map.begin()))
          print_sc(map.begin()->second.children());
        typename Dictionary::iterator it;
        for (it = map.begin()+1; it != map.end(); ++it)
          {
            std::cout << "," << it->first;
            /*if (map.begin()->second.children() == root())
              std::cout << "Sweet potato";*/
            if (has_children(it))
              print_sc(it->second.children());
          }
      }
    std::cout << ")";
  }

  template <typename KNearestNeighbours>
  bool all_faces_in(KNearestNeighbours &knn, int witness_id, int k, VertexHandle inserted_vertex)
  {
    std::cout << "All face in with the landmark " << inserted_vertex << std::endl;
    std::vector< VertexHandle > facet;
    //VertexHandle curr_vh = curr_sh->first;
    // CHECK ALL THE FACETS
    for (int i = 0; i != k+1; ++i)
      {
        if (knn[witness_id][i] != inserted_vertex)
          {
            facet = {};
            for (int j = 0; j != k+1; ++j)
              {
                if (j != i)
                  {
                    facet.push_back(knn[witness_id][j]);
                  }
              }//endfor
            if (find(facet) == null_simplex())
              return false;
            std::cout << "++++ finished loop safely\n";
          }//endif
      } //endfor
      return true;
  }

/**
 * \brief Permutes the vector in such a way that the landmarks appear first
 * \arg W is the vector of points which will be the witnesses
 * \arg nbP is the number of witnesses
 * \arg nbL is the number of landmarks
 * \arg WL is the matrix of the nearest landmarks with respect to witnesses (output)
 */

  template <typename KNearestNeighbours>
    void furthestPoints(Point_Vector &W, int nbP, int nbL, KNearestNeighbours &WL)
  //void furthestPoints(Point_Vector &W, int nbP, int nbL, Point_Vector &L)
  {
    std::cout << "Enter furthestPoints "<< std::endl;
    // What is density and why we need it.
    // basically it is the stop indicator for "no more landmarks"
  
    //double density = 5.;
    std::vector< std::vector<double> > wit_land_dist(nbP,std::vector<double>()); // distance matrix witness x landmarks
    std::vector< int >                 chosen_landmarks;    // landmark list
    WL = KNearestNeighbours(nbP,std::vector<int>()); //nbP copies of empty vectors
    int current_number_of_landmarks=0;
    double curr_max_dist = 0;
    double curr_dist;
    //double infty = std::numeric_limits<double>::infinity();
    // double mindist = infty;
    int curr_max_w=0;
    int j;
    int temp_swap_int;
    double temp_swap_double;

    //CHOICE OF THE FIRST LANDMARK
    std::cout << "Enter the first landmark stage\n";
    srand(354698);
    int rand_int = rand()% nbP;
    curr_max_w = rand_int;

    for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
      {
        chosen_landmarks.push_back(curr_max_w); // first landmark is random
        for (auto v: WL)
          v.push_back(current_number_of_landmarks);
        for (unsigned int i = 0; i < wit_land_dist.size(); ++i)
          {
            curr_dist = euclidean_distance(W[i],W[chosen_landmarks[current_number_of_landmarks]]);
            wit_land_dist[i].push_back(curr_dist);
            if (curr_dist > curr_max_dist)
              {
                curr_max_dist = curr_dist;
                curr_max_w = i;
              }
            j = current_number_of_landmarks;
            while (j > 0 && wit_land_dist[i][j-1] > wit_land_dist[i][j])
              {
                temp_swap_int = WL[i][j];
                WL[i][j] = WL[i][j-1];
                WL[i][j-1] = temp_swap_int;
                temp_swap_double = wit_land_dist[i][j];
                wit_land_dist[i][j] = wit_land_dist[i][j-1];
                wit_land_dist[i][j-1] = temp_swap_double;
              }
          }
      }
    /*
    while (1) {
      curr_w = 0;
      curr_max_dist = -1;
      for(Point_Vector::iterator itW = W.begin(); itW != W.end(); itW++) {
        //compute distance from w and L
        mindist = infty;
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
      //density = sqrt(curr_max_dist);
      //std::cout << "[" << current_number_of_landmarks << ":" << density <<"] ";
      if(L.size() == nbL) break;
    }
    */
    //std::cout << endl;
    //return L;
  }
  
}; //class Witness_complex

  
} // namespace Guhdi

#endif
