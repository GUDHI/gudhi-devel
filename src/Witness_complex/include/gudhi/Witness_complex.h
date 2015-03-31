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


  /** \addtogroup simplex_tree
   *  Witness complex is a simplicial complex defined on two sets of points in \f$\mathbf{R}^D\f$:
   *  \f$W\f$ set of witnesses and \f$L \subseteq W\f$ set of landmarks. The simplices are based on points in \f$L\f$
   *  and a simplex belongs to the witness complex if and only if it is witnessed (there exists a point \f$w \in W\f$ such that
   *  w is closer to the vertices of this simplex than others) and all of its faces are witnessed as well. 
   */
  template<typename FiltrationValue = double,
           typename SimplexKey      = int,
           typename VertexHandle    = int>
  class Witness_complex: public Simplex_tree<> {

  private:
    
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
  
    
    /** \brief Type for the vertex handle.
     *
     * Must be a signed integer type. It admits a total order <. */
    typedef VertexHandle Vertex_handle;
    
    /* Type of node in the simplex tree. */ 
    typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;
    /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
    typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;
    typedef typename Dictionary::iterator Simplex_handle;
  
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
        counter++;
        vv = {i};
        /* TODO Filtration */
        returnValue = insert_simplex(vv, Filtration_value(0.0));
        /* TODO Error if not inserted : normally no need here though*/
      }
      //std::cout << "Successfully added landmarks" << std::endl;
      // PRINT2
      //print_sc(root()); std::cout << std::endl;
      int u,v;     // two extremities of an edge
      if (nbL > 1) // if the supposed dimension of the complex is >0
        {
          for (int i=0; i != nbW; ++i)
            {
              // initial fill of active witnesses list
              u = knn[i][0];
              v = knn[i][1];
              vv = {u,v};
              returnValue = this->insert_simplex(vv,Filtration_value(0.0));
              //print_sc(root()); std::cout << std::endl;
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
      //std::cout << "Successfully added edges" << std::endl;
      while (!active_w.empty() && k < nbL )
        {
          //std::cout << "Started the step k=" << k << std::endl;
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
            }
          k++;
        }
      //print_sc(root()); std::cout << std::endl;
    }

    /** \brief Construction of witness complex from points given explicitly
     */
  void witness_complex_from_file(Point_Vector point_vector, int nbL)
  {
    std::vector<std::vector< int > > WL;
    furthestPoints(point_vector, point_vector.size(), nbL, WL);
    witness_complex(WL);
  }
  
private:

    /** \brief Print functions
    */
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
          if (has_children(map.begin()))
            print_sc(map.begin()->second.children());
          typename Dictionary::iterator it;
          for (it = map.begin()+1; it != map.end(); ++it)
            {
              std::cout << "," << it->first;
              if (has_children(it))
                print_sc(it->second.children());
            }
        }
      std::cout << ")";
    }
    
    template <typename KNearestNeighbours>
    bool all_faces_in(KNearestNeighbours &knn, int witness_id, int k, VertexHandle inserted_vertex)
    {
      //std::cout << "All face in with the landmark " << inserted_vertex << std::endl;
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
              //std::cout << "++++ finished loop safely\n";
            }//endif
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
    
    template <typename T>
    void print_vvector(std::vector< std::vector <T> > vv)
    {
      std::cout << "[";
      if (!vv.empty())
        {
          print_vector(*(vv.begin()));
          for (auto it = vv.begin()+1; it != vv.end(); ++it)
            {
              std::cout << ",";
              print_vector(*it);
            }
        }
      std::cout << "]\n";
    }
  
/**
 * \brief Landmark choice strategy by iteratively adding the landmark the furthest from the
 * current landmark set
 * \arg W is the vector of points which will be the witnesses
 * \arg nbP is the number of witnesses
 * \arg nbL is the number of landmarks
 * \arg WL is the matrix of the nearest landmarks with respect to witnesses (output)
 */

  template <typename KNearestNeighbours>
    void furthestPoints(Point_Vector &W, int nbP, int nbL, KNearestNeighbours &WL)
  //void furthestPoints(Point_Vector &W, int nbP, int nbL, Point_Vector &L)
  {
    //std::cout << "Enter furthestPoints "<< std::endl;
    // What is density and why we need it.
    // basically it is the stop indicator for "no more landmarks"
    //std::cout << "W="; print_vvector(W);
    //double density = 5.;
    std::vector< std::vector<double> > wit_land_dist(nbP,std::vector<double>()); // distance matrix witness x landmarks
    std::vector< int >                 chosen_landmarks;    // landmark list

    WL = KNearestNeighbours(nbP,std::vector<int>()); //nbP copies of empty vectors
    int current_number_of_landmarks=0;
    double curr_max_dist = 0;
    double curr_dist;
    double infty = std::numeric_limits<double>::infinity();
    std::vector< double > dist_to_L(nbP,infty);
    // double mindist = infty;
    int curr_max_w=0;
    int j;
    int temp_swap_int;
    double temp_swap_double;

    //CHOICE OF THE FIRST LANDMARK
    //std::cout << "Enter the first landmark stage\n";
    srand(354698);
    int rand_int = rand()% nbP;
    curr_max_w = rand_int;

    for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
      {
        chosen_landmarks.push_back(curr_max_w);
        curr_max_dist = 0;
        //std::cout << "**********Entered loop with current number of landmarks = " << current_number_of_landmarks << std::endl;
        //std::cout << "WL="; print_vvector(WL);
        //std::cout << "WLD="; print_vvector(wit_land_dist);
        //std::cout << "landmarks="; print_vector(chosen_landmarks); std::cout << std::endl;
        for (auto v: WL)
          v.push_back(current_number_of_landmarks);
        for (int i = 0; i < nbP; ++i)
          {
            //std::cout << "In the loop with i=" << i << " and landmark=" << chosen_landmarks[current_number_of_landmarks] << std::endl;
            //std::cout << "W[i]="; print_vector(W[i]); std::cout << " W[landmark]="; print_vector(W[chosen_landmarks[current_number_of_landmarks]]); std::cout << std::endl;
            if (i != chosen_landmarks[current_number_of_landmarks])
              curr_dist = euclidean_distance(W[i],W[chosen_landmarks[current_number_of_landmarks]]);
            else
              curr_dist = 0;
            //std::cout << "The problem is not in distance function\n";
            wit_land_dist[i].push_back(curr_dist);
            WL[i].push_back(current_number_of_landmarks);
            //std::cout << "Push't back\n";
            if (curr_dist < dist_to_L[i])
              dist_to_L[i] = curr_dist;
            if (dist_to_L[i] > curr_max_dist)
              {
                curr_max_dist = curr_dist;
                curr_max_w = i;
              }
            j = current_number_of_landmarks;
            //std::cout << "First half complete\n";
            while (j > 0 && wit_land_dist[i][j-1] > wit_land_dist[i][j])
              {
                temp_swap_int = WL[i][j];
                WL[i][j] = WL[i][j-1];
                WL[i][j-1] = temp_swap_int;
                temp_swap_double = wit_land_dist[i][j];
                wit_land_dist[i][j] = wit_land_dist[i][j-1];
                wit_land_dist[i][j-1] = temp_swap_double;
                --j;
              }
            //std::cout << "result WL="; print_vvector(WL);
            //std::cout << "result WLD="; print_vvector(wit_land_dist);
            //std::cout << "result distL="; print_vector(dist_to_L); std::cout << std::endl;
            //std::cout << "End loop\n";
          }
      }
    //std::cout << endl;
  }
  
}; //class Witness_complex

  
} // namespace Guhdi

#endif
