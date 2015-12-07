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


  /**  Witness complex is a simplicial complex defined on two sets of points in \f$\mathbf{R}^D\f$:
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
    
  private:
    int nbL;                   // Number of landmarks
    double density;            // Desired density

  public:

    /** \brief Set number of landmarks to nbL_
     */
    void setNbL(int nbL_)
    {
      nbL = nbL_;
    }

    /** \brief Set density to density_
     */
    void setDensity(double density_)
    {
      density = density_;
    }
    
    /**
     * /brief Iterative construction of the witness complex basing on a matrix of k nearest neighbours of the form {witnesses}x{landmarks}.
     * Landmarks are supposed to be in [0,nbL-1]
     */
    
    template< typename KNearestNeighbours >
    void witness_complex(KNearestNeighbours & knn)
    //void witness_complex(std::vector< std::vector< Vertex_handle > > & knn)
    {
      std::cout << "**Start the procedure witness_complex" << std::endl;
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
        returnValue = insert_simplex(vv, Filtration_value(0.0));
        /* TODO Error if not inserted : normally no need here though*/
      }
      int k=1; /* current dimension in iterative construction */
      for (int i=0; i != nbW; ++i)
        active_w.push_back(i);
      std::cout << "k=0, active witnesses: " << active_w.size() << std::endl;
      //std::cout << "Successfully added edges" << std::endl;
      int D = knn[0].size();
      while (!active_w.empty() && k < D )
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
                  returnValue = insert_simplex(simplex_vector,0.0);
                  it++;
                }
              else
                active_w.erase(it++); //First increase the iterator and then erase the previous element
            }
	  std::cout << "k=" << k << ", active witnesses: " << active_w.size() << std::endl;
          k++;
        }
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

  public:
    /** \brief Print functions
     */

    void st_to_file(std::ofstream& out_file)
    {
      sc_to_file(out_file, root());
    }

  private:
    void sc_to_file(std::ofstream& out_file, Siblings * sibl)
    {
      assert(sibl);
      children_to_file(out_file, sibl->members_);
    }
    
    void children_to_file(std::ofstream& out_file, Dictionary& map)
    {
      out_file << "(" << std::flush;
      if (!map.empty())
        {
          out_file << map.begin()->first << std::flush;
          if (has_children(map.begin()))
            sc_to_file(out_file, map.begin()->second.children());
          typename Dictionary::iterator it;
          for (it = map.begin()+1; it != map.end(); ++it)
            {
              out_file << "," << it->first << std::flush;
              if (has_children(it))
                sc_to_file(out_file, it->second.children());
            }
        }
      out_file << ")" << std::flush;
    }


    /** \brief Check if the facets of the k-dimensional simplex witnessed 
     *  by witness witness_id are already in the complex.
     *  inserted_vertex is the handle of the (k+1)-th vertex witnessed by witness_id
     */
    template <typename KNearestNeighbours>
    bool all_faces_in(KNearestNeighbours &knn, int witness_id, int k)
    {
      //std::cout << "All face in with the landmark " << inserted_vertex << std::endl;
      std::vector< VertexHandle > facet;
      //VertexHandle curr_vh = curr_sh->first;
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
          if (find(facet) == null_simplex())
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

  public:
/**
 * \brief Landmark choice strategy by iteratively adding the landmark the furthest from the
 * current landmark set
 * \arg W is the vector of points which will be the witnesses
 * \arg nbP is the number of witnesses
 * \arg nbL is the number of landmarks
 * \arg WL is the matrix of the nearest landmarks with respect to witnesses (output)
 */

  template <typename KNearestNeighbours>
    void landmark_choice_by_furthest_points(Point_Vector &W, int nbP, KNearestNeighbours &WL)
  {
    Point_Vector wit_land_dist(nbP,std::vector<double>());    // distance matrix witness x landmarks
    typeVectorVertex  chosen_landmarks;                       // landmark list

    WL = KNearestNeighbours(nbP,std::vector<int>());                             
    int current_number_of_landmarks=0;                        // counter for landmarks 
    double curr_max_dist = 0;                                 // used for defining the furhest point from L
    double curr_dist;                                         // used to stock the distance from the current point to L
    double infty = std::numeric_limits<double>::infinity();   // infinity (see next entry)
    std::vector< double > dist_to_L(nbP,infty);               // vector of current distances to L from points
    int curr_max_w=0;                                         // the point currently furthest from L 
    int j;
    int temp_swap_int;                                        
    double temp_swap_double;

    //CHOICE OF THE FIRST LANDMARK
    std::cout << "Enter the first landmark stage\n";
    srand(354698);
    int rand_int = rand()% nbP;
    curr_max_w = rand_int; //For testing purposes a pseudo-random number is used here

    for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
      {
        //curr_max_w at this point is the next landmark
        chosen_landmarks.push_back(curr_max_w);
        for (auto v: WL)
          v.push_back(current_number_of_landmarks);
        for (int i = 0; i < nbP; ++i)
          {
            curr_dist = euclidean_distance(W[i],W[chosen_landmarks[current_number_of_landmarks]]);
            wit_land_dist[i].push_back(curr_dist);
            WL[i].push_back(current_number_of_landmarks);
            if (curr_dist < dist_to_L[i])
              dist_to_L[i] = curr_dist;
            j = current_number_of_landmarks;
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
          }
        curr_max_dist = 0;
        for (int i = 0; i < nbP; ++i) {
          if (dist_to_L[i] > curr_max_dist)
            {
              curr_max_dist = dist_to_L[i];
              curr_max_w = i;
            }
        }
      }
  }

    /** \brief Landmark choice strategy by taking random vertices for landmarks.
     *
     */

    //    template <typename KNearestNeighbours>
    void landmark_choice_by_random_points(Point_Vector &W, int nbP, std::set<int> &L)
    {
      int current_number_of_landmarks=0;                        // counter for landmarks 

      srand(24660);
      int chosen_landmark = rand()%nbP;
      for (current_number_of_landmarks = 0; current_number_of_landmarks != nbL; current_number_of_landmarks++)
        {
          while (L.find(chosen_landmark) != L.end())
            {
              srand((int)clock());
              chosen_landmark = rand()% nbP;
            }
          L.insert(chosen_landmark);
        }
    }

    
    /** \brief Construct the matrix |W|x(D+1) of D+1 closest landmarks
     *  where W is the set of witnesses and D is the ambient dimension
     */
    template <typename KNearestNeighbours>
    void nearest_landmarks(Point_Vector &W, std::set<int> &L, KNearestNeighbours &WL)
    {
      int D = W[0].size();
      int nbP = W.size();
      WL = KNearestNeighbours(nbP,std::vector<int>());
      typedef std::pair<double,int> dist_i;
      typedef bool (*comp)(dist_i,dist_i);
      for (int W_i = 0; W_i < nbP; W_i++)
        {
          std::priority_queue<dist_i, std::vector<dist_i>, comp> l_heap([&](dist_i j1, dist_i j2){return j1.first > j2.first;});
          std::set<int>::iterator L_it;
          int L_i;
          for (L_it = L.begin(), L_i=0; L_it != L.end(); L_it++, L_i++)
            {
              dist_i dist = std::make_pair(euclidean_distance(W[W_i],W[*L_it]), L_i);
              l_heap.push(dist);
            }
          for (int i = 0; i < D+1; i++)
            {
              dist_i dist = l_heap.top();
              WL[W_i].push_back(dist.second);
              l_heap.pop();
            }
        }
    }


    
  public:
    /** \brief Verification if every simplex in the complex is witnessed
     */
    template< class KNearestNeighbors >
    bool is_witness_complex(KNearestNeighbors WL)
    {
      //bool final_result = true;
      for (Simplex_handle sh: complex_simplex_range())
        {
          bool is_witnessed = false;
          typeVectorVertex simplex;
          int nbV = 0; //number of verticed in the simplex
          for (int v: simplex_vertex_range(sh))
            simplex.push_back(v);
          nbV = simplex.size();
          for (typeVectorVertex w: WL)
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
                  std::cout << "The simplex ";
                  print_vector(simplex);
                  std::cout << " is witnessed by the witness ";
                  print_vector(w);
                  std::cout << std::endl;
                  break;
                }
            }
          if (!is_witnessed)
            {
              std::cout << "The following simplex is not witnessed ";
              print_vector(simplex);
              std::cout << std::endl;
              assert(is_witnessed);
              return false;
            }
        }
      return true; // Arrive here if the not_witnessed check failed all the time
    }

    
}; //class Witness_complex


  
} // namespace Guhdi

#endif
