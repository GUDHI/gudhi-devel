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

#ifndef GOOD_LINKS_H_
#define GOOD_LINKS_H_

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>
#include <utility>
#include "gudhi/reader_utils.h"
#include "gudhi/distance_functions.h"
#include <gudhi/Dim_list_iterator.h>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <limits>
#include <math.h>
#include <ctime>
#include <iostream>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

// Needed for the adjacency graph in bad link search
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace Gudhi {
  
template < typename Simplicial_complex >
class Good_links {
  
  typedef typename Simplicial_complex::Simplex_handle Simplex_handle;
  typedef typename Simplicial_complex::Vertex_handle Vertex_handle;
  typedef std::vector<Vertex_handle> Vertex_vector;
  
  // graph is an adjacency list
  typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Adj_graph;
  // map that gives to a certain simplex its node in graph and its dimension
  //typedef std::pair<boost::vecS,Color> Reference;
  typedef boost::graph_traits<Adj_graph>::vertex_descriptor Vertex_t;
  typedef boost::graph_traits<Adj_graph>::edge_descriptor Edge_t;
  typedef boost::graph_traits<Adj_graph>::adjacency_iterator Adj_it;
  typedef std::pair<Adj_it, Adj_it> Out_edge_it;
  
  typedef boost::container::flat_map<Simplex_handle, Vertex_t> Graph_map;
  typedef boost::container::flat_map<Vertex_t, Simplex_handle> Inv_graph_map;

public:
  Good_links(Simplicial_complex& sc): sc_(sc)
  {
    int dim = 0;
    for (auto sh: sc_.complex_simplex_range())
      if (sc_.dimension(sh) > dim)
        dim = sc_.dimension(sh);
    dimension = dim;
    count_good = Vertex_vector(dim);
    count_bad = Vertex_vector(dim);
  }
  
private:
  
  Simplicial_complex& sc_;
  unsigned dimension;
  std::vector<int> count_good;
  std::vector<int> count_bad;
  
  void add_vertices_to_link_graph(Vertex_vector& star_vertices,
                                  typename Vertex_vector::iterator curr_v,
                                  Adj_graph& adj_graph,
                                  Graph_map& d_map,
                                  Graph_map& f_map,
                                  Vertex_vector& curr_simplex,
                                  int curr_d,
                                  int link_dimension)
  {
    Simplex_handle sh;
    Vertex_t vert;
    typename Vertex_vector::iterator it;
    //std::pair<typename Graph_map::iterator,bool> resPair;
    //typename Graph_map::iterator resPair;
    //Add vertices
    //std::cout << "Entered add vertices\n";
    for (it = curr_v; it != star_vertices.end(); ++it) {
      curr_simplex.push_back(*it);               //push next vertex in question
      curr_simplex.push_back(star_vertices[0]);  //push the center of the star
      /*
        std::cout << "Searching for ";
        print_vector(curr_simplex);
        std::cout << " curr_dim " << curr_d << " d " << dimension << "";
      */
      Vertex_vector curr_simplex_copy(curr_simplex);
      sh = sc_.find(curr_simplex_copy);                   //a simplex of the star
      curr_simplex.pop_back();                   //pop the center of the star
      curr_simplex_copy = Vertex_vector(curr_simplex);
      if (sh != sc_.null_simplex()) {
        //std::cout << " added\n";
        if (curr_d == link_dimension) {
          sh = sc_.find(curr_simplex_copy);               //a simplex of the link
          assert(sh != sc_.null_simplex());          //ASSERT!
          vert = boost::add_vertex(adj_graph);
          d_map.emplace(sh,vert);
        }
        else {        
          if (curr_d == link_dimension-1) {
            sh = sc_.find(curr_simplex_copy);               //a simplex of the link
            assert(sh != sc_.null_simplex());
            vert = boost::add_vertex(adj_graph);
            f_map.emplace(sh,vert);
          }
          
          //delete (&curr_simplex_copy); //Just so it doesn't stack
          add_vertices_to_link_graph(star_vertices,
                                     it+1,
                                     adj_graph,
                                     d_map,
                                     f_map,
                                     curr_simplex,
                                     curr_d+1, link_dimension);
        }
      }
      /*
        else
        std::cout << "\n";
      */
      curr_simplex.pop_back();                           //pop the vertex in question
    }
  }

  void add_edges_to_link_graph(Adj_graph& adj_graph,
                               Graph_map& d_map,
                               Graph_map& f_map)
    {
      Simplex_handle sh;
      // Add edges
      //std::cout << "Entered add edges:\n";
      typename Graph_map::iterator map_it;
      for (auto d_map_pair : d_map) {
        //std::cout << "*";
        sh = d_map_pair.first;
        Vertex_t d_vert = d_map_pair.second;
        for (auto facet_sh : sc_.boundary_simplex_range(sh))
          //for (auto f_map_it : f_map)
          {
            //std::cout << "'"; 
            map_it = f_map.find(facet_sh);
            //We must have all the facets in the graph at this point
            assert(map_it != f_map.end());
            Vertex_t f_vert = map_it->second;
            //std::cout << "Added edge " << sh->first << "-" << map_it->first->first << "\n";
            boost::add_edge(d_vert,f_vert,adj_graph);
          }
      }
    }


  
  /* \brief Verifies if the simplices formed by vertices given by link_vertices 
   * form a pseudomanifold.
   * The idea is to make a bipartite graph, where vertices are the d- and (d-1)-dimensional
   * faces and edges represent adjacency between them.
   */
  bool link_is_pseudomanifold(Vertex_vector& star_vertices,
                              int dimension)
  {
    Adj_graph adj_graph;
    Graph_map d_map, f_map;
    // d_map = map for d-dimensional simplices
    // f_map = map for its facets
    Vertex_vector init_vector = {};
    add_vertices_to_link_graph(star_vertices,
                               star_vertices.begin()+1,
                               adj_graph,
                               d_map,
                               f_map,
                               init_vector,
                               0, dimension);
    //std::cout << "DMAP_SIZE: " << d_map.size() << "\n";
    //std::cout << "FMAP_SIZE: " << f_map.size() << "\n";
    add_edges_to_link_graph(adj_graph,
                            d_map,
                            f_map);
    for (auto f_map_it : f_map) {
      //std::cout << "Degree of " << f_map_it.first->first << " is " << boost::out_degree(f_map_it.second, adj_graph) << "\n";
      if (boost::out_degree(f_map_it.second, adj_graph) != 2){
        count_bad[dimension]++;
        return false;
      }
    }
    // At this point I know that all (d-1)-simplices are adjacent to exactly 2 d-simplices
    // What is left is to check the connexity
    //std::vector<int> components(boost::num_vertices(adj_graph));
    return true; //Forget the connexity
    //return (boost::connected_components(adj_graph, &components[0]) == 1);
  }

  int star_dim(Vertex_vector& star_vertices,
               typename Vertex_vector::iterator curr_v,
               int curr_d,
               Vertex_vector& curr_simplex,
               typename std::vector< int >::iterator curr_dc)
  {
    //std::cout << "Entered star_dim for " << *(curr_v-1) << "\n";
    Simplex_handle sh;
    int final_d = curr_d;
      typename Vertex_vector::iterator it;
      typename Vertex_vector::iterator dc_it;
      //std::cout << "Current vertex is " <<  
      for (it = curr_v, dc_it = curr_dc; it != star_vertices.end(); ++it, ++dc_it)
        {
          curr_simplex.push_back(*it);
          Vertex_vector curr_simplex_copy(curr_simplex);
          /*
          std::cout << "Searching for ";
          print_vector(curr_simplex);
          std::cout << " curr_dim " << curr_d << " final_dim " << final_d;
          */
          sh = sc_.find(curr_simplex_copy); //Need a copy because find sorts the vector and I want star center to be the first
          if (sh != sc_.null_simplex())
            {
              //std::cout << " -> " << *it << "\n";
              int d = star_dim(star_vertices,
                               it+1,
                               curr_d+1,
                               curr_simplex,
                               dc_it);
              if (d >= final_d)
                {
                  final_d = d;
                  //std::cout << d << " ";
                  //print_vector(curr_simplex);
		  //std::cout << std::endl;
                }
              if (d >= *dc_it)
                *dc_it = d;
            }
          /*
          else
            std::cout << "\n";
          */
          curr_simplex.pop_back();
        }
      return final_d;
    }

public:
  
  /** \brief Returns true if the link is good
   */
  bool has_good_link(Vertex_handle v)
  {
    typedef Vertex_vector typeVectorVertex;
    typeVectorVertex star_vertices;
    // Fill star_vertices
    star_vertices.push_back(v);
    for (auto u: sc_.complex_vertex_range())
      {
        typeVectorVertex edge = {u,v};
        if (u != v && sc_.find(edge) != sc_.null_simplex())   
          star_vertices.push_back(u);
      }
    // Find the dimension
    typeVectorVertex init_simplex = {star_vertices[0]};
    bool is_pure = true;
    std::vector<int> dim_coface(star_vertices.size(), 1);
    int d = star_dim(star_vertices,
                     star_vertices.begin()+1,
                     0,
                     init_simplex,
                     dim_coface.begin()+1) - 1; //link_dim = star_dim - 1
    
    assert(init_simplex.size() == 1);
    // if (!is_pure)
    //   std::cout << "Found an impure star around " << v << "\n";
    for (int dc: dim_coface)
      is_pure = (dc == dim_coface[0]);
    /*
      if (d == count_good.size())
      {
      std::cout << "Found a star of dimension " << (d+1) << " around " << v << "\nThe star is ";
      print_vector(star_vertices); std::cout << std::endl;
      }
    */
    //if (d == -1) count_bad[0]++;
    bool b= (is_pure && link_is_pseudomanifold(star_vertices,d));
    if (d != -1) {if (b) count_good[d]++; else count_bad[d]++;}
    if (!is_pure) count_bad[0]++;
    return (d != -1 && b && is_pure);
  }

private:
  void add_max_simplices_to_graph(Vertex_vector& star_vertices,
                                  typename Vertex_vector::iterator curr_v,
                                  Adj_graph& adj_graph,
                                  Graph_map& d_map,
                                  Graph_map& f_map,
                                  Inv_graph_map& inv_d_map,
                                  Vertex_vector& curr_simplex,
                                  int curr_d,
                                  int link_dimension)
  {
    Simplex_handle sh;
    Vertex_t vert;
    typename Vertex_vector::iterator it;
    //std::pair<typename Graph_map::iterator,bool> resPair;
    //typename Graph_map::iterator resPair;
    //Add vertices
    //std::cout << "Entered add vertices\n";
    for (it = curr_v; it != star_vertices.end(); ++it) {
      curr_simplex.push_back(*it);               //push next vertex in question
      //curr_simplex.push_back(star_vertices[0]);  //push the center of the star
      /*
        std::cout << "Searching for ";
        print_vector(curr_simplex);
        std::cout << " curr_dim " << curr_d << " d " << dimension << "";
      */
      Vertex_vector curr_simplex_copy(curr_simplex);
      sh = sc_.find(curr_simplex_copy);                   //a simplex of the star
      //curr_simplex.pop_back();                   //pop the center of the star
      curr_simplex_copy = Vertex_vector(curr_simplex);
      if (sh != sc_.null_simplex()) {
        //std::cout << " added\n";
        if (curr_d == link_dimension) {
          sh = sc_.find(curr_simplex_copy);               //a simplex of the link
          assert(sh != sc_.null_simplex()); //ASSERT!
          vert = boost::add_vertex(adj_graph);
          d_map.emplace(sh,vert);
          inv_d_map.emplace(vert,sh);
        }
        else {      
          if (curr_d == link_dimension-1) {
              sh = sc_.find(curr_simplex_copy);               //a simplex of the link
              assert(sh != sc_.null_simplex());
              vert = boost::add_vertex(adj_graph);
              f_map.emplace(sh,vert);
          }        
          //delete (&curr_simplex_copy); //Just so it doesn't stack
          add_max_simplices_to_graph(star_vertices,
                                     it+1,
                                     adj_graph,
                                     d_map,
                                     f_map,
                                     inv_d_map,
                                     curr_simplex,
                                     curr_d+1,
                                     link_dimension);
        }
      }
      /*
        else
        std::cout << "\n";
      */
      curr_simplex.pop_back();                           //pop the vertex in question
    }
  }

public:
  bool complex_is_pseudomanifold()
  {
    Adj_graph adj_graph;
    Graph_map d_map, f_map;
    // d_map = map for d-dimensional simplices
    // f_map = map for its facets
    Inv_graph_map inv_d_map;
    Vertex_vector init_vector = {};
    std::vector<int> star_vertices;
    for (int v: sc_.complex_vertex_range())
      star_vertices.push_back(v);
    add_max_simplices_to_graph(star_vertices,
                               star_vertices.begin(),
                               adj_graph,
                               d_map,
                               f_map,
                               inv_d_map,
                               init_vector,
                               0, dimension);
    //std::cout << "DMAP_SIZE: " << d_map.size() << "\n";
    //std::cout << "FMAP_SIZE: " << f_map.size() << "\n";
    add_edges_to_link_graph(adj_graph,
                            d_map,
                            f_map);
    for (auto f_map_it : f_map) {
      //std::cout << "Degree of " << f_map_it.first->first << " is " << boost::out_degree(f_map_it.second, adj_graph) << "\n";
      if (boost::out_degree(f_map_it.second, adj_graph) != 2) {
        // if (boost::out_degree(f_map_it.second, adj_graph) >= 3) {		 
        //   // std::cout << "This simplex has 3+ cofaces: ";
        //   // for(auto v : sc_.simplex_vertex_range(f_map_it.first))
        //   //   std::cout << v << " ";
        //   // std::cout << std::endl;
        //   Adj_it ai, ai_end; 
        //   for (std::tie(ai, ai_end) = boost::adjacent_vertices(f_map_it.second, adj_graph); ai != ai_end; ++ai) {
        //     auto it = inv_d_map.find(*ai);
        //     assert (it != inv_d_map.end());
        //     typename Simplicial_complex::Simplex_handle sh = it->second;
        //     for(auto v : sc_.simplex_vertex_range(sh))
        //       std::cout << v << " ";
        //     std::cout << std::endl;
        //   }
        // }
        count_bad[dimension]++;
        return false;
      }
    }
    // At this point I know that all (d-1)-simplices are adjacent to exactly 2 d-simplices
    // What is left is to check the connexity
    //std::vector<int> components(boost::num_vertices(adj_graph));
    return true; //Forget the connexity
    //return (boost::connected_components(adj_graph, &components[0]) == 1);
  }

  int number_good_links(int dim)
  {
    return count_good[dim];
  }
  
  int number_bad_links(int dim)
  {
    return count_bad[dim];
  }
  
};
    
} // namespace Gudhi

#endif
