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
#include "gudhi/Simplex_tree.h"
#include <vector>
#include <math.h>

namespace Gudhi {

  /*
template<typename FiltrationValue = double,
         typename SimplexKey      = int,
         typename VertexHandle    = int>
class Simplex_tree;
  */  

    /*  
template<typename FiltrationValue = double,
         typename SimplexKey      = int,
         typename VertexHandle    = int>
class Witness_complex: public Simplex_tree {
  //class Witness_complex: public Simplex_tree<FiltrationValue, SimplexKey, VertexHandle> {
public:

// typedef int                      Simplex_handle; //index in vector complex_
 
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
//  typedef VertexHandle Vertex_handle;

  /* Type of node in the simplex tree. */
//  typedef Simplex_tree_node_explicit_storage<Simplex_tree> Node;
  /* Type of dictionary Vertex_handle -> Node for traversing the simplex tree. */
//  typedef typename boost::container::flat_map<Vertex_handle, Node> Dictionary;

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

/*
  typedef std::vector< double > Point_t;
 typedef std::vector< Point_t > Point_Vector;

Witness_complex(int number_of_landmarks, std::string file_name)
{
  /*
  Point_Vector & points;
  read_points(file_name, points);
  landmark_extraction(points);
  */
/*  
}

private:


/**
 * \brief Permutes the vector in such a way that the landmarks appear first
 */

/*
void landmark_extraction(int nbP, Point_Vector & points, int inL)
{
  int i,j;
  for (i = 0; i != nbP; ++i)
    {
      for (j = 0; j != nbP; )
    }
}
*/

  /*
  
double distPoints(Point_t &p, Point_t &q) {
	double dist = 0.;
        Point_t::iterator itp = p.begin();
        Point_t::iterator itq = q.begin();
	while(itp!=p.end()) {
		dist += (*itp - *itq)*(*itp - *itq);
		itp++; itq++;}
	return dist;
}

void furthestPoints(Point_Vector &W, int nbP, std::string file_land, int dim, int nbL, Point_Vector &L) {
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
        curr_dist = distPoints(*itW,*itL);
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
  */
}; //class Witness_complex
*/
  
} // namespace Guhdi

#endif
