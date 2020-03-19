/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam 
 *
 *    Copyright (C) 2018 INRIA Sophia Antipolis (France)
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
#pragma once

#include <gudhi/Rips_edge_list.h>
#include <boost/functional/hash.hpp>
// #include <boost/graph/adjacency_list.hpp>


#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_map>
#include <tuple>
#include <list>
#include <algorithm>
#include <chrono>

#include <ctime>
#include <fstream>

#include <Eigen/Sparse>


typedef std::size_t Vertex;
using Edge 					= std::pair<Vertex,Vertex>;  // This is an ordered pair, An edge is stored with convention of the first element being the smaller i.e {2,3} not {3,2}. However this is at the level of row indices on actual vertex lables
using EdgeFilt              = std::pair<Edge, double>;
using edge_list 			= std::vector<Edge>;

using MapVertexToIndex 	  	= std::unordered_map<Vertex, std::size_t>;
using Map 				  	= std::unordered_map<Vertex,Vertex>;

using sparseRowMatrix  	 	= Eigen::SparseMatrix<double, Eigen::RowMajor> ;
using rowInnerIterator 		= sparseRowMatrix::InnerIterator;

using intVector 	   = std::vector<int>;
using doubleVector 	   = std::vector<double>;
using vertexVector     = std::vector<Vertex>;
using boolVector       = std::vector<bool>;

using doubleQueue 	   = std::queue<double>;
using edgeQueue		   = std::queue<Edge>;

using EdgeFiltQueue    = std::queue<EdgeFilt>;
using EdgeFiltVector   = std::vector<EdgeFilt>;

typedef std::vector< std::tuple< double, Vertex, Vertex > > Filtered_sorted_edge_list;
typedef std::unordered_map<Edge, bool, boost::hash< Edge > > u_edge_map;
typedef std::unordered_map<Edge, std::size_t, boost::hash< Edge > > u_edge_to_idx_map;


//!  Class SparseMsMatrix 
/*!
  The class for storing the Vertices v/s MaxSimplices Sparse Matrix and performing collapses operations using the N^2() Algorithm.
*/
class FlagComplexSpMatrix
{
	private:

  	std::unordered_map<int,Vertex> rowToVertex;

  	// Vertices strored as an unordered_set
  	std::unordered_set<Vertex> vertices; 
  //! Stores the 1-simplices(edges) of the original Simplicial Complex.
   	edge_list oneSimplices;
	
	//Unordered set of  removed edges. (to enforce removal from the matrix)
	std::unordered_set<Edge, boost::hash<Edge>> u_set_removed_redges;

	//Unordered set of  dominated edges. (to inforce removal from the matrix)
	std::unordered_set<Edge, boost::hash<Edge>> u_set_dominated_redges; 
 

	//Map from egde to its index
	u_edge_to_idx_map edge_to_index_map;
	//Boolean vector to indicate if the index is critical or not.
	boolVector critical_edge_indicator; 	// critical indicator
	
	//Boolean vector to indicate if the index is critical or not.
	boolVector dominated_edge_indicator; 	// domination indicator

	//! Stores the Map between vertices<B>rowToVertex  and row indices <B>rowToVertex -> row-index</B>.
    /*!
      \code
      MapVertexToIndex = std::unordered_map<Vertex,int>
      \endcode
      So, if the original simplex tree had vertices 0,1,4,5 <br>
      <B>rowToVertex</B> would store : <br>
      \verbatim
      Values =  | 0 | 1 | 4 | 5 | 
      Indices =   0   1   2   3
      \endverbatim
      And <B>vertexToRow</B> would be a map like the following : <br>
      \verbatim
      0 -> 0
      1 -> 1
      4 -> 2
      5 -> 3
      \endverbatim
    */
	MapVertexToIndex vertexToRow;

	//! Stores the number of vertices in the original Simplicial Complex.
    /*!
      This stores the count of vertices (which is also the number of rows in the Matrix).
    */
	std::size_t rows;

	std::size_t numOneSimplices;

	std::size_t numDomEdge;

	//! Stores the Sparse matrix of double values representing the Original Simplicial Complex.
    /*!
      \code
      sparseRowMatrix   = Eigen::SparseMatrix<double, Eigen::RowMajor> ;
      \endcode
      ;
	*/

	sparseRowMatrix* sparse_colpsd_adj_Matrix;       // Stores the collapsed sparse matrix representaion.
	sparseRowMatrix sparseRowAdjMatrix; 			// This is row-major version of the same sparse-matrix, to facilitate easy access to elements when traversing the matrix row-wise.

	

	//! Stores <I>true</I> for dominated rows and  <I>false</I> for undominated rows. 
    /*!
      Initialised to a vector of length equal to the value of the variable <B>rows</B> with all <I>false</I> values.
      Subsequent removal of dominated vertices is reflected by concerned entries changing to <I>true</I> in this vector.
    */
  	boolVector vertDomnIndicator;  //(domination indicator)
	

	boolVector activeIndicator; 	// active indicator
	boolVector contractionIndicator; //(contraction indicator)
	
	//! Stores the indices of the rows to-be checked for domination in the current iteration. 
    /*!
      Initialised with all rows for the first iteration.
      Subsequently once a dominated row is found, its non-dominated neighbhour indices are inserted.
    */
	//doubleQueue rowIterator;

  	doubleQueue rowIterator;

  	//! Stores the indices-pair of the edges to-be checked for domination in the current iteration. 
    /*!
      Initialised with all egdes for the first iteration.
      Subsequently once a dominated row is found, its non-dominated neighbhour indices are inserted. // To be clarified.
    */
	//doubleQueue rowIterator;

  	// Queue of filtered edges, for edge-collapse, the indices of the edges are the row-indices.
  	EdgeFiltQueue filteredEgdeIter;

  	// Vector of filtered edges, for edge-collapse, the indices of the edges are the row-indices.
  	EdgeFiltVector fEgdeVector;

  	// List of non-dominated edges,  the indices of the edges are the vertex lables!!.
  	Filtered_sorted_edge_list criticalCoreEdges;
  	// Stores the indices from the sorted filtered edge vector.
  	// std::set<std::size_t> recurCriticalCoreIndcs; 


	//! Stores <I>true</I> if the current row is inserted in the queue <B>rowIterator<B> otherwise its value is <I>false<I>. 
    /*!
      Initialised to a boolean vector of length equal to the value of the variable <B>rows</B> with all <I>true</I> values.
      Subsequent removal/addition of a row from <B>rowIterator<B> is reflected by concerned entries changing to <I>false</I>/<I>true</I> in this vector.
    */
	boolVector rowInsertIndicator;  //(current iteration row insertion indicator)

	
	//! Map that stores the current status of the edges after the vertex-collapse has been performed. .
    /*!
      \code
      u_edge_map = std::unordered_map<Edge, bool>
      \endcode
     The values an edge can take are true, false; 
     true -> Inserted in filteredEgdeIter;
     false -> Not inserted in filteredEgdeIter;
    */
	u_edge_map edgeStatusMap; 

    //! Map that stores the Reduction / Collapse of vertices.
    /*!
      \code
      Map = std::unordered_map<Vertex,Vertex>
      \endcode
      This is empty to begin with. As and when collapses are done (let's say from dominated vertex <I>v</I> to dominating vertex <I>v'</I>) : <br>
      <B>ReductionMap</B>[<I>v</I>] = <I>v'</I> is entered into the map. <br>
      <I>This does not store uncollapsed vertices. What it means is that say vertex <I>x</I> was never collapsed onto any other vertex. Then, this map <B>WILL NOT</B> have any entry like <I>x</I> -> <I>x</I>.
      Basically, it will have no entry corresponding to vertex <I>x</I> at all. </I> 
    */
	Map ReductionMap;

	bool vertexCollapsed;
	bool edgeCollapsed;
	//Variable to indicate if filtered-edge-collapse has to be performed. 
	int expansion_limit;

	void init()
	{
		rowToVertex.clear();
    	vertexToRow.clear();
    	oneSimplices.clear();
    	ReductionMap.clear();
    	
		vertDomnIndicator.clear();
		rowInsertIndicator.clear();
		rowIterator.push(0);
		rowIterator.pop();

		filteredEgdeIter.push({{0,0},0});
		filteredEgdeIter.pop();
		fEgdeVector.clear();

		rows = 0;
		numDomEdge = 0;
		
		numOneSimplices = 0;
		expansion_limit = 2;

  		vertexCollapsed = false;
  		edgeCollapsed   = false;
	}
	
	//!	Function for computing the sparse-matrix corresponding to the core of the complex. It also prepares the working list filteredEgdeIter for edge collapses
   	void after_vertex_collapse()
	{
     	sparse_colpsd_adj_Matrix   =  new sparseRowMatrix(rows,rows); // Just for debugging purpose.
    	oneSimplices.clear();  
    	if(not filteredEgdeIter.empty())
    		std::cout << "Working list for edge collapses are not empty before the edge-collapse." << std::endl;

		for(int rw = 0 ; rw < rows ; ++rw)
		{
			if(not vertDomnIndicator[rw]) 				//If the current column is not dominated
			{
				auto nbhrs_to_insert = closed_neighbours_row_index(rw); // returns row indices of the non-dominated vertices.
        		for(auto & v: nbhrs_to_insert) {
          			sparse_colpsd_adj_Matrix->insert(rw, v) = 1; // This creates the full matrix
          			if(rw < v) {
          				oneSimplices.push_back({rowToVertex[rw],rowToVertex[v]});
          				filteredEgdeIter.push({{rw,v},1}) ; 
          				// if(rw == v)
          				// std::cout << "Pushed the edge {" << rw << ", " << v <<  "} " << std::endl;
          				edgeStatusMap[{rw,v}] = true;
          			}
        		}
			}			
		}
		// std::cout << "Total number of non-zero elements before domination check are: " << sparse_colpsd_adj_Matrix->nonZeros() << std::endl;
		// std::cout << "Total number of edges for domination check are: " << filteredEgdeIter.size() << std::endl;
    	// std::cout << *sparse_colpsd_adj_Matrix << std::endl;
		return ;
	}

	//! Function to fully compact a particular vertex of the ReductionMap.
    /*!
      It takes as argument the iterator corresponding to a particular vertex pair (key-value) stored in the ReductionMap. <br>
	  It then checks if the second element of this particular vertex pair is present as a first element of some other key-value pair in the map.
	  If no, then the first element of the vertex pair in consideration is fully compact. 
	  If yes, then recursively call fully_compact_this_vertex() on the second element of the original pair in consideration and assign its resultant image as the image of the first element of the original pair in consideration as well. 
    */
	void fully_compact_this_vertex(Map::iterator iter)
	{
		Map::iterator found = ReductionMap.find(iter->second);
		if ( found == ReductionMap.end() )
			return;

		fully_compact_this_vertex(found);
		iter->second = ReductionMap[iter->second];
	}

	//! Function to fully compact the Reduction Map.
    /*!
      While doing strong collapses, we store only the immediate collapse of a vertex. Which means that in one round, vertex <I>x</I> may collapse to vertex <I>y</I>.
      And in some later round it may be possible that vertex <I>y</I> collapses to <I>z</I>. In which case our map stores : <br>
      <I>x</I> -> <I>y</I> and also <I>y</I> -> <I>z</I>. But it really should store :
      <I>x</I> -> <I>z</I> and <I>y</I> -> <I>z</I>. This function achieves the same. <br>
      It basically calls fully_compact_this_vertex() for each entry in the map.
    */
	void fully_compact()
	{
		Map::iterator it = ReductionMap.begin();
		while(it != ReductionMap.end())
		{
			fully_compact_this_vertex(it);
			it++;
		}
	}

	void sparse_strong_vertex_collapse()
	{
 		complete_vertex_domination_check(rowIterator, rowInsertIndicator, vertDomnIndicator); 		// Complete check for rows in rowIterator, rowInsertIndicator is a list of boolean indicator if a vertex is already inserted in the working row_queue (rowIterator)
  		if( not rowIterator.empty())
			sparse_strong_vertex_collapse();
		else
			return ;
  	}

	void complete_vertex_domination_check (doubleQueue& iterator, boolVector& insertIndicator, boolVector& domnIndicator)
	{
	  	double k;
	  	doubleVector nonZeroInnerIdcs;
	    while(not iterator.empty())       // "iterator" contains list(FIFO) of rows to be considered for domination check 
	    { 
	      	k = iterator.front();
	      	iterator.pop();
	      	insertIndicator[k] = false;
	    	if( not domnIndicator[k]) 				// Check if is  already dominated
	    	{ 
		        nonZeroInnerIdcs  = closed_neighbours_row_index(k); 					     
		        for (doubleVector::iterator it = nonZeroInnerIdcs.begin(); it!=nonZeroInnerIdcs.end(); it++) 
		        {
		       		int checkDom = vertex_domination_check(k, *it);   	// "true" for row domination comparison
		        	if( checkDom == 1)                                  	// row k is dominated by *it, k <= *it;
			        {
			            setZero(k, *it);
			            break ;
			        }
		          	else if(checkDom == -1)                 				// row *it is dominated by k, *it <= k;
		            	setZero(*it, k);
			    }         
	    	}
	    }
	}

	bool check_edge_domination(Edge e) // Edge e is the actual edge (u,v). Not the row ids in the matrixs
	{
		auto u = std::get<0>(e);
		auto v = std::get<1>(e);

		auto rw_u = vertexToRow[u];	  		
		auto rw_v = vertexToRow[v];
		auto rw_e = std::make_pair(rw_u,rw_v);
		// std::cout << "The edge {" << u << ", " << v <<  "} is going for domination check." << std::endl;
	    auto commonNeighbours  = closed_common_neighbours_row_index(rw_e); 
	    // std::cout << "And its common neighbours are." << std::endl;
	    // for (doubleVector::iterator it = commonNeighbours.begin(); it!=commonNeighbours.end(); it++) {
	    	// std::cout << rowToVertex[*it] << ", " ;
	    // }
	    //std::cout<< std::endl;
	    if(commonNeighbours.size() > 2) {
	    	if (commonNeighbours.size() == 3)
	    		return true;
	    	else	
		        for (doubleVector::iterator it = commonNeighbours.begin(); it!=commonNeighbours.end(); it++) {	
		        	auto rw_c = *it; // Typecasting
		       		if(rw_c != rw_u and rw_c != rw_v) {
		       			auto neighbours_c = closed_neighbours_row_index(rw_c);
						if(std::includes(neighbours_c.begin(), neighbours_c.end(), commonNeighbours.begin(), commonNeighbours.end())) // If neighbours_c contains the common neighbours.
							return true;
					}         
		       	}
		}
		return false;
	}

	bool check_domination_indicator(Edge e) // The edge should be sorted by the indices and indices are original
	{
		return dominated_edge_indicator[edge_to_index_map[e]];
	}

	std::set<std::size_t> three_clique_indices(std::size_t  crit) {
		std::set<std::size_t> edge_indices;
		
		EdgeFilt fe = fEgdeVector.at(crit);
		Edge e = std::get<0>(fe);
		Vertex u = std::get<0>(e);
		Vertex v = std::get<1>(e);

		// std::cout << "The  current critical edge to re-check criticality with filt value is : {" << u << "," << v << "}; "<< std::get<1>(fe) << std::endl;
	  	auto rw_u = vertexToRow[u];
	  	auto rw_v = vertexToRow[v];
	  	auto rw_critical_edge = std::make_pair(rw_u,rw_v);

	  	doubleVector commonNeighbours  = closed_common_neighbours_row_index(rw_critical_edge); 

  	    if(commonNeighbours.size() > 2) {					     
	        for (doubleVector::iterator it = commonNeighbours.begin(); it!=commonNeighbours.end(); it++) {
	        	auto rw_c = *it;
	        	if(rw_c != rw_u and rw_c != rw_v) {
		        	auto e_with_new_nbhr_v = std::minmax(u,rowToVertex[rw_c]);
			       	auto e_with_new_nbhr_u = std::minmax(v,rowToVertex[rw_c]);
			       	edge_indices.emplace(edge_to_index_map[e_with_new_nbhr_v]);
			       	edge_indices.emplace(edge_to_index_map[e_with_new_nbhr_u]);
		       }
	        }
	    }
	    return edge_indices;

	}

  	void set_edge_critical(std::size_t indx, double filt)
  	{
		// std::cout << "The curent index  with filtration value " << indx << ", " << filt << " is primary critical" << std::endl;
		std::set<std::size_t> effectedIndcs = three_clique_indices(indx);  
   		if(effectedIndcs.size() > 0){
	     	for(auto idx = indx-1; idx > 0 ; idx--) {
				EdgeFilt fec = fEgdeVector.at(idx);
				Edge e = std::get<0>(fec);
				Vertex u = std::get<0>(e);
				Vertex v = std::get<1>(e);
				if ( not critical_edge_indicator.at(idx) ) {	// If idx is not critical so it should be proceses, otherwise it stays in the graph // prev code : recurCriticalCoreIndcs.find(idx) == recurCriticalCoreIndcs.end()
					if( effectedIndcs.find(idx) != effectedIndcs.end()) { 					// If idx is affected
						if(not check_edge_domination(e)) {
							// std::cout << "The curent index is became critical " << idx  << std::endl;
							critical_edge_indicator.at(idx) = true;
							criticalCoreEdges.push_back({filt,u,v});   	
							std::set<std::size_t> inner_effected_indcs = three_clique_indices(idx);  
							for(auto inr_idx = inner_effected_indcs.rbegin(); inr_idx != inner_effected_indcs.rend(); inr_idx++ )
					   		{
					   			if(*inr_idx < idx)
					   				effectedIndcs.emplace(*inr_idx);
					   		}
					   		inner_effected_indcs.clear();							
							// std::cout << "The following edge is critical with filt value: {" << std::get<0>(e) << "," << std::get<1>(e) << "}; " << filt << std::endl;
						}
						else
							u_set_dominated_redges.emplace(std::minmax(vertexToRow[u],vertexToRow[v]));
					}
					else 			// Idx is not affected hence dominated.
						u_set_dominated_redges.emplace(std::minmax(vertexToRow[u],vertexToRow[v]));

				}
			}	
	   			
	   	}	
   		effectedIndcs.clear();	
   		u_set_dominated_redges.clear();

  	} 

  	void critical_core_edges()
	{
	  	std::size_t totEdges = fEgdeVector.size();

	  	std::size_t endIdx = 0;

	  	u_set_removed_redges.clear();
	  	u_set_dominated_redges.clear();
	  	critical_edge_indicator.clear();

	    while( endIdx < totEdges)
	    {
	    	EdgeFilt fec = fEgdeVector.at(endIdx);

		    insert_new_edges(std::get<0>(std::get<0>(fec)), std::get<1>(std::get<0>(fec)),std::get<1>(fec)); // Inserts the edge in the sparse matrix to update the graph (G_i)
		    // cfiltVal = std::get<1>(fEgdeVector.at(endIdx));
		    // std::cout << "The current processing index is " << endIdx << std::endl;

	      	Edge e = std::get<0>(fec);
	      	Vertex u = std::get<0>(e);
			Vertex v = std::get<1>(e);
			edge_to_index_map.emplace(std::minmax(u,v), endIdx);
			critical_edge_indicator.push_back(false);
			dominated_edge_indicator.push_back(false);

			if ( not check_edge_domination(e) )
			{
    			critical_edge_indicator.at(endIdx) = true;
    			dominated_edge_indicator.at(endIdx) = false;
    			criticalCoreEdges.push_back({std::get<1>(fec),u,v});
    			if(endIdx > 1)
    				set_edge_critical(endIdx, std::get<1>(fec));

			}
			else
				dominated_edge_indicator.at(endIdx) = true;
		    endIdx++;
		}

		std::cout << "The total number of critical edges is: " << criticalCoreEdges.size() << std::endl;
		std::cout << "The total number of non-critical edges is: " << totEdges - criticalCoreEdges.size() << std::endl;
}


	int vertex_domination_check( double i, double j) // True for row comparison, false for column comparison
	{
		if(i != j)
		{
			doubleVector Listi = closed_neighbours_row_index(i);
			doubleVector Listj = closed_neighbours_row_index(j);
			if(Listj.size() <= Listi.size())
			{
      			if(std::includes(Listi.begin(), Listi.end(), Listj.begin(), Listj.end())) // Listj is a subset of Listi
					return -1;
			}
			
			else 
				if(std::includes(Listj.begin(), Listj.end(), Listi.begin(), Listi.end())) // Listi is a subset of Listj
					return 1;
		}
		return 0;	
	}

	doubleVector closed_neighbours_row_index(double indx) 	// Returns list of non-zero columns of the particular indx. 
	{													
	  	doubleVector nonZeroIndices; 
	  	Vertex u = indx;
	  	Vertex v;    
	  	// std::cout << "The neighbours of the vertex: " << rowToVertex[u] << " are. " << std::endl;
	  	if(not vertDomnIndicator[indx]) {
	      	for (rowInnerIterator it(sparseRowAdjMatrix, indx); it; ++it) {             // Iterate over the non-zero columns
	      		v = it.index();
	        	if(not vertDomnIndicator[v] and u_set_removed_redges.find(std::minmax(u, v)) == u_set_removed_redges.end() and u_set_dominated_redges.find(std::minmax(u, v)) == u_set_dominated_redges.end()) { // If the vertex v is not dominated and the edge {u,v} is still in the matrix
	            	nonZeroIndices.push_back(it.index());  // inner index, here it is equal to it.columns()
	            	// std::cout << rowToVertex[it.index()] << ", " ;
	        	}
	        }
	        // std::cout << std::endl;
	    }
      	return nonZeroIndices;
	}

	doubleVector closed_common_neighbours_row_index(Edge e) // Returns the list of closed neighbours of the edge :{u,v}.
	{													
	  	doubleVector common; 
	  	doubleVector nonZeroIndices_u; 
		doubleVector nonZeroIndices_v;
	  	double u = std::get<0>(e) ;
	  	double v = std::get<1>(e) ;
		
		nonZeroIndices_u = closed_neighbours_row_index(u);
		nonZeroIndices_v = closed_neighbours_row_index(v);
		std::set_intersection(nonZeroIndices_u.begin(), nonZeroIndices_u.end(), nonZeroIndices_v.begin(), nonZeroIndices_v.end(), std::inserter(common, common.begin()));

	   	return common;
	}

	void setZero(double dominated, double dominating)
	{
  		for(auto & v: closed_neighbours_row_index(dominated))	
	      if(not rowInsertIndicator[v]) // Checking if the row is already inserted	
	      {  
	        rowIterator.push(v); 
	        rowInsertIndicator[v] = true;
	      }
	    vertDomnIndicator[dominated] = true;
  		ReductionMap[rowToVertex[dominated]]    = rowToVertex[dominating];
  		
  		vertexToRow.erase(rowToVertex[dominated]);
  		vertices.erase(rowToVertex[dominated]);
  		rowToVertex.erase(dominated);  
	}
 	
	vertexVector closed_neighbours_vertex_index(double rowIndx) // Returns list of non-zero "vertices" of the particular colIndx. the difference is in the return type
	{
		vertexVector colmns ; 
  		for(auto & v: closed_neighbours_row_index(rowIndx)) // Iterate over the non-zero columns				          
      		colmns.push_back(rowToVertex[v]); 	
      	std::sort(colmns.begin(), colmns.end());
  		return colmns;
	}

	vertexVector vertex_closed_active_neighbours(double rowIndx) // Returns list of all non-zero "vertices" of the particular colIndx which are currently active. the difference is in the return type.
	{
		vertexVector colmns ; 
  		for(auto & v: closed_neighbours_row_index(rowIndx))  // Iterate over the non-zero columns
     		if(not contractionIndicator[v])  					      // Check if the row corresponds to a contracted vertex
      			colmns.push_back(rowToVertex[v]); 			        
      	std::sort(colmns.begin(), colmns.end());	
  		return colmns;
	}
	
	vertexVector closed_all_neighbours_row_index(double rowIndx) // Returns list of all non-zero "vertices" of the particular colIndx whether dominated or not. the difference is in the return type.
	{
		vertexVector colmns ; 
  		for (rowInnerIterator itCol(sparseRowAdjMatrix,rowIndx); itCol; ++itCol)  // Iterate over the non-zero columns
     		colmns.push_back(rowToVertex[itCol.index()]); 			// inner index, here it is equal to it.row()
     	std::sort(colmns.begin(), colmns.end());
  		return colmns;
	}

	void swap_rows(const Vertex & v, const Vertex & w) {	// swap the rows of v and w. Both should be members of the skeleton
		if(membership(v) && membership(w)){
			auto rw_v = vertexToRow[v];
			auto rw_w = vertexToRow[w];
			vertexToRow[v] = rw_w;
			vertexToRow[w] = rw_v;
			rowToVertex[rw_v] = w;
			rowToVertex[rw_w] = v;
		}
	}

public:

	//! Default Constructor
    /*!
      Only initialises all Data Members of the class to empty/Null values as appropriate.
      One <I>WILL</I> have to create the matrix using the Constructor that has an object of the Simplex_tree class as argument.
    */
	
	FlagComplexSpMatrix()
	{
		init();
	}

	FlagComplexSpMatrix(std::size_t expRows)
	{
		init();
		sparseRowAdjMatrix  = sparseRowMatrix(expansion_limit*expRows, expansion_limit*expRows);    	// Initializing sparseRowAdjMatrix, This is a row-major sparse matrix.  
	}

	//! Main Constructor
    /*!
      Argument is an instance of Filtered_sorted_edge_list. <br>
      This is THE function that initialises all data members to appropriate values. <br>
      <B>rowToVertex</B>, <B>vertexToRow</B>, <B>rows</B>, <B>cols</B>, <B>sparseRowAdjMatrix</B> are initialised here.
      <B>vertDomnIndicator</B>, <B>rowInsertIndicator</B> ,<B>rowIterator<B> are initialised by init() function which is called at the begining of this. <br>
    */
	FlagComplexSpMatrix(const size_t & num_vertices, const Filtered_sorted_edge_list  & edge_t) {
			init();
	 		sparseRowAdjMatrix  = sparseRowMatrix(num_vertices, num_vertices);    	// Initializing sparseRowAdjMatrix, This is a row-major sparse matrix.  
			
			for(size_t bgn_idx = 0; bgn_idx < edge_t.size(); bgn_idx++) {
		  		// std::vector<size_t>  s = {std::get<1>(edge_t.at(bgn_idx)), std::get<2>(edge_t.at(bgn_idx))};
		  		// insert_new_edges(std::get<1>(edge_t.at(bgn_idx)), std::get<2>(edge_t.at(bgn_idx)), 1);
		  		fEgdeVector.push_back({{std::get<1>(edge_t.at(bgn_idx)), std::get<2>(edge_t.at(bgn_idx))},std::get<0>(edge_t.at(bgn_idx))}); 
		  	}	
		
		  	// sparseRowAdjMatrix.makeCompressed(); 
	}

	//!	Destructor.
    /*!
      Frees up memory locations on the heap.
    */
	~FlagComplexSpMatrix()
	{
	}

	//!	Function for performing strong collapse.
    /*!
      calls sparse_strong_vertex_collapse(), and
      Then, it compacts the ReductionMap by calling the function fully_compact().
    */
	double strong_vertex_collapse() {
		auto begin_collapse  = std::chrono::high_resolution_clock::now();
		sparse_strong_vertex_collapse();
		vertexCollapsed = true;
		auto end_collapse  = std::chrono::high_resolution_clock::now();
	
		auto collapseTime = std::chrono::duration<double, std::milli>(end_collapse- begin_collapse).count();
		// std::cout << "Time of Collapse : " << collapseTime << " ms\n" << std::endl;
		
		// Now we complete the Reduction Map
		fully_compact();
		//Post processing...
		after_vertex_collapse();
		return collapseTime;
	}
	
	// Performs edge collapse in a decreasing sequence of the filtration value.
	Filtered_sorted_edge_list filtered_edge_collapse() {
		auto begin_collapse  = std::chrono::high_resolution_clock::now();
		critical_core_edges();
		vertexCollapsed = false;
		edgeCollapsed	= true;
		auto end_collapse  = std::chrono::high_resolution_clock::now();
	
		auto collapseTime = std::chrono::duration<double, std::milli>(end_collapse- begin_collapse).count();
		std::cout << "Time of filtered edge Collapse : " << collapseTime << " ms\n" << std::endl;
		//Post processing...
		// after_edge_collapse();
		// std::cout << sparseRowAdjMatrix << std::endl;
		// for(auto idx = criticalCoreEdges.begin(); idx != criticalCoreEdges.end(); idx++ )
		// {	
		// 	std::cout << 
		// }
		return criticalCoreEdges;
	}

	// double strong_vertex_edge_collapse() {
	// 	auto begin_collapse  = std::chrono::high_resolution_clock::now();
	// 	strong_vertex_collapse();
	// 	strong_edge_collapse();
	// 	// strong_vertex_collapse();
	// 	auto end_collapse = std::chrono::high_resolution_clock::now();
	
	// 	auto collapseTime = std::chrono::duration<double, std::milli>(end_collapse- begin_collapse).count();
	// 	return collapseTime;
	// }

	bool membership(const Vertex & v) {
		auto rw = vertexToRow.find(v);
		if(rw != vertexToRow.end())	
			return true;
		else
			return false;
	}
	 
	bool membership(const Edge & e) {  
    	auto u = std::get<0>(e);
        auto v = std::get<1>(e);
	    if(membership(u) && membership(v)) {
	    	auto rw_u = vertexToRow[u];
	    	auto rw_v = vertexToRow[v];
	    	if(rw_u <= rw_v)
	    		for( auto x : closed_neighbours_row_index(rw_v)){ // Taking advantage of sorted lists.
	    			if(rw_u == x)
	    				return true;	
	    			else if(rw_u < x)
	    			 	return false;
	    		}	 
	    	else
	    		for( auto x : closed_neighbours_row_index(rw_u)){ // Taking advantage of sorted lists.
	    			if(rw_v == x)
	    				return true;	
	    			else if(rw_v < x)
	    			 	return false;		
	    		}	 
	    }
		return false;

	}	 
	void insert_vertex(const Vertex & vertex, double filt_val)
	{
		auto rw = vertexToRow.find(vertex);
		if(rw == vertexToRow.end()) {	 
			sparseRowAdjMatrix.insert(rows,rows) = filt_val; 		// Initializing the diagonal element of the adjency matrix corresponding to rw_b.
			vertDomnIndicator.push_back(false);
			rowInsertIndicator.push_back(true);
			contractionIndicator.push_back(false);
			rowIterator.push(rows);
			vertexToRow.insert(std::make_pair(vertex, rows));
			rowToVertex.insert(std::make_pair(rows, vertex)); 
			vertices.emplace(vertex);
			rows++;
		}
    } 

	void insert_new_edges(const Vertex & u, const Vertex & v, double filt_val) // The edge must not be added before, it should be a new edge.
	{	
		insert_vertex(u, filt_val);
		if( u != v) {
			insert_vertex(v, filt_val);
			// std::cout << "Insertion of the edge begins " << u <<", " << v << std::endl;

			auto rw_u = vertexToRow.find(u);
			auto rw_v = vertexToRow.find(v);
			// std::cout << "Inserting the edge " << u <<", " << v << std::endl;
			sparseRowAdjMatrix.insert(rw_u->second,rw_v->second) 	  = filt_val;
			sparseRowAdjMatrix.insert(rw_v->second,rw_u->second) 	  = filt_val;
			oneSimplices.emplace_back(u, v);
			numOneSimplices++;
		}	
	    // else
		// 	std::cout << "Already a member simplex,  skipping..." << std::endl;
       		
	}

			

    std::size_t num_vertices() const {
       return vertices.size();
	}

	//!	Function for returning the ReductionMap.
    /*!
      This is the (stl's unordered) map that stores all the collapses of vertices. <br>
      It is simply returned.
    */
	
	Map reduction_map() const {
		return ReductionMap;
	}
    std::unordered_set<Vertex> vertex_set() const {
    	return vertices;
    }
    sparseRowMatrix collapsed_matrix() const {
    	return *sparse_colpsd_adj_Matrix;
    }

    sparseRowMatrix uncollapsed_matrix() const {
    	return sparseRowAdjMatrix;
    }
    
    edge_list all_edges() const {
    	return oneSimplices;
    }

    vertexVector active_neighbors(const Vertex & v) {  
    	vertexVector nb;
    	auto rw_v = vertexToRow.find(v);
    	if(rw_v != vertexToRow.end())  		
    		nb = vertex_closed_active_neighbours(rw_v->second);
	   	return nb;
    }

    vertexVector neighbors(const Vertex & v) {  
    	vertexVector nb;
    	auto rw_v = vertexToRow.find(v);
    	if(rw_v != vertexToRow.end())  		
    		nb = closed_neighbours_vertex_index(rw_v->second);
    	
    	return nb;
    }

    vertexVector active_relative_neighbors(const Vertex & v, const Vertex & w){
    	std::vector<Vertex> diff; 
    	if(membership(v) && membership(w)){
    		auto nbhrs_v = active_neighbors(v);
    		auto nbhrs_w = active_neighbors(w);
 			std::set_difference(nbhrs_v.begin(), nbhrs_v.end(), nbhrs_w.begin(), nbhrs_w.end(), std::inserter(diff, diff.begin()));	
        }	
        return diff;
    }

  
    void contraction(const Vertex & del, const Vertex & keep){	
		if(del != keep){
			bool del_mem  = membership (del);
			bool keep_mem = membership(keep);
			if( del_mem && keep_mem)
			{
				doubleVector  del_indcs, keep_indcs, diff;
				auto row_del = vertexToRow[del];
				auto row_keep = vertexToRow[keep];
				del_indcs = closed_neighbours_row_index(row_del);
				keep_indcs = closed_neighbours_row_index(row_keep);
				std::set_difference(del_indcs.begin(), del_indcs.end(), keep_indcs.begin(), keep_indcs.end(), std::inserter(diff, diff.begin()));
				for (auto & v : diff) {
					if( v != row_del){
	        			sparseRowAdjMatrix.insert(row_keep,v) = 1;
	        			sparseRowAdjMatrix.insert(v, row_keep) = 1;
	        		}		
				}	
				vertexToRow.erase(del);
				vertices.erase(del);
				rowToVertex.erase(row_del); 
				//setZero(row_del->second, row_keep->second);
			}
			else if(del_mem && not keep_mem)
			{
				vertexToRow.insert(std::make_pair(keep, vertexToRow.find(del)->second));
				rowToVertex[vertexToRow.find(del)->second] = keep;
				vertices.emplace(keep);
				vertices.erase(del);
				vertexToRow.erase(del);

			}
			else
			{
				std::cerr << "The first vertex entered in the method contraction() doesn't exist in the skeleton." <<std::endl;	
				exit(-1); 	
			}
		}	
	}

     void relable(const Vertex & v, const Vertex & w){ // relable v as w.
     	if(membership(v) and v != w){
     		auto rw_v = vertexToRow[v];
      		rowToVertex[rw_v] = w;
      		vertexToRow.insert(std::make_pair(w, rw_v));	 
			vertices.emplace(w);
			vertexToRow.erase(v);
      		vertices.erase(v);
      		// std::cout<< "Re-named the vertex " << v << " as " << w << std::endl;
     	}
     } 

    //Returns the contracted edge. along with the contracted vertex in the begining of the list at {u,u} or {v,v}
	
    void active_strong_expansion(const Vertex & v, const Vertex & w, double filt_val){
		if(membership(v) && membership(w) && v!= w){
			// std::cout << "Strong expansion of the vertex " << v << " and " << w << " begins. " << std::endl;
			auto active_list_v_w = active_relative_neighbors(v,w);
			auto active_list_w_v = active_relative_neighbors(w,v);
			if(active_list_w_v.size() < active_list_v_w.size()){ 	// simulate the contraction of w by expanding the star of v
				for (auto & x : active_list_w_v){
					active_edge_insertion(v,x, filt_val);
					// std::cout << "Inserted the edge " << v << " , " << x  << std::endl;
				}
				swap_rows(v,w);
				// std::cout << "A swap of the vertex " << v << " and " << w << "took place." << std::endl;
			}
			else {  					
				for (auto & y : active_list_v_w){
					active_edge_insertion(w,y,filt_val);
					// std::cout << "Inserted the edge " << w << ", " << y  << std::endl;
				}
			}
			auto rw_v = vertexToRow.find(v);
			contractionIndicator[rw_v->second] = true;
		}
		if(membership(v) && !membership(w)){
			relable(v,w);	
		}
	}
	void active_edge_insertion(const Vertex & v, const Vertex & w, double filt_val){
		insert_new_edges(v,w, filt_val);
		//update_active_indicator(v,w);
	}	

	void print_sparse_skeleton(){
		std::cout << sparseRowAdjMatrix << std::endl;
	}	

};