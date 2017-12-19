/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Divyansh Pareek
 *
 *    Copyright (C) 2017 INRIA Sophia Antipolis (France)
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

#include <gudhi/Simplex_tree.h>
#include <iostream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <list>
#include <algorithm>

#include <ctime>
#include <fstream>

using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex = Simplex_tree::Vertex;
using Filtration_value = Simplex_tree::Filtration_value;

using vertexVector = std::vector< Vertex >;
using boolVector = std::vector<bool>;
using MapVertexToIndex = std::unordered_map<Vertex,int>;

using Tuple = std::tuple<int,int,bool>;
using List = std::list<Tuple>;

using Map = std::unordered_map<Vertex,Vertex>;
using Vector = std::vector<int>;
using typePairSimplexBool = std::pair< Simplex_tree::Simplex, bool >;

bool tuplecompare(const Tuple &lhs, const Tuple &rhs) // this is kept like this because we want to sort in decreasing order always
{
	return std::get<0>(lhs) > std::get<0>(rhs);
}

//!  Class MsMatrix 
/*!
  The class for storing the Vertices v/s MaxSimplices Matrix and doing collapse operations using that matrix.
*/
class MsMatrix
{
private:
	//! Stores the vertices of the original Simplicial Complex.
    /*!
      \code
      vertexVector = std::vector< Vertex >
      \endcode
      So basically this is a vector that stores all the vertices of the Original Simplicial Complex. <br>
      So, if the original simplex tree had vertices 0,1,4,5 <br>
      This would store : <br>
      \verbatim
      Values =  | 0 | 1 | 4 | 5 | 
      Indices =   0   1   2   3
      \endverbatim
    */
	vertexVector vertex_list;

	//! Stores the Reverse Map between indices and values of the vector <B>vertex_list</B>.
    /*!
      \code
      MapVertexToIndex = std::unordered_map<Vertex,int>
      \endcode
      So, if the original simplex tree had vertices 0,1,4,5 <br>
      <B>vertex_list</B> would store : <br>
      \verbatim
      Values =  | 0 | 1 | 4 | 5 | 
      Indices =   0   1   2   3
      \endverbatim
      And <B>reverse_map</B> would be a map like the following : <br>
      \verbatim
      0 -> 0
      1 -> 1
      4 -> 2
      5 -> 3
      \endverbatim
    */
	MapVertexToIndex reverse_map;

	//! Stores the number of vertices in the original Simplicial Complex.
    /*!
      This stores the count of vertices (which is also the number of rows in the Matrix).
    */
	int rows;

	//! Stores the Matrix of bool values representing the Original Simplicial Complex.
    /*!
      \code
      boolVector = std::vector<bool>
      \endcode
      So after counting the number of rows, this is initialised as : <br>
      \code
      MxSimplices = new boolVector[rows];
      \endcode
      And filled with columns by the Constructor with a Simplex tree as an argument.
    */
	boolVector* MxSimplices;

	//! Stores the number of Maximal Simplices in the original Simplicial Complex.
    /*!
      This stores the count of Maximal Simplices (which is also the number of columns in the Matrix).
    */
	int cols;

	//! Stores <I>true</I> for active rows and  <I>false</I> for deactivated rows. 
    /*!
      Initialised to an array of length equal to the value of the variable <B>rows</B> with all <I>true</I> values.
      Subsequent removal of dominated vertices is reflected by concerned entries changing to <I>false</I> in this array.
    */
	bool* active_rows;

	//! Stores <I>true</I> for active columns and  <I>false</I> for deactivated columns. 
    /*!
      Initialised to an array of length equal to the value of the variable <B>cols</B> with all <I>true</I> values.
      Subsequent removal of Maximal Simplices (caused by removal of vertices) is reflected by concerned entries changing to <I>false</I> in this array.
    */
	bool* active_cols;

	//! Stores the list of vertices in an "appropriate" manner convenient for collapse.
    /*!
      \code
      List = std::list<Tuple>
      Tuple = std::tuple<int,int,bool>
      \endcode
      Stores the list of vertices that is essential for executing the collapse. The tuple represents : <br>
      <I>(no. of MxSimp of which this vertex is a part of , row number in the matrix , does this have to be checked in this step)</I> <br>
      The bool values are all assigned to <I>true</I> in the initialisation steps. (All vertices are candidates for collapse to begin with)
    */
	List vert_indices;

	//! Stores the list of simplices in an "appropriate" manner convenient for collapse.
    /*!
      \code
      List = std::list<Tuple>
      Tuple = std::tuple<int,int,bool>
      \endcode
      Stores the list of simplices that is essential for executing the collapse. The tuple represents : <br>
      <I>(no. of vertices that are a part of this MaxSimplex , column number in the matrix , does this have to be checked in this step)</I> <br>
      The bool values are all assigned to <I>false</I> in the initialisation steps. (All MaxSimplices are not candidates for collapse to begin with because by definition they are maximal)
    */
	List simp_indices;

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

	//! Initialisation of vertices' and maximal simplices' lists.
    /*!
      <I>Assumption : Assumes that the 2D matrix is formed.</I><SUP>#</SUP> <br>
      Initialises the two lists of vertices and maximal simplices. <br>
      These lists are important because the collapse operations are mainly performed on these only.<br>
      Lists are of 3-tuples of the type <I>(int , int , bool)</I>.<br>
      For vertices :<br>
          <I>(no. of MxSimp of which this vertex is a part of , row number in the matrix , true)</I><br>
      For MaximalSimplices :<br>
          <I>(no. of vertices that are a part of this MxSimp , column number in the matrix , false)</I>
    */
	void init_lists()
	{
		vert_indices.clear();
		for(int r = 0 ; r < rows ; ++r)
		{
			int no_simp = 0;
			for(int c = 0 ; c < cols ; ++c)
			{
				if(MxSimplices[r][c])
				{
					no_simp++;
				}
			}
			Tuple add_this_vert(no_simp,r,true);
			vert_indices.push_back(add_this_vert);
		}

		simp_indices.clear();
		for(int c = 0 ; c < cols ; ++c)
		{
			int no_vert = 0;
			for(int r = 0 ; r < rows ; ++r)
			{
				if(MxSimplices[r][c])
				{
					no_vert++;
				}
			}
			Tuple add_this_simp(no_vert,c,false);
			simp_indices.push_back(add_this_simp);
		}

		ReductionMap.clear();
	}

	//! Checks if index1 contains index2.
    /*!
      Checks if the <I>ticks</I> (true bool values) of row/column<SUP>1</SUP> numbered <B>index2</B> are a subset of those of row/column<SUP>1</SUP> numbered <B>index1</B>.<br>
      <SUP>1</SUP> Does this for rows if the value of argument <B>which</B> is <I>true</I> and for columns if value of the argument <B>which</B> is <I>false</I>.<br>
      Returns <I>true</I> if they are indeed a subset, else <I>false</I>.
    */
	bool check(int index1, int index2, bool which)
	{
		bool ret = true;

		if(which)
		{
			for(int c = 0 ; c < cols ; ++c)
			{
				if( active_cols[c] and MxSimplices[index2][c] and not MxSimplices[index1][c] )
				{
					ret = false;
					break;
				}
			}
		}
		else
		{
			for(int r = 0 ; r < rows ; ++r)
			{
				if( active_rows[r] and MxSimplices[r][index2] and not MxSimplices[r][index1] )
				{
					ret = false;
					break;
				}
			}
		}
		return ret;
	}

	//! Analyses a list for possible collapses/contractions.
    /*!
      <I>Assumes that the list is sorted on the first integer values of the tuples in decreasing order.</I><SUP>#</SUP> <br>
      Does this for the list of vertices (<B>vert_indices</B>) if the value of the argument <B>which</B> is <I>true</I>, otherwise does this for the list of simplices (<B>simp_indices</B>). <br>
      Iterates over the list and for each tuple (say <I>tup</I>) where the bool value is <I>true</I> (meaning that this is a possible candidate for removal),
      it goes over all the entries that have appeared before <I>tup</I> in the list (because only entries before this <I>tup</I> are candidates for dominating it : sorted condition),
      and calls the function check() to see if there is actually a dominating-dominated pair. <br>
      <br>
      If yes, then :<br>
      (1) Remove <I>tup</I> from the list. <SUP>**</SUP> <I>This was the reason of using a list : constant time removal of entries</I> <SUP>**</SUP> <br>
      (2) Append the row number of <I>tup</I> (as contained in the second <I>int</I> element of the tuple) to the vector that we will return. <br>
      (3) If we are operating on <B>vert_indices</B> (ie, if <B>which</B> is <I>true</I>) : Add this pair to the <B>ReductionMap</B>. <br>
      <br>
      If no, then change that <I>true</I> value to <I>false</I> in <I>tup</I>.
    */
	Vector analyse_list(bool which)
	{
		List* ptr;
		if(which)
		{
			ptr = &vert_indices;
		}
		else
		{
			ptr = &simp_indices;
		}

		Vector to_return;
		to_return.clear();

		List::iterator iter = ptr->begin();
		while(iter != ptr->end() )
		{
			if( std::get<2>(*iter) )
			{
				bool remove_iter = false;
				int store_iter_index = std::get<1>(*iter);
				for(List::iterator before = ptr->begin() ; before != iter ; ++before)
				{
					int store_before_index = std::get<1>(*before);
					if( check( store_before_index , store_iter_index , which) )
					{
						remove_iter = true;
						if(which)
						{
							Vertex vert_dominating = vertex_list[ store_before_index ];
							Vertex vert_dominated = vertex_list[ store_iter_index ];
							ReductionMap[vert_dominated] = vert_dominating;
						}
						break;
					}
				}
				if(remove_iter)
				{
					to_return.push_back( store_iter_index );
					iter = ptr->erase(iter);
				}
				else
				{
					Tuple new_tuple( std::get<0>(*iter) , store_iter_index , false );
					*iter = new_tuple;
					++iter;
				}
			}
			else
			{
				++iter;
			}
		}
		
		ptr = NULL;
		return to_return;
	}

	//! Deactivates rows/columns that were identified by analyse_list().
    /*!
      Does this for rows if the value of argument <B>which</B> is <I>true</I> and for columns if value of the argument <B>which</B> is <I>false</I>.<br>
      Deactivation means for all elements (say <I>elem</I>) of the vector <B>indices_to_off</B>, it does : <br>
      <B>active_rows</B>[<I>elem</I>] = <I>false</I> or <B>active_cols</B>[<I>elem</I>] = <I>false</I> as indicated by the value of the argument variable <B>which</B>.
    */
	void turn_off(Vector indices_to_off, bool which)
	{
		if(which)
		{
			for(auto to_off : indices_to_off)
			{
				active_rows[to_off] = false;
			}
		}
		else
		{
			for(auto to_off : indices_to_off)
			{
				active_cols[to_off] = false;
			}
		}
	}

	//! Modifies the other list because of changes identified by analyse_list().
    /*!
      Does this for <B>vert_indices</B> if the value of argument <B>which</B> is <I>true</I> and for <B>simp_indices</B> if value of the argument <B>which</B> is <I>false</I>.<br>
      Essentially, analyse_list() identifies and does some deletions from one of the lists. This potentially causes changes in the other list. This function does those changes.<br>
      Let's say analyse_list() analysed the list of vertices(<B>vert_indices</B>) and removed some of them. Now, MaxSimplices that contained one or more of those vertices (which were removed) would be affected.<br>
      So this function would go over all the vertices (say <I>v</I>) that were removed (as contained in the argument vector <B>indices_to_off</B>) and go over all MaxSimplices (say <I>mxs</I>) 
      (contained in the data member <br> <B>simp_indices</B>) and if 
      <B>MxSimplices</B>[<I>v</I>][<I>mxs</I>] {ie, this MaxSimplex named <I>mxs</I> contains <I>v</I>}, then change the tuple corresponding to this <I>mxs</I> from : <br>
      <I>(num , index , true/false) to (num-1 , index , true)</I>.<br>
      <br>
      Similarly, if the list of simplices(<B>simp_indices</B>) would have been changed by analyse_list(), this function would appropriately modify the list of vertices(<B>vert_indices</B>) {and the argument <B>which</B> would be <I>true</I>}.
    */
	void modify(Vector indices_to_off, bool which)
	{
		if(which)
		{
			for(auto to_off : indices_to_off)
			{
				for (List::iterator trav = vert_indices.begin(); trav != vert_indices.end(); ++trav)
				{
					int this_index = std::get<1>(*trav);
					if( MxSimplices[this_index][to_off] )
					{
						Tuple change( std::get<0>(*trav) - 1, this_index , true );
						*trav = change;
					}
				}
			}
		}
		else
		{
			for(auto to_off : indices_to_off)
			{
				for (List::iterator trav = simp_indices.begin(); trav != simp_indices.end(); ++trav)
				{
					int this_index = std::get<1>(*trav);
					if( MxSimplices[to_off][this_index] )
					{
						Tuple change( std::get<0>(*trav) - 1, this_index , true );
						*trav = change;
					}
				}
			}
		}
	}

	//! Performs one step of the entire strong collapse.
    /*!
      <I>init_lists() should be done before this.</I><SUP>#</SUP> <br>
      Returns <I>true</I> if this one step actually caused some reduction. Else, returns <I>false</I>, which indicates that we have reached the core of the complex. <br>
      What it does : <br>
      1. Sort <B>vert_indices</B> on the first integer values of the tuples in decreasing order.
      2. Apply analyse_list() function on the <B>vert_indices</B>.
      3. If no reduction happens (ie analyse_list() returns an empty vector), then it means that we have reached the core of the complex. Hence, return <I>false</I>.
      4. Otherwise, some reduction has happened (and analyse_list() returns a non-empty vector of indices to remove >> call this vector <I>vec_vert</I>).
      5. So then, apply turn_off() on <B>vert_indices</B> using this <I>vec_vert</I> and apply modify() on <br> <B>simp_indices</B> using this <I>vec_vert</I>.
      6. Sort <B>simp_indices</B> on the first integer values of the tuples in decreasing order.
      7. Apply analyse_list() function on the <B>simp_indices</B> ; say it returns a vector <I>vec_simp</I> that it removed.
      8. Then, apply turn_off() on <B>simp_indices</B> using this <I>vec_simp</I> and apply modify() on <br> <B>vert_indices</B> using this <I>vec_simp</I>.
    */
	bool one_step_collapse()
	{
		vert_indices.sort(tuplecompare);

		Vector vert_off = analyse_list(true);
		if(vert_off.empty()) // terminate if no row(vertex) is 'dominated' by some other row(vertex)
		{
			return false;
		}
		turn_off(vert_off,true);

		modify(vert_off,false);

		simp_indices.sort(tuplecompare);

		Vector simp_off = analyse_list(false);
		turn_off(simp_off,false);

		modify(simp_off,true);

		return true;
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

public:

	//! Default Constructor
    /*!
      Only initialises all Data Members of the class to empty/Null values as appropriate.
      One <I>WILL</I> have to create the matrix using the Constructor that has an object of the Simplex_tree class as argument.
    */
	MsMatrix()
	{
		vertex_list.clear();
		reverse_map.clear();

		rows = 0;
		cols = 0;
		MxSimplices = NULL;

		active_rows = NULL;
		active_cols = NULL;

		vert_indices.clear();
		simp_indices.clear();

		ReductionMap.clear();
	}

	//! Main Constructor
    /*!
      Argument is an instance of Simplex_tree. <br>
      This is THE function that initialises all data members to appropriate values. <br>
      <B>vertex_list</B>, <B>reverse_map</B>, <B>rows</B>, <B>cols</B>, <B>MxSimplices</B>, <B>active_rows</B> and <B>active_cols</B> are initialised here.
      <B>vert_indices</B> and <B>simp_indices</B> are initialised by init_lists() function which is called at the end of this. <br>
      What this does:
      	1. Populate <B>vertex_list</B> and <B>reverse_map</B> by going over through the vertices of the Simplex_tree and assign the variable <B>rows</B> = no. of vertices
      	2. Initialise the variable <B>cols</B> to zero and allocate memory from the heap to <B>MxSimplices</B> by doing <br>
      	        <I>MxSimplices = new boolVector[rows];</I>
      	3. Iterate over all simplices [Depth-First-Search fashion] (using Gudhi's <I>complex_simplex_range()</I> ) and for all leaf nodes of the tree [candidates for Maximal Simplices] :<br>
      	Check if there is already a maximal simplex inserted into the matrix that is a coface of the current simplex in concern.<br>
      	If not, insert this simplex as a maximal simplex into the Matrix [candidacy confirmed] and increment the variable <B>cols</B> by one.<br>
      	Else, don't, because it is confirmed that this is not a maximal simplex.
      	4. Initialise <B>active_rows</B> to an array of length equal to the value of the variable <B>rows</B> and all values assigned true. [All vertices are there in the simplex to begin with]
      	5. Initialise <B>active_cols</B> to an array of length equal to the value of the variable <B>cols</B> and all values assigned true. [All maximal simplices are maximal to begin with]
      	6. Calls the private function init_lists().
    */
	MsMatrix(Simplex_tree st)
	{
		vertex_list.clear();
		reverse_map.clear();

		rows = 0;
		for(auto vertex:st.complex_vertex_range() )
		{
			vertex_list.push_back(vertex);
			reverse_map[vertex] = rows;
			rows++;
		}
		// vertex_list and reverse_map have been formed

		MxSimplices = new boolVector[rows];
		cols = 0; // maintains the count of maximal simplices

		for(auto simplex:st.complex_simplex_range() )
		{
			if( not st.has_children(simplex) )
			{
				std::vector<int> ptrs_to_check; // pointers to rows of MxSimplices that have to be checked
				boolVector insert_mxs;
				for(int init = 0 ; init < rows ; ++init) // initialised to all false values
				{
					insert_mxs.push_back(false);
				}
				int count_check = 0;
				for(auto vertex : st.simplex_vertex_range(simplex))
				{
					ptrs_to_check.push_back( reverse_map[vertex] );
					insert_mxs[ reverse_map[vertex] ] = true;
					count_check++;
				}

				bool insert_this = true;
				int j = 0;
				while(j < cols)
				{
					int k;
					for(k = 0 ; k < count_check ; ++k)
					{
						if( not MxSimplices[ ptrs_to_check[k] ][j] )
						{
							break;
						}
					}
					if(k == count_check)
					{
						insert_this = false;
						break;
					}
					j++;
				}

				if(insert_this)
				{
					for(int ins = 0 ; ins < rows ; ++ins)
					{
						MxSimplices[ins].push_back(insert_mxs[ins]);
					}
					cols++;
				}
			}
		}
		// MxSimplices is formed appropriately

		active_rows = new bool[rows];
		for(int dum = 0 ; dum < rows ; ++dum)
		{
			active_rows[dum] = true;
		}
		// initialise all rows to active

		active_cols = new bool[cols];
		for(int dum = 0 ; dum < cols ; ++dum)
		{
			active_cols[dum] = true;
		}
		// initialise all columns to active
		init_lists();
	}

	//!	Destructor.
    /*!
      Frees up memory locations on the heap.
      Specifically, does delete[ ] on :
      	1. <B>active_rows</B>
      	2. <B>active_cols</B>
      	3. <B>MxSimplices</B>
    */
	~MsMatrix()
	{
		delete[] active_rows;
		delete[] active_cols;

		delete[] MxSimplices;
	}

	//!	Function for performing strong collapse.
    /*!
      While one step collapses are possible, it does them and stops when the matrix has reached to the core of the initial complex. <br>
      Then, it completes the ReductionMap by calling the function fully_compact().
    */
	void strong_collapse()
	{
		while( one_step_collapse() )
		{}
		// Now we complete the Reduction Map
		fully_compact();
	}
	//!	Function for computing the Simplex_tree corresponding to the core of the complex.
    /*!
      First calls strong_collapse(), and then computes the Simplex_tree of the core using the Matrix that we have.
      How does it compute the simplex tree ? <br>
      Goes over all the columns (remaining MaximalSimplices) and for each of them, inserts that simplex <br>
      ['that simplex' means the maximal simplex with all the (remaining) vertices] with all subfaces using the <br>
      <I>insert_simplex_and_subfaces()</I> function from Gudhi's Simplex_tree.
    */
	Simplex_tree collapsed_tree()
	{
		strong_collapse();

		Simplex_tree new_st;
		for(int co = 0 ; co < cols ; ++co)
		{
			if(active_cols[co])
			{
				vertexVector mx_simplex_to_insert;
				for(int iter = 0 ; iter < rows ; ++iter)
				{
					if(active_rows[iter] and MxSimplices[iter][co])
					{
						mx_simplex_to_insert.push_back(vertex_list[iter]);
					}
				}
				new_st.insert_simplex_and_subfaces(mx_simplex_to_insert); // might find the filtration value here from the original tree
			}
		}
		return new_st;
	}

	//!	Function for returning the ReductionMap.
    /*!
      This is the (stl's unordered) map that stores all the collapses of vertices. <br>
      It is simply returned.
    */
	Map reduction_map()
	{
		return ReductionMap;
	}
};