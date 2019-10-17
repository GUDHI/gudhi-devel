/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
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

#include <gudhi/Fake_simplex_tree.h>
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
// #include <Eigen/Dense>

using Fake_simplex_tree 	= Gudhi::Fake_simplex_tree ;
using Vertex 	            = Fake_simplex_tree::Vertex;
using Simplex             = Fake_simplex_tree::Simplex;

using MapVertexToIndex 	  = std::unordered_map<Vertex,int>;
using Map 				  	    = std::unordered_map<Vertex,Vertex>;

using sparseMatrix 		 	= Eigen::SparseMatrix<double> ;
using sparseRowMatrix   = Eigen::SparseMatrix<double, Eigen::RowMajor> ;

using rowInnerIterator 			= sparseRowMatrix::InnerIterator;
using columnInnerIterator 	= sparseMatrix::InnerIterator;

using intVector 	     = std::vector<int>;
using doubleVector 	   = std::vector<double>;
using vertexVector     = std::vector<Vertex>;
using simplexVector    = std::vector<Simplex>;
using boolVector       = std::vector<bool>;

using doubleQueue 	   = std::queue<double>;

//!  Class SparseMsMatrix 
/*!
  The class for storing the Vertices v/s MaxSimplices Sparse Matrix and performing collapses operations using the N^2() Algorithm.
*/
class SparseMsMatrix
{
private:
	//! Stores the vertices of the original Simplicial Complex in unordered fashion.
    /*!
      \code
      vertexVector = std::vector< Vertex >
      \endcode
      This is a vector that stores all the vertices of the Original Simplicial Complex. <br>
      So, if the original simplex tree had vertices 0,1,4,5 <br>
      This would store : <br>
      \verbatim
      Values =  | 0 | 1 | 4 | 5 | 
      Indices =   0   1   2   3
      \endverbatim
    */
  	vertexVector vertex_list;
  	std::unordered_set<Vertex> vertices; // Vertices strored as an unordered_set.

  //! Stores the maximal simplices of the original Simplicial Complex.
    /*!
      \code
      simplexVector = std::vector< Simplex >
      \endcode
      This is a vector that stores all the maximal simplices of the Original Simplicial Complex. <br>
      \endverbatim
    */
 	simplexVector maximal_simplices;
 
	//! Stores the Map between vertices<B>vertex_list  and row indices <B>vertex_list -> row-index</B>.
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
	int rows;
	//int numColpsdRows; // Total number of collapsed rows.
	//! Stores the number of Maximal Simplices in the original Simplicial Complex.
    /*!
      This stores the count of Maximal Simplices (which is also the number of columns in the Matrix).
    */
	int cols;
	//int numColpsdMxml; // Total number of final maximal simplices.

	//! Stores the Sparse matrix of double values representing the Original Simplicial Complex.
    /*!
      \code
      sparseMatrix   = Eigen::SparseMatrix<double> ;
      \endcode
      So after counting the number of rows and num of Maximal simplices, this is initialised as : <br>
      \code
      sparseMxSimplices =  sparseMatrix(rows,numMaxSimplices);
      \endcode
      And filled with columns by the Constructor with a Fake Simplex tree as an argument.
    ;
	sparseMatrix* Sparse*/

	sparseMatrix sparseMxSimplices;
	sparseRowMatrix sparseRowMxSimplices; // This is row-major version of the same sparse-matrix, to facilitate easy access to elements when traversing the matrix row-wise.

	//! Stores <I>true</I> for dominated rows and  <I>false</I> for undominated rows. 
    /*!
      Initialised to a vector of length equal to the value of the variable <B>rows</B> with all <I>false</I> values.
      Subsequent removal of dominated vertices is reflected by concerned entries changing to <I>true</I> in this vector.
    */
  	boolVector vertDomnIndicator;  //(domination indicator)
	//! Stores <I>true</I> for maximal simplex(dominated) columns and  <I>false</I> for a non-maximal(non-dominated) columns. 
    /*!
      Initialised to an vector of length equal to the value of the variable <B>cols</B> with all <I>false</I> values.
      Subsequent removal of Maximal Simplices (caused by removal of vertices) is reflected by concerned entries changing to <I>false</I> in this array.
    */
  	boolVector simpDomnIndicator; //(domination indicator)

	//! Stores the indices of the rows to-be checked for domination in the current iteration. 
    /*!
      Initialised to a queue with all row-indices inserted.
      Subsequently once the row is checked for dominated the row-index is poped out from the queue. A row-index is inserted once again if it is a non-zero element of a dominated column.
    */
  	doubleQueue rowIterator;
  	//! Stores the indices of the colums to-be checked for domination in the current iteration. 
    /*!
      Initialised to an empty queue.
      Subsequently once a dominated row is found, its non-zero column indices are inserted.
    */
	doubleQueue columnIterator;

	//! Stores <I>true</I> if the current row is inserted in the queue <B>rowIterator<B> otherwise its value is <I>false<I>. 
    /*!
      Initialised to a boolean vector of length equal to the value of the variable <B>rows</B> with all <I>true</I> values.
      Subsequent removal/addition of a row from <B>rowIterator<B> is reflected by concerned entries changing to <I>false</I>/<I>true</I> in this vector.
    */
	boolVector rowInsertIndicator;  //(current iteration row insertion indicator)

  //! Stores <I>true</I> if the current column is inserted in the queue <B>columnIterator<B> otherwise its value is <I>false<I>. 
    /*!
      Initialised to a boolean vector of length equal to the value of the variable <B>cols</B> with all <I>false</I> values.
      Subsequent addition/removal of a column in  <B>columnIterator<B> is reflected by concerned entries changing to <I>true</I>/<I>false</I> in this vector.
    */
	boolVector colInsertIndicator; //(current iteration column insertion indicator)
	
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
      <I>Assumption : Assumes that the sparse matrix is formed.</I><SUP>#</SUP> <br>
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
		vertDomnIndicator.clear();
		rowInsertIndicator.clear();
		//rowIterator.clear();
		for (double r=0; r< sparseMxSimplices.innerSize(); ++r) //innerSize is number of rows for column-major sparse Matrix.
		{
			vertDomnIndicator.push_back(false);
			rowInsertIndicator.push_back(true);
			rowIterator.push(r);
		}
		simpDomnIndicator.clear();
		colInsertIndicator.clear();
		//columnIterator.clear();
		for (double c=0; c< sparseMxSimplices.outerSize(); ++c) //outerSize is number of columns for column-major sparse Matrix.
		{
			simpDomnIndicator.push_back(false);
			colInsertIndicator.push_back(false);
			columnIterator.push(0);  //  A naive way to initialize, might be reduntant.
 			columnIterator.pop(); 
		}
		ReductionMap.clear();
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

	void sparse_strong_collapse()
	{
 		complete_domination_check(rowIterator, rowInsertIndicator, vertDomnIndicator, true); 		// Complete check for rows, by setting last argument "true" 
  	complete_domination_check(columnIterator, colInsertIndicator, simpDomnIndicator, false); 	// Complete check for columns, by setting last argument "false"
		if( not rowIterator.empty())
			sparse_strong_collapse();
		else
			return ;
  }

	void complete_domination_check (doubleQueue& iterator, boolVector& insertIndicator, boolVector& domnIndicator, bool which)
	{
  	double k;
  	doubleVector nonZeroInnerIdcs;
    doubleVector nonZeroOuterIdcs;
    while(not iterator.empty())       // "iterator" contains list(FIFO) of rows/columns to be considered for domination check 
    { 
      	k = iterator.front();
      	iterator.pop();
      	insertIndicator[k] = false;

    	if( not domnIndicator[k]) 				// Check only if not already dominated
    	{ 
	        nonZeroInnerIdcs  = read(k, which); 					          // If which == "true", returns the non-zero columns of row k, otherwise if which == "false" returns the non-zero rows of column k.
	        nonZeroOuterIdcs  = read(nonZeroInnerIdcs[0], !which); 	// If which == 'true', returns the non-zero rows of the first non-zero column from previous columns and vice versa... 
	        for (doubleVector::iterator it = nonZeroOuterIdcs.begin(); it!=nonZeroOuterIdcs.end(); it++) 
	        {
	       		int checkDom = pair_domination_check(k, *it, which);   	// "true" for row domination comparison
	        	if( checkDom == 1)                                  	// row k is dominated by *it, k <= *it;
		        {
		            setZero(k, *it, which);
		            break ;
		        }
	          	else if(checkDom == -1)                 				// row *it is dominated by k, *it <= k;
	            	setZero(*it, k, which);
		    }         
    	}
    }
	}

	int pair_domination_check( double i, double j, bool which) // True for row comparison, false for column comparison
	{
		if(i != j)
		{
			doubleVector Listi = read(i, which);
			doubleVector Listj = read(j, which);
			
      if(std::includes(Listi.begin(), Listi.end(), Listj.begin(), Listj.end())) // Listj is a subset of Listi
				return -1;
			else if(std::includes(Listj.begin(), Listj.end(), Listi.begin(), Listi.end())) // Listi is a subset of Listj
				return 1;
		}
		return 0;	
	}

	doubleVector read(double indx, bool which) // Returns list of non-zero rows(which = true)/columns(which = false) of the particular indx.
	{
  	doubleVector nonZeroIndices;     
  	if(which)
  		nonZeroIndices = read<rowInnerIterator, sparseRowMatrix>(sparseRowMxSimplices, simpDomnIndicator,indx);
  	else
  		nonZeroIndices = read<columnInnerIterator, sparseMatrix>(sparseMxSimplices, vertDomnIndicator,indx);
  	return nonZeroIndices;
	}

	void setZero(double dominated, double dominating, bool which)
	{
  	if (which)
  	{	
  		vertDomnIndicator[dominated] = true;
  		ReductionMap[vertex_list[dominated]]    = vertex_list[dominating];
  		setZero<rowInnerIterator, sparseRowMatrix>(sparseRowMxSimplices, simpDomnIndicator, colInsertIndicator, columnIterator, dominated);
     	}
  	else
  	{
   		simpDomnIndicator[dominated] = true;
   		setZero<columnInnerIterator, sparseMatrix>(sparseMxSimplices, vertDomnIndicator, rowInsertIndicator, rowIterator, dominated);
  	}
	}
 	
 	vertexVector readColumn(double colIndx) // Returns list of non-zero vertices of the particular colIndx.
	{
		vertexVector rows ; 
  		for (sparseMatrix::InnerIterator itRow(sparseMxSimplices,colIndx); itRow; ++itRow)  // Iterate over the non-zero columns
     		if(not vertDomnIndicator[itRow.index()])  // Check if the row corresponds to a dominated vertex
      			rows.push_back(vertex_list[itRow.index()]); // inner index, here it is equal to it.row()

  		return rows;
	}

	template<typename type, typename matrix>
  void setZero(const matrix& m, boolVector& domnIndicator, boolVector& insertIndicator, doubleQueue& iterator, double indx)
	{
	    for (type it(m,indx); it; ++it)  // Iterate over the non-zero rows/columns
	      if(not domnIndicator[it.index()] && not insertIndicator[it.index()]) // Checking if the row/column is already dominated(set zero) or inserted	
	      {  
	        iterator.push(it.index());
	        insertIndicator[it.index()] = true;
	      }
	} 
	template <typename type, typename matrix>
	doubleVector read(const matrix& m, const boolVector& domnIndicator, double indx)
	{
		doubleVector nonZeroIndices; 
      	for (type it(m, indx); it; ++it)             // Iterate over the non-zero rows/columns
        	if(not domnIndicator[it.index()])
            	nonZeroIndices.push_back(it.index());  // inner index, here it is equal to it.row()/it.columns()
    
      return nonZeroIndices;
	}
  	
public:

	//! Default Constructor
    /*!
      Only initialises all Data Members of the class to empty/Null values as appropriate.
      One <I>WILL</I> have to create the matrix using the Constructor that has an object of the Simplex_tree class as argument.
    */
	SparseMsMatrix()
	{
		vertex_list.clear();
    maximal_simplices.clear();
		vertexToRow.clear();

		rows = 0;
		cols = 0;
		ReductionMap.clear();
	  
	  vertDomnIndicator.clear();
	  simpDomnIndicator.clear();

	}

	//! Main Constructor
    /*!
      Argument is an instance of Fake_simplex_tree. <br>
      This is THE function that initialises all data members to appropriate values. <br>
      <B>vertex_list</B>, <B>vertexToRow</B>, <B>rows</B>, <B>cols</B>, <B>sparseMxSimplices</B> are initialised here.
      <B>vertDomnIndicator</B>, <B>rowInsertIndicator</B> ,<B>rowIterator<B>,<B>simpDomnIndicator<B>,<B>colInsertIndicator<B> and <B>columnIterator<B> are initialised by init_lists() function which is called at the end of this. <br>
      What this does:
      	1. Populate <B>vertex_list</B> and <B>vertexToRow</B> by going over through the vertices of the Fake_simplex_tree and assign the variable <B>rows</B> = no. of vertices
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
	SparseMsMatrix(Fake_simplex_tree st)
	{
    vertex_list.clear();
		vertexToRow.clear();
  	maximal_simplices.clear();
  	rows = 0;

		auto begin = std::chrono::high_resolution_clock::now();
		
  	maximal_simplices = st.max_simplices();

    //Extracting the unique vertices using unordered map.
    for(const Simplex& s : maximal_simplices)
      for (Vertex v : s)
        vertices.emplace(v);

    //Creating the map from vertex to row-index.
    for(auto v: vertices)
    {
      vertex_list.push_back(v);  // vertex_list is implicitely a map as well, from rows to vertices
      vertexToRow[v] = rows;
      rows++;
    }     
	 
    // vertex_list and vertexToRow have been formed 	
  
	  int numMaxSimplices   = maximal_simplices.size();
		sparseMxSimplices     = sparseMatrix(rows,numMaxSimplices); 			  // Initializing sparseMxSimplices, This is a column-major sparse matrix.  
		sparseRowMxSimplices  = sparseRowMatrix(rows,numMaxSimplices);    	// Initializing sparseRowMxSimplices, This is a row-major sparse matrix.  
    //We have two sparse matrices to have flexibility of iteratating through non-zero rows and columns efficiently.

		cols = 0; 												        // Maintains the current count of maximal simplices
		for(auto simplex: maximal_simplices) 			// Adding each maximal simplices iteratively.
		{
			for(auto vertex : st.simplex_vertex_range(simplex))  // Transformes the current maximal simplex to row-index representation and adds it as a column in the sparse matrices.
			{
        sparseMxSimplices.insert(vertexToRow[vertex],cols) = 1;
        sparseRowMxSimplices.insert(vertexToRow[vertex],cols) = 1;
			}
			cols++;
		}
    
    auto end = std::chrono::high_resolution_clock::now();
		std::cout << "Created the sparse Matrix" << std::endl;
		std::cout << "Time for sparse Matrix formation is : " <<  std::chrono::duration<double, std::milli>(end- begin).count()
              << " ms\n" << std::endl;
		std::cout << "Total number of Initial maximal simplices are: " << cols << std::endl;
		
		sparseMxSimplices.makeCompressed(); 	             //Optional for memory saving
	  sparseRowMxSimplices.makeCompressed(); 
       	
    init_lists();
		
    // std::cout << sparseMxSimplices << std::endl;
	}

	//!	Destructor.
    /*!
      Frees up memory locations on the heap.
    */
	~SparseMsMatrix()
	{
	}

	//!	Function for performing strong collapse.
    /*!
      calls sparse_strong_collapse(), and
      Then, it completes the ReductionMap by calling the function fully_compact().
    */
	void strong_collapse()
	{
		sparse_strong_collapse();
		// Now we complete the Reduction Map
		fully_compact();
	}


	//!	Function for computing the Fake Simplex_tree corresponding to the core of the complex.
    /*!
      First calls strong_collapse(), and then computes the Fake Simplex_tree of the core using the Sparse matrix that we have.
      How does it compute the Fake simplex tree ? <br>
      Goes over all the columns (remaining MaximalSimplices) and for each of them, inserts that simplex <br>
      ['that simplex' means the maximal simplex with all the (remaining) vertices] with all subfaces using the <br>
      <I>insert_simplex_and_subfaces()</I> function from Gudhi's Fake_simplex_tree.
    */
	Fake_simplex_tree collapsed_tree()
	{
		strong_collapse();

		Fake_simplex_tree new_st;
     	sparseMatrix sparseColpsdMxSimplices     =  sparseMatrix(rows,cols); // Just for debudding purpose.
    	int j = 0;
		for(int co = 0 ; co < simpDomnIndicator.size() ; ++co)
		{
			if(not simpDomnIndicator[co]) 				//If the current column is not dominated
			{
				vertexVector mx_simplex_to_insert = readColumn(co); 
				new_st.insert_simplex_and_subfaces(mx_simplex_to_insert); // might find the filtration value here from the original tree	
        		// for(auto v:mx_simplex_to_insert)
          // 			sparseColpsdMxSimplices.insert(vertexToRow[v],j) = 1;
			  	j++; 
      }		
		}
    	// std::cout << sparseColpsdMxSimplices << std::endl;
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

	// template<typename Set> std::set<Set> powerset(const Set& s, size_t n)
	// {
	//     typedef typename Set::const_iterator SetCIt;
	//     typedef typename std::set<Set>::const_iterator PowerSetCIt;
	//     std::set<Set> res;
	//     if(n > 0) {
	//         std::set<Set> ps = powerset(s, n-1);
	//         for(PowerSetCIt ss = ps.begin(); ss != ps.end(); ss++)
	//             for(SetCIt el = s.begin(); el != s.end(); el++) {
	//                 Set subset(*ss);
	//                 subset.insert(*el);
	//                 res.insert(subset);
	//             }
	//         res.insert(ps.begin(), ps.end());
	//     } else
	//         res.insert(Set());
	//     return res;
	// }
	// template<typename Set> std::set<Set> powerset(const Set& s)
	// {
	//     return powerset(s, s.size());
	// }
};