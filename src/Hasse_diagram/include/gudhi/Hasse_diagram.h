/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017 Swansea University UK
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
 
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include <gudhi/Hasse_diagram_cell.h>

#ifndef HASSE_DIAGRAM_H
#define HASSE_DIAGRAM_H


namespace Gudhi {

namespace Hasse_diagram {
	

template <typename Cell_type> class is_before_in_dimension;
	
	
	

/**
 * \class Hasse_diagram
 * \brief Data structure to store Hasse diagrams.
 *
 * \ingroup Hasse_diagram
 *
 * \details 
 * This is a data structure to store a Hasse diagrams. It allows to store a general 
 * chain complexes in a form of Hasse diagrams. It implements insertion and 
 * removal of cells. It allows to store and read Hasse diagrams from files. 
 * 
 * Please note that this class is not suitable to use with Gudhi engine to compute
 * persistent homology. For that purpose, please use the derived class provided in 
 * Hasse_diagram_persistence.h
 *
 * Please refer to \ref Hasse_diagram for examples.
 *
 * The complex is a template class requiring the following parameters:
 * Cell_type - a parameter describing a cell of Hasse diagram. Please refer to Hasse_diagram_cell.h for further details.
 *
 */	
template < typename Cell_type >
class Hasse_diagram
{
public:
	using Cell_range = std::vector< Cell_type* >;

	/**
	 * Default constructor.
	**/ 
	Hasse_diagram(){}
	
	/**
	 * Creating Hasse diagram from a file. The file format is the following:
	 * Number of cells
	 * cell dimension
	 * ids of cell boundary elements followed by the incidence coefficient.
	 * the two lines above are repeated for each cell. 
	 * It is assumed that the id of a cell is its position in the file. 	 
	**/ 	
    Hasse_diagram( const char* filename );
    
    /**
	 * Constructor to create a Hasse diagram from a vector of cells. It is assumed 
	 * that all the cells have boundaries set up. Setting up the coboundaries will
	 * be done in the constructor based on the information about boundaries. 
	**/ 
    Hasse_diagram( const Cell_range& cells_ ):cells(cells_),number_of_deleted_cells(0)
    {
		this->set_up_positions();
		this->set_up_coboundaries();
	}		

    
	/**
	 * After many operation of deleting cells, this->cells vector may became
	 * very fragmented. Also, the complexity of operation using all the iterators
	 * depends on the actual size of a structure (where the deleted elements are still
	 * stored. This procedure remove permanently all the deleted elements. Ideally, 
	 * it should be initialized when the proportion of deleted elements is larger than a 
	 * predefined constant. 
	**/     
    void clean_up_the_structure()
    {
		bool dbg = false	;
		if ( this->number_of_deleted_cells == 0 )return;
		
		if ( dbg )std::cout << "Calling clean_up_the_structure() procedure. \n";	
		//count the number of not deleted cells:
		size_t number_of_non_deleted_cells = this->cells.size() - this->number_of_deleted_cells;
		
		//create a new vector to store the undeleted cells:
		std::vector< Cell_type* > new_cells;
		new_cells.reserve( number_of_non_deleted_cells );
		//fill the new vector in and adjust the new positions.
		//In the same time make sure that the boundary and coboundary vectors 
		//in every cell are valid.
		size_t counter = 0;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			if ( !this->cells[i]->deleted() )
			{							
				new_cells.push_back( this->cells[i] );
				this->cells[i]->position = static_cast<unsigned>( counter );
				this->cells[i]->remove_deleted_elements_from_boundary_and_coboundary();
				++counter;							
			}
			else
			{
				delete this->cells[i];
			}
		}
		this->cells.swap(new_cells);
		this->number_of_deleted_cells = 0;
		if ( dbg )std::cout << "Done with the clean_up_the_structure() procedure. \n";	
	} 
	
	/**
	 * Procedure that allow to add a cell into the structure. This procedure 
	 * automatically fill in coboundaries of boundary elements, so do not 
	 * duplicate it.
	**/ 
	void add_cell( Cell_type* cell )
	{
		cell->position = static_cast<unsigned>( this->cells.size() );
		this->cells.push_back( cell );		
		//we still need to check if cobounadies of boundary elements of this 
		//cell are set up in the correct way:
		for ( size_t bd = 0 ; bd != cell->boundary.size() ; ++bd ) 
		{
			cell->boundary[bd].first->coBoundary.push_back( std::make_pair( cell,cell->boundary[bd].second ) );
		}
	}
	
	/**
	 * Procedure that allow to remove a cell into the structure.
	**/ 
	void remove_cell( Cell_type* cell )
	{
		//if the flag enable_checking_validity_of_complex is set to true,
		//we will check if the cell that is to be deleted do not have 
		//a non deleted cell in the coboundary and if this is the case, we 
		//will print out the warning, since this can potentially be an error. 
		if ( enable_checking_validity_of_complex )
		{
			for ( size_t cbd = 0 ; cbd != cell->coBoundary.size() ; ++cbd )
			{
				if ( !cell->coBoundary[cbd].first->deleted() )
				{
					std::cout << "Warning, you are deleting cell which have non-deleted cells in the coboundary. This may lead inconsistencies in the data structure.\n";
					break;
				}
			}
		}
		
		cell->delete_cell();
		this->number_of_deleted_cells++;
		
		//in case the structure gets too fragmented, we are calling the 
		//to clean it up.
		if ( this->number_of_deleted_cells/(static_cast<double>( this->cells.size() )) > 
			this->proportion_of_removed_cells_that_triggers_reorganization_of_structure )
		{
			this->clean_up_the_structure();
		}
	}
	
	/**
	 * A procedure writng Hasse diagram to file. The Hasse diagram can be later
	 * reconstructed using Hasse_diagram( const char* filename ) constructor.
	**/ 
	void write_to_file( const char* filename );
	
	/**
	 * Writing to a stream operator.
	**/ 	
	friend std::ostream& operator<<( std::ostream& out, const Hasse_diagram< Cell_type >& c )
	{
		for ( size_t i = 0 ; i != c.cells.size() ; ++i )
		{
			//if the cell is deleted, ignore it.
			if ( c.cells[i]->deleted() )continue;
			out << *(c.cells[i]);
		}
		return out;
	}
	
	friend class is_before_in_dimension<Cell_type>;
	
	/**
	 * A basic iterator that iterate through all the cells in the structure. It is the 
	 * user's responsibility to check if the cell is deleted or not. 
	**/ 
	typedef typename std::vector<Cell_type*>::iterator Simple_all_cells_iterator;
	typedef typename std::vector<Cell_type*> Simple_all_cells_iterator_range;
	Simple_all_cells_iterator_range simple_all_cells_iterator_range(){return this->cells;}
	
	/**
	 * Procedure that retuns a cell in the position pos in the vector of cells.
	 * Note that this cell will change after calling clean_up_the_structure()
	 * procedure.
	**/ 
	inline Cell_type* give_me_cell_at_position( size_t pos )
	{
		if ( pos < this->cells.size() )
		{
			return this->cells[pos];
		}
		else
		{
			std::cerr << "Wrong position of a cell in the give_me_cell_at_position function.\n";
			throw "Wrong position of a cell in the give_me_cell_at_position function.\n";
		}		
	}
	
	/**
	 * Function that display a string being a signature of a structure. 
	 * Used mainly for debugging purposes. 
	**/ 
	std::string full_signature_of_the_structure()
	{
		std::string result;
		for ( size_t i = 0 ; i != this->cells.size() ; ++i )
		{
			result += this->cells[i]->full_signature_of_the_structure();
		}		
		return result;
	}	
	
protected:	
	Cell_range cells;
	
	//to check how fragmented the data structure is (as a result of removing cells).
	size_t number_of_deleted_cells; 
	
	/**
	 * This procedure assumes that the boundaries are already set up for all
	 * the cells, and set up the coboundaries based on them.
	**/  
	void set_up_coboundaries();
	
	
	/**
	 * When cells are not constructed by the class, but given from elsewhere,
	 * they may not have the positions being set up. This procedure sets them up. 
	**/ 
	void set_up_positions();
	
	static double proportion_of_removed_cells_that_triggers_reorganization_of_structure;
	
	/**
	* This variable indicate if a warning should be given anytime a cell is 
	* deleted that have nondeleted cell in the coboundary.
	**/ 
	static bool enable_checking_validity_of_complex;
};//Hasse_diagram

template < typename Cell_type >
bool Hasse_diagram<Cell_type>::enable_checking_validity_of_complex = true;


template <typename Cell_type>
double Hasse_diagram<Cell_type>::proportion_of_removed_cells_that_triggers_reorganization_of_structure = 0.5;


template < typename Cell_type >
Hasse_diagram<Cell_type>::Hasse_diagram( const char* filename )
{
	
	//We assume that the cells in the file are enumerated in increasing order. 
	//The idea of the i-th cell in the file is by default i (starting from zero). 
	//Moreover, the cells are provided from low dimensional to high dimensiona.
	//By doing so, we know that when constructing a given cell, all its boundary
	//has already been constructed. 
	//Here is the format of a file:
	//Number of cells
	//cell dimension
	//ids of cell boundary elements followed by the incidence coefficient.
	//Note that coboundary vector will be computed based on boundary vector.
	bool dbg = false;
	std::string line;
	
	this->number_of_deleted_cells = 0;
	
	std::ifstream in( filename );
	if ( !in.good() )
	{
		std::cout << "The file do not exist, program will now terminate.\n";
		throw "The file do not exist, program will now terminate.\n";
	}
	
	std::getline(in, line);
	while ( line[0] == '#' )
	{\
		std::getline(in, line);
	}
	std::stringstream iss(line);
	
	
	unsigned number_of_cells;
	iss >> number_of_cells;
	this->cells.reserve( number_of_cells );	
	
	//create all the cells:
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{		
		this->cells.push_back( new Cell_type() );		
	}
	
	
	std::getline(in, line);
	
	if ( dbg )
	{
		std::cout << "Number of cells : " << number_of_cells << std::endl;		
	}
	
	size_t size_of_last_boundary = 10;//to initially reserve a vector for bounary elements.
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		Cell_type* new_cell = this->cells[i];		
		while ( line[0] == '#' )
		{
			std::getline(in, line);	
		}
		
		iss.str("");
		iss.clear();
		iss << line;		
		iss >> new_cell->position >> new_cell->dimension;
		
		if ( dbg )std::cout << "Position and dimension of the cell : " << new_cell->position << " , " << new_cell->dimension << std::endl;
		
		if ( new_cell->position != i )
		{
			std::cerr << "Wrong numeration of cells in the file. Cell number : " << i << " is marked as : " << new_cell->position << " in the file." << std::endl;
			throw "Wrong numeration of cells in the file.";
		}
		if ( iss.good() )
		{
			//in this case we still have a filtration value to be read
			//from the file.
			iss >> new_cell->filtration;
			if ( dbg )std::cout << "Filtration of the cell : " << new_cell->filtration << std::endl;
		}
		else
		{
			new_cell->filtration = 0;
		}
		
		std::getline(in, line);
		while ( line[0] == '#' )
		{
			std::getline(in, line);
		}		

		iss.str("");
		iss.clear();
		iss << line;
		unsigned cell_id;
		typename Cell_type::Incidence_type incidence_coef;
		std::vector< std::pair< unsigned,typename Cell_type::Incidence_type > > bdry;
		bdry.reserve( size_of_last_boundary );
		if ( dbg )std::cout << "Here are the boundary elements of the cell.\n";
		while ( iss.good() )
		{			
			iss >> cell_id;
			if ( !iss.good() )continue;
						
			if ( cell_id >= i )
			{
				std::cerr << "Wrong format of a file. THe cell number : " << i << " contain in a boundary a cell that has not been introduced yet.\n";
			}
			iss >> incidence_coef;			
			if ( dbg )std::cout << "( " <<  cell_id << " , " << incidence_coef << " ), ";
			bdry.push_back( std::pair< unsigned,typename Cell_type::Incidence_type >(cell_id,incidence_coef) );
		}
				
		size_of_last_boundary = bdry.size();						
		new_cell->boundary.reserve( size_of_last_boundary );
		for ( size_t bd = 0 ; bd != size_of_last_boundary ; ++bd )
		{			
			new_cell->boundary.push_back( std::make_pair(this->cells[ bdry[bd].first ] , bdry[bd].second) );
		}				
		if ( dbg )
		{
			std::cout << "new_cell->boundary.size() : " << new_cell->boundary.size() << std::endl;		
			std::cout << "Done with this cell. \n";
			getchar();
		}
		
		std::getline(in, line);		
		while ( line[0] == '#' )
		{
			std::getline(in, line);			
		}		
	}
	//now once the boundaries are set, we are to set up the coboundaries.	
	this->set_up_coboundaries();	
}


template < typename Cell_type >
void Hasse_diagram<Cell_type>::set_up_coboundaries()
{	
	//first we check the number of coboundary elements for each cell:
	size_t number_of_cells = this->cells.size();
	std::vector< unsigned > sizes_of_coboundary( number_of_cells , 0 );
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		std::vector< std::pair<Cell_type*,typename Cell_type::Incidence_type> > bdry = this->cells[i]->get_boundary();
		for ( size_t bd = 0 ; bd != bdry.size() ; ++bd )
		{
			sizes_of_coboundary[ bdry[bd].first->get_position() ]++;
		}
	}
	
	//now we reserve the space for all coboundaries
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		this->cells[i]->coBoundary.reserve( sizes_of_coboundary[i] );
	}
	
	//and now we set up the coboundaries.
	for ( size_t i = 0 ; i != number_of_cells ; ++i )
	{
		for ( size_t bd = 0 ; bd != this->cells[i]->boundary.size() ; ++bd )
		{			
			this->cells[ this->cells[i]->boundary[bd].first->position ]
			->
			coBoundary.push_back
			                    ( 
			                    std::make_pair(this->cells[i], 
								this->cells[i]->boundary[bd].second)
			                    );
		}
	}
}


template < typename Cell_type >
void Hasse_diagram<Cell_type>::set_up_positions()
{
	for ( size_t i = 0 ; i != this->cells.size() ; ++i )
	{
		this->cells[i]->get_position() = static_cast<unsigned>(i);
	}	
}//set_up_positions

template < typename Cell_type >
void Hasse_diagram<Cell_type>::write_to_file( const char* filename )
{		
	std::ofstream out( filename );
	//If there are any deleted cells, then we need to clean up the structure 
	//first before writing it to a file. The reason for that is because the 
	//file format assumes a continuous enumeration of cells (from zero to the
	//number od cells). This is not satisfied if we have in the structure any
	//deleted elements.	
	if ( this->number_of_deleted_cells != 0 )
	{
		this->clean_up_the_structure();
	}
	//now we can write to a file the number of (non deleted) cells in the structure. 
	out << this->cells.size() << std::endl;
	//and then the rest of the Hasse diagram.
	out << *this;		
	out.close();
}//template < typename Cell_type >









/**
 * This is a function that take any representation that implements Hasse_complex 
 * interface and return vector of Cell_type* based on it. It is used to construct 
 * objects of class Hasse_diagram and Hasse_diagram_persistence
**/ 

template <typename Complex_type , typename Cell_type>
std::vector<Cell_type*> convert_to_vector_of_Cell_type( Complex_type& cmplx )
{
	bool dbg = false;
	
	if ( dbg )
	{
		std::cout << "cmplx.num_simplices() : " << cmplx.num_simplices() << std::endl;
	}
	
	//create vector of cells of suitable length:
	std::vector< Cell_type* > cells_of_Hasse_diag( cmplx.num_simplices() );
	for ( size_t i = 0 ; i != cmplx.num_simplices() ; ++i ) 
	{
		cells_of_Hasse_diag[i] = new Cell_type();
	}
	
	//First we need to assign keys. It is not neccessary for cubical complexes,
	//but in simplex tree this is not done by default.
	std::vector< typename Complex_type::Simplex_key > boundary;
	typename Complex_type::Filtration_simplex_range range = cmplx.filtration_simplex_range();
	size_t pos = 0;
	for ( typename Complex_type::Filtration_simplex_iterator it = range.begin() ; it != range.end() ; ++it )
	{
		cmplx.assign_key( *it, pos );
		++pos;
	}
	
	
	size_t counter = 0;
	for ( typename Complex_type::Filtration_simplex_iterator it = range.begin() ; it != range.end() ; ++it )
	{
		if ( dbg )
		{
			std::cout << "This is cell number : " << counter << std::endl;
		}		
		Cell_type* this_cell = cells_of_Hasse_diag[counter];
		
		this_cell->get_dimension() = static_cast<int>( cmplx.dimension(*it) );		
		this_cell->get_filtration() = static_cast<typename Cell_type::Filtration_type>( cmplx.filtration(*it) );
		
		if ( dbg )
		{
			std::cout << "this_cell->get_dimension() : " << this_cell->get_dimension() << std::endl;
			std::cout << "this_cell->get_filtration() : " << this_cell->get_filtration() << std::endl;
		}
	
		//get the boundary:
		boundary.clear();
		boundary.reserve(10);
		typename Complex_type::Boundary_simplex_range bd_range = cmplx.boundary_simplex_range( *it );
		for ( typename Complex_type::Boundary_simplex_iterator bd = bd_range.begin() ; bd != bd_range.end() ; ++bd )
		{
			//std::cerr << "cmplx.key(*bd) : " << cmplx.key(*bd) << std::endl;
			boundary.push_back( cmplx.key(*bd) );
		}
		
		if ( dbg )
		{
			std::cerr << "boundary.size() : " << boundary.size() << std::endl << "And here are the boundary elements : " << std::endl;
			for ( size_t bd = 0 ; bd != boundary.size() ; ++bd )
			{	
				std::cout << boundary[bd] << " ";
			}
			std::cout << std::endl;
		}
		
		//get the boundary in the Hasse diagram format:
		this_cell->boundary.reserve( boundary.size() );
		typename Cell_type::Incidence_type incidence = 1;
		for ( size_t bd = 0 ; bd != boundary.size() ; ++bd )
		{						
			this_cell->boundary.push_back
			( std::make_pair( cells_of_Hasse_diag[ boundary[bd]  ] , incidence ) );
			incidence *= -1;
		}
		
		if ( dbg )
		{
			std::cout << "this_cell->boundary.size() : " << this_cell->boundary.size() << std::endl;
			std::cerr << "Set up for this cell \n";
			getchar();
		}
		++counter;		
	}	
	return cells_of_Hasse_diag;	
}//convert_to_vector_of_Cell_type



/**
 * This is a function to convert any representation that implements Hasse_complex interface
 * into Hasse diagram
**/ 
template <typename Complex_type , typename Cell_type>
Hasse_diagram<Cell_type>* convert_to_Hasse_diagram( Complex_type& cmplx )
{	
	return new Hasse_diagram<Cell_type>( convert_to_vector_of_Cell_type(cmplx) );	
}//convert_to_Hasse_diagram



}//namespace Hasse_diagram
}//namespace Gudhi

#endif //HASSE_DIAGRAM_H
