/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2017  Swansea University UK
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


#include <gudhi/reader_utils.h>
#include <gudhi/Hasse_diagram.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main() 
{    
    //In this example we will show how to insert and remove elements from 
    //a Hasse diagram.
    //We will start from the following one dimensional complex:
    //   
    //  |f|---o---|g|---p---|h|
    //   |                   |
    //   m                   n 
    //   |                   |
    //  |d|                 |e|
    //   |                   |   
    //   k                   l 
    //   |                   |
    //  |a|---i---|b|---j---|c|
    //
    // Then we will remove all the cells except from 
    //  |a|---i---|b|---j---|c|
    //
    // And starting from this we will construct a new complex:
    //
    //  |a|---i---|b|---j---|c|
    //   |                   | 
    //   t                   w 
    //   |                   |
    //  |r|-------u---------|s|   
    //
    // At evey step we will display the data structure and compute 
    // its homology.
    
    typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;	
    typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> Hasse_diag;
    typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;	
	
	Cell* a = new Cell( 0,0. );
	Cell* b = new Cell( 0,0. );
	Cell* c = new Cell( 0,0. );
	Cell* d = new Cell( 0,0. );
	Cell* e = new Cell( 0,0. );
	Cell* f = new Cell( 0,0. );
	Cell* g = new Cell( 0,0. );
	Cell* h = new Cell( 0,0. );
	
	Cell* i = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_i = i->get_boundary();
    boundary_of_i.push_back( std::pair<Cell*,int>( a,1 ) );
    boundary_of_i.push_back( std::pair<Cell*,int>( b,1 ) );
  
	Cell* j = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_j = j->get_boundary();
    boundary_of_j.push_back( std::pair<Cell*,int>( c,1 ) );
    boundary_of_j.push_back( std::pair<Cell*,int>( b,1 ) );
	
	
	Cell* k = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_k = k->get_boundary();
    boundary_of_k.push_back( std::pair<Cell*,int>( a,1 ) );
    boundary_of_k.push_back( std::pair<Cell*,int>( d,1 ) );
    
	Cell* l = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_l = l->get_boundary();
    boundary_of_l.push_back( std::pair<Cell*,int>( c,1 ) );
    boundary_of_l.push_back( std::pair<Cell*,int>( e,1 ) );
	
	Cell* m = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_m = m->get_boundary();
    boundary_of_m.push_back( std::pair<Cell*,int>( d,1 ) );
    boundary_of_m.push_back( std::pair<Cell*,int>( f,1 ) );
    
	Cell* n = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_n = n->get_boundary();
    boundary_of_n.push_back( std::pair<Cell*,int>( e,1 ) );
    boundary_of_n.push_back( std::pair<Cell*,int>( h,1 ) );
	
	Cell* o = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_o = o->get_boundary();
    boundary_of_o.push_back( std::pair<Cell*,int>( f,1 ) );
    boundary_of_o.push_back( std::pair<Cell*,int>( g,1 ) );
	
	Cell* p = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_p = p->get_boundary();
    boundary_of_p.push_back( std::pair<Cell*,int>( g,1 ) );
    boundary_of_p.push_back( std::pair<Cell*,int>( h,1 ) );
	
	std::vector< Cell* > vector_of_initial_cells = {a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p};
	
	//And now we create the Hasse diagram:
	Hasse_diag hd( vector_of_initial_cells );  
	std::cout << "Here is the initial Hasse diagam : " << std::endl << hd << std::endl;

	//Let us now compute homology of this diagram:	
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology
	<Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell>, Field_Zp> Persistent_cohomology;

	Persistent_cohomology pcoh(hd);  
	int field_characteristic = 11;
	double min_persistence = 0;

	std::cout << "And here are the persistence of the Hasse diagram : " << std::endl;

	pcoh.init_coefficients(field_characteristic);    
	pcoh.compute_persistent_cohomology(min_persistence);
	pcoh.output_diagram();
	
	std::cout << "Now we are removing some cells of the Hasse diagram. " << std::endl;
	
	//Let us now remove the cells: k,l,m,n,o,p,d,e,f,g,h.
	hd.remove_cell( k );	
	hd.remove_cell( l );	
	hd.remove_cell( m );	
	hd.remove_cell( n );	
	hd.remove_cell( o );
	hd.remove_cell( p );
	hd.remove_cell( d );	
	hd.remove_cell( e );	
	hd.remove_cell( f );
	hd.remove_cell( g );
	hd.remove_cell( h );
	
	//Please note that deletion of the cell is already done by the clean_up_the_structure() procedure
	//This procedure is called every time the number of deleted cells is above a certain fraction.
	//This procedure also have to be invoken before computing persistence and before further operations
	//on the data structure:
	
	hd.clean_up_the_structure();
	
	//Here is the complex after the removal:
	std::cout << "Here is the Hasse diagram after the removal: " << std::endl << hd << std::endl;
	
	//let us compute its homology:
	std::cout << "Here is the homology of the structure after the removal operations. \n";
	Persistent_cohomology pcoh1(hd);  
	pcoh1.init_coefficients(field_characteristic);    
	pcoh1.compute_persistent_cohomology(min_persistence);
	pcoh1.output_diagram();
	
	
	//Let us now add the new cells to the structure. First we create them:
	Cell* r = new Cell( 0,0. );
	Cell* s = new Cell( 0,0. );
	
	Cell* t = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_t = t->get_boundary();
    boundary_of_t.push_back( std::pair<Cell*,int>( a,1 ) );
    boundary_of_t.push_back( std::pair<Cell*,int>( r,1 ) );
  
	Cell* u = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_u = u->get_boundary();
    boundary_of_u.push_back( std::pair<Cell*,int>( r,1 ) );
    boundary_of_u.push_back( std::pair<Cell*,int>( s,1 ) );
  
    Cell* w = new Cell( 1,0. );
	std::vector< std::pair<Cell*,int> >& boundary_of_w = w->get_boundary();
    boundary_of_w.push_back( std::pair<Cell*,int>( s,1 ) );
    boundary_of_w.push_back( std::pair<Cell*,int>( c,1 ) );
    
    //And then add them to the sturcture:
    hd.add_cell(r);
    hd.add_cell(s);
    hd.add_cell(t);
    hd.add_cell(u);
    hd.add_cell(w);
    
    //After adding up new cells to the structures some arrays using to compute 
    //(persistent) homology may not be valid. To make them valid we call set_up_the_arrays.
    //Note that this is only the case after addition. If there are removals, please 
    //call clean_up_the_structure(), which automatically call set_up_the_arrays().
    hd.set_up_the_arrays();
    
    
    //Let us display the new Hasse diagram:
    std::cout << "Here is the Hasse diagram after addition of the cells : " << std::endl;
    std::cout << hd << std::endl;
  
	//Let us finally compute the homology of it:
	std::cout << "Here is the homology of the structure after the add_cell operations. \n";
	Persistent_cohomology pcoh2(hd);  
	pcoh2.init_coefficients(field_characteristic);    
	pcoh2.compute_persistent_cohomology(min_persistence);
	pcoh2.output_diagram();	
	
	delete a;
	delete b;
	delete c;
	delete i;
	delete j;
	delete r;
	delete s;
	delete t;
	delete u;
	delete w;	

    return 0;
}

