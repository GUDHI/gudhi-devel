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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Hasse_diagram"
#include <boost/test/unit_test.hpp>

#include <gudhi/reader_utils.h>
#include <gudhi/Hasse_diagram.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>

// standard stuff
#include <iostream>
#include <sstream>
#include <vector>


typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;	
typedef Gudhi::Hasse_diagram::Hasse_diagram<Cell> Hasse_diag;
typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> Hasse_diag_persistence;
typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;	
typedef Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diag_persistence, Field_Zp> Persistent_cohomology;


BOOST_AUTO_TEST_CASE(check_construction) 
{
	//test of a Hasse diagram from one of the examples:
	Cell* A = new Cell(0,0.0);
	Cell* B = new Cell(1,1.0);
	Cell* C = new Cell(0,1.0);
	Cell* D = new Cell(1,1.0);
	Cell* E = new Cell(1,1.0);
	Cell* F = new Cell(2,2.0);
	Cell* G = new Cell(1,1.0);
	Cell* H = new Cell(2,2.0);
	Cell* I = new Cell(0,1.0);
	Cell* J = new Cell(1,1.0);
	Cell* K = new Cell(0,1.0);
	Cell* L = new Cell(1,1.0);
	Cell* M = new Cell(1,1.0);
	Cell* N = new Cell(2,2.0);
	Cell* O = new Cell(1,1.0);
	Cell* P = new Cell(2,2.0);

	//setting up boundaries and coboundaries of cells:
	//for cell A:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell B:  	  
	std::vector< std::pair<Cell*,int> >& boundary_of_B = B->get_boundary();
	boundary_of_B.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_B.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell C:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell D
	std::vector< std::pair<Cell*,int> >& boundary_of_D = D->get_boundary();
	boundary_of_D.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_D.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell E
	std::vector< std::pair<Cell*,int> >& boundary_of_E = E->get_boundary();
	boundary_of_E.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_E.push_back( std::pair<Cell*,int>( I,1 ) );

	//for cell F
	std::vector< std::pair<Cell*,int> >& boundary_of_F = F->get_boundary();
	boundary_of_F.push_back( std::pair<Cell*,int>( B,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( G,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( J,1 ) );  

	//for cell G	
	std::vector< std::pair<Cell*,int> >& boundary_of_G = G->get_boundary();
	boundary_of_G.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_G.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell H
	std::vector< std::pair<Cell*,int> >& boundary_of_H = H->get_boundary();
	boundary_of_H.push_back( std::pair<Cell*,int>( D,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( G,1 ) );

	//for cell I:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell J
	std::vector< std::pair<Cell*,int> >& boundary_of_J = J->get_boundary();
	boundary_of_J.push_back( std::pair<Cell*,int>( I,1 ) );
	boundary_of_J.push_back( std::pair<Cell*,int>( K,1 ) );

	//for cell K:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell L	
	std::vector< std::pair<Cell*,int> >& boundary_of_L = L->get_boundary();
	boundary_of_L.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_L.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell M
	std::vector< std::pair<Cell*,int> >& boundary_of_M = M->get_boundary();
	boundary_of_M.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_M.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell N
	std::vector< std::pair<Cell*,int> >& boundary_of_N = N->get_boundary();
	boundary_of_N.push_back( std::pair<Cell*,int>( J,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( B,1 ) );  

	//for cell O
	std::vector< std::pair<Cell*,int> >& boundary_of_O = O->get_boundary();
	boundary_of_O.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_O.push_back( std::pair<Cell*,int>( C,1 ) );    

	//for cell P
	std::vector< std::pair<Cell*,int> >& boundary_of_P = P->get_boundary();
	boundary_of_P.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( D,1 ) );   

	std::vector< Cell* > vect_of_cells = {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P};

	//Here is a construction of a standard Hasse diagram:
	Hasse_diag hd( vect_of_cells );  
	
	//std::cout << hd.full_signature_of_the_structure();
	
	std::string reference_output = 
"dimension: 0 filtration: 0.000000 position: 0 deleted_: 0 boundary  coBoundary \
( 1 1) ( 3 1) ( 4 1) ( 12 1) dimension: 1 filtration: 1.000000 position: 1 \
deleted_: 0 boundary ( 0 1) ( 2 1)  coBoundary ( 5 1) ( 13 1) dimension: 0 \
filtration: 1.000000 position: 2 deleted_: 0 boundary  coBoundary ( 1 1) ( 3 1) \
( 6 1) ( 14 1) dimension: 1 filtration: 1.000000 position: 3 deleted_: 0 boundary \
( 0 1) ( 2 1)  coBoundary ( 7 1) ( 15 1) dimension: 1 filtration: 1.000000 position: \
4 deleted_: 0 boundary ( 0 1) ( 8 1)  coBoundary ( 5 1) ( 7 1) dimension: 2 filtration: \
2.000000 position: 5 deleted_: 0 boundary ( 1 1) ( 4 1) ( 6 1) ( 9 1)  coBoundary dimension: \
1 filtration: 1.000000 position: 6 deleted_: 0 boundary ( 10 1) ( 2 1)  coBoundary ( 5 1) \
( 7 1) dimension: 2 filtration: 2.000000 position: 7 deleted_: 0 boundary ( 3 1) ( 4 1) \
( 11 1) ( 6 1)  coBoundary dimension: 0 filtration: 1.000000 position: 8 deleted_: 0 boundary  \
coBoundary ( 4 1) ( 9 1) ( 11 1) ( 12 1) dimension: 1 filtration: 1.000000 position: 9 \
deleted_: 0 boundary ( 8 1) ( 10 1)  coBoundary ( 5 1) ( 13 1) dimension: 0 filtration: \
1.000000 position: 10 deleted_: 0 boundary  coBoundary ( 6 1) ( 9 1) ( 11 1) ( 14 1) dimension: \
1 filtration: 1.000000 position: 11 deleted_: 0 boundary ( 10 1) ( 8 1)  coBoundary ( 7 1) ( 15 1) \
dimension: 1 filtration: 1.000000 position: 12 deleted_: 0 boundary ( 0 1) ( 8 1)  coBoundary \
( 13 1) ( 15 1) dimension: 2 filtration: 2.000000 position: 13 deleted_: 0 boundary ( 9 1) \
( 12 1) ( 14 1) ( 1 1)  coBoundary dimension: 1 filtration: 1.000000 position: 14 deleted_: 0 \
boundary ( 10 1) ( 2 1)  coBoundary ( 13 1) ( 15 1) dimension: 2 filtration: 2.000000 position: \
15 deleted_: 0 boundary ( 11 1) ( 14 1) ( 12 1) ( 3 1)  coBoundary ";
	
	BOOST_CHECK( reference_output.compare( hd.full_signature_of_the_structure() ) == 0 );
	
	delete A;
	delete B;
	delete C;
	delete D;
	delete E;
	delete F;
	delete G;
	delete H;
	delete I;
	delete J;
	delete K;
	delete L;
	delete M;
	delete N;
	delete O;
	delete P;
}


BOOST_AUTO_TEST_CASE(Hasse_diagram_from_file_test) 
{
	Hasse_diag hd( "cw_decomposition_of_torus.hasse" );  
	
	//std::cout << hd.full_signature_of_the_structure();
	
	std::string reference_output = 
"dimension: 0 filtration: 0.000000 position: 0 deleted_: 0 boundary  coBoundary ( 4 1) \
( 5 1) ( 8 1) ( 9 1) dimension: 0 filtration: 0.000000 position: 1 deleted_: 0 boundary  \
coBoundary ( 6 1) ( 7 1) ( 8 1) ( 9 1) dimension: 0 filtration: 0.000000 position: 2 deleted_: \
0 boundary  coBoundary ( 4 1) ( 5 1) ( 10 1) ( 11 1) dimension: 0 filtration: 0.000000 position: 3 \
deleted_: 0 boundary  coBoundary ( 6 1) ( 7 1) ( 10 1) ( 11 1) dimension: 1 filtration: 0.000000 position: \
4 deleted_: 0 boundary ( 0 1) ( 2 1)  coBoundary ( 12 1) ( 13 1) dimension: 1 filtration: 0.000000 position: \
5 deleted_: 0 boundary ( 0 1) ( 2 1)  coBoundary ( 14 1) ( 15 1) dimension: 1 filtration: 0.000000 position: \
6 deleted_: 0 boundary ( 1 1) ( 3 1)  coBoundary ( 12 1) ( 13 1) dimension: 1 filtration: 0.000000 position: \
7 deleted_: 0 boundary ( 1 1) ( 3 1)  coBoundary ( 14 1) ( 15 1) dimension: 1 filtration: 0.000000 position: \
8 deleted_: 0 boundary ( 0 1) ( 1 1)  coBoundary ( 12 1) ( 14 1) dimension: 1 filtration: 0.000000 position: \
9 deleted_: 0 boundary ( 0 1) ( 1 1)  coBoundary ( 13 1) ( 15 1) dimension: 1 filtration: 0.000000 position: \
10 deleted_: 0 boundary ( 2 1) ( 3 1)  coBoundary ( 12 1) ( 14 1) dimension: 1 filtration: 0.000000 position: \
11 deleted_: 0 boundary ( 2 1) ( 3 1)  coBoundary ( 13 1) ( 15 1) dimension: 2 filtration: 0.000000 position: \
12 deleted_: 0 boundary ( 4 1) ( 6 1) ( 8 1) ( 10 1)  coBoundary dimension: 2 filtration: 0.000000 position: \
13 deleted_: 0 boundary ( 4 1) ( 6 1) ( 9 1) ( 11 1)  coBoundary dimension: 2 filtration: 0.000000 position: \
14 deleted_: 0 boundary ( 5 1) ( 7 1) ( 8 1) ( 10 1)  coBoundary dimension: 2 filtration: 0.000000 position: \
15 deleted_: 0 boundary ( 5 1) ( 7 1) ( 9 1) ( 11 1)  coBoundary ";

	BOOST_CHECK( reference_output.compare( hd.full_signature_of_the_structure() ) == 0 );
}


BOOST_AUTO_TEST_CASE(Hasse_diagram_insertion_and_removal_test) 
{
	//in this test we will start from one Hasse diagram, we will add a few elements to
	//it, and then remove them. At the end we will check if we have obtained the same 
	//Hasse diagram that we started from.	
	
	//This is our initial structure:
	//	
	//    |C|---G---|D|
	//     |         |
	//     H         F
	//     |         |
	//    |A|---E---|B|
	//
	Cell* A = new Cell(0,0.0);
	Cell* B = new Cell(0,0.0);
	Cell* C = new Cell(0,0.0);
	Cell* D = new Cell(0,0.0);
	
	Cell* E = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_E = E->get_boundary();
	boundary_of_E.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_E.push_back( std::pair<Cell*,int>( B,1 ) ); 
	
	Cell* F = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_F = F->get_boundary();
	boundary_of_F.push_back( std::pair<Cell*,int>( B,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( D,1 ) ); 
	
	Cell* G = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_G = G->get_boundary();
	boundary_of_G.push_back( std::pair<Cell*,int>( C,1 ) );
	boundary_of_G.push_back( std::pair<Cell*,int>( D,1 ) ); 
	
	Cell* H = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_H = H->get_boundary();
	boundary_of_H.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( C,1 ) ); 
	
	std::vector< Cell* > cells = {A,B,C,D,E,F,G,H};
	Hasse_diag hd( cells );
	
	std::string initial_structure = hd.full_signature_of_the_structure();
	
	//now we will add, and subsequently remove a few additional cells. This
	//is the structure we aim for:
	//	
	//    |C|---G---|D|---N---|K|
	//     |         |         |
	//     H    I    F    O    M
	//     |         |         | 
	//    |A|---E---|B|---L---|J|
	//	
	
	Cell* I = new Cell(2,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_I = I->get_boundary();
	boundary_of_I.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_I.push_back( std::pair<Cell*,int>( H,1 ) ); 
	boundary_of_I.push_back( std::pair<Cell*,int>( F,1 ) ); 
	boundary_of_I.push_back( std::pair<Cell*,int>( G,1 ) );	
	hd.add_cell( I );
	
	Cell* J = new Cell(0,0.0);	
	hd.add_cell( J );
	
	Cell* K = new Cell(0,0.0);	
	hd.add_cell( K );
	
	Cell* L = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_L = L->get_boundary();
	boundary_of_L.push_back( std::pair<Cell*,int>( B,1 ) );
	boundary_of_L.push_back( std::pair<Cell*,int>( J,1 ) ); 
	hd.add_cell( L );
	
	Cell* M = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_M = M->get_boundary();
	boundary_of_M.push_back( std::pair<Cell*,int>( J,1 ) );
	boundary_of_M.push_back( std::pair<Cell*,int>( K,1 ) ); 
	hd.add_cell( M );
	
	Cell* N = new Cell(1,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_N = N->get_boundary();
	boundary_of_N.push_back( std::pair<Cell*,int>( D,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( K,1 ) ); 
	hd.add_cell( N );
	
	Cell* O = new Cell(2,0.0);
	std::vector< std::pair<Cell*,int> >& boundary_of_O = O->get_boundary();
	boundary_of_O.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_O.push_back( std::pair<Cell*,int>( F,1 ) ); 
	boundary_of_O.push_back( std::pair<Cell*,int>( M,1 ) ); 
	boundary_of_O.push_back( std::pair<Cell*,int>( N,1 ) );	
	hd.add_cell( O );
	
	
	//and now we will remove the cells I, J, K, L, M, N, O
	hd.remove_cell( O );
	hd.remove_cell( N );
	hd.remove_cell( M );
	hd.remove_cell( L );
	hd.remove_cell( K );
	hd.remove_cell( J );
	hd.remove_cell( I );
	
	//now we need to call the clean_up_the_structure() procedure:
	hd.clean_up_the_structure();
	
	//and now compare results:
	std::string final_structure = hd.full_signature_of_the_structure();
	
	BOOST_CHECK( final_structure.compare( initial_structure ) == 0 );
	
	
	delete A;
	delete B;
	delete C;
	delete D;
	delete E;
	delete F;
	delete G;
	delete H;	
}


BOOST_AUTO_TEST_CASE( write_to_file_test ) 
{
	//test of a Hasse diagram from one of the examples:
	Cell* A = new Cell(0,0.0);
	Cell* B = new Cell(1,1.0);
	Cell* C = new Cell(0,1.0);
	Cell* D = new Cell(1,1.0);
	Cell* E = new Cell(1,1.0);
	Cell* F = new Cell(2,2.0);
	Cell* G = new Cell(1,1.0);
	Cell* H = new Cell(2,2.0);
	Cell* I = new Cell(0,1.0);
	Cell* J = new Cell(1,1.0);
	Cell* K = new Cell(0,1.0);
	Cell* L = new Cell(1,1.0);
	Cell* M = new Cell(1,1.0);
	Cell* N = new Cell(2,2.0);
	Cell* O = new Cell(1,1.0);
	Cell* P = new Cell(2,2.0);

	//setting up boundaries and coboundaries of cells:
	//for cell A:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell B:  	  
	std::vector< std::pair<Cell*,int> >& boundary_of_B = B->get_boundary();
	boundary_of_B.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_B.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell C:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell D
	std::vector< std::pair<Cell*,int> >& boundary_of_D = D->get_boundary();
	boundary_of_D.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_D.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell E
	std::vector< std::pair<Cell*,int> >& boundary_of_E = E->get_boundary();
	boundary_of_E.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_E.push_back( std::pair<Cell*,int>( I,1 ) );

	//for cell F
	std::vector< std::pair<Cell*,int> >& boundary_of_F = F->get_boundary();
	boundary_of_F.push_back( std::pair<Cell*,int>( B,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( G,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( J,1 ) );  

	//for cell G	
	std::vector< std::pair<Cell*,int> >& boundary_of_G = G->get_boundary();
	boundary_of_G.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_G.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell H
	std::vector< std::pair<Cell*,int> >& boundary_of_H = H->get_boundary();
	boundary_of_H.push_back( std::pair<Cell*,int>( D,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( G,1 ) );

	//for cell I:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell J
	std::vector< std::pair<Cell*,int> >& boundary_of_J = J->get_boundary();
	boundary_of_J.push_back( std::pair<Cell*,int>( I,1 ) );
	boundary_of_J.push_back( std::pair<Cell*,int>( K,1 ) );

	//for cell K:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell L	
	std::vector< std::pair<Cell*,int> >& boundary_of_L = L->get_boundary();
	boundary_of_L.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_L.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell M
	std::vector< std::pair<Cell*,int> >& boundary_of_M = M->get_boundary();
	boundary_of_M.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_M.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell N
	std::vector< std::pair<Cell*,int> >& boundary_of_N = N->get_boundary();
	boundary_of_N.push_back( std::pair<Cell*,int>( J,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( B,1 ) );  

	//for cell O
	std::vector< std::pair<Cell*,int> >& boundary_of_O = O->get_boundary();
	boundary_of_O.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_O.push_back( std::pair<Cell*,int>( C,1 ) );    

	//for cell P
	std::vector< std::pair<Cell*,int> >& boundary_of_P = P->get_boundary();
	boundary_of_P.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( D,1 ) );   

	std::vector< Cell* > vect_of_cells = {A,C,I,K,B,D,E,G,J,L,M,O,F,H,N,P};

	//Here is a construction of a standard Hasse diagram:
	Hasse_diag hd( vect_of_cells );  	
	hd.write_to_file( "test_intermediate_output" );
	Hasse_diag hd1( "test_intermediate_output" );  
	
	std::string first_structure = hd.full_signature_of_the_structure();
	std::string second_structure = hd1.full_signature_of_the_structure();

	
	BOOST_CHECK( first_structure.compare( second_structure ) == 0 );
	
	delete A;
	delete B;
	delete C;
	delete D;
	delete E;
	delete F;
	delete G;
	delete H;
	delete I;
	delete J;
	delete K;
	delete L;
	delete M;
	delete N;
	delete O;
	delete P;
}



BOOST_AUTO_TEST_CASE( computations_of_persistence_test ) 
{
	//test of a Hasse diagram from one of the examples:
	Cell* A = new Cell(0,0.0);
	Cell* B = new Cell(1,1.0);
	Cell* C = new Cell(0,1.0);
	Cell* D = new Cell(1,1.0);
	Cell* E = new Cell(1,1.0);
	Cell* F = new Cell(2,2.0);
	Cell* G = new Cell(1,1.0);
	Cell* H = new Cell(2,2.0);
	Cell* I = new Cell(0,1.0);
	Cell* J = new Cell(1,1.0);
	Cell* K = new Cell(0,1.0);
	Cell* L = new Cell(1,1.0);
	Cell* M = new Cell(1,1.0);
	Cell* N = new Cell(2,2.0);
	Cell* O = new Cell(1,1.0);
	Cell* P = new Cell(2,2.0);

	//setting up boundaries and coboundaries of cells:
	//for cell A:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell B:  	  
	std::vector< std::pair<Cell*,int> >& boundary_of_B = B->get_boundary();
	boundary_of_B.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_B.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell C:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell D
	std::vector< std::pair<Cell*,int> >& boundary_of_D = D->get_boundary();
	boundary_of_D.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_D.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell E
	std::vector< std::pair<Cell*,int> >& boundary_of_E = E->get_boundary();
	boundary_of_E.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_E.push_back( std::pair<Cell*,int>( I,1 ) );

	//for cell F
	std::vector< std::pair<Cell*,int> >& boundary_of_F = F->get_boundary();
	boundary_of_F.push_back( std::pair<Cell*,int>( B,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( G,1 ) );
	boundary_of_F.push_back( std::pair<Cell*,int>( J,1 ) );  

	//for cell G	
	std::vector< std::pair<Cell*,int> >& boundary_of_G = G->get_boundary();
	boundary_of_G.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_G.push_back( std::pair<Cell*,int>( C,1 ) );

	//for cell H
	std::vector< std::pair<Cell*,int> >& boundary_of_H = H->get_boundary();
	boundary_of_H.push_back( std::pair<Cell*,int>( D,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( E,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_H.push_back( std::pair<Cell*,int>( G,1 ) );

	//for cell I:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell J
	std::vector< std::pair<Cell*,int> >& boundary_of_J = J->get_boundary();
	boundary_of_J.push_back( std::pair<Cell*,int>( I,1 ) );
	boundary_of_J.push_back( std::pair<Cell*,int>( K,1 ) );

	//for cell K:
	//Nothing needs to be done, coboundaries will be set up automatically.

	//for cell L	
	std::vector< std::pair<Cell*,int> >& boundary_of_L = L->get_boundary();
	boundary_of_L.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_L.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell M
	std::vector< std::pair<Cell*,int> >& boundary_of_M = M->get_boundary();
	boundary_of_M.push_back( std::pair<Cell*,int>( A,1 ) );
	boundary_of_M.push_back( std::pair<Cell*,int>( I,1 ) );  

	//for cell N
	std::vector< std::pair<Cell*,int> >& boundary_of_N = N->get_boundary();
	boundary_of_N.push_back( std::pair<Cell*,int>( J,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_N.push_back( std::pair<Cell*,int>( B,1 ) );  

	//for cell O
	std::vector< std::pair<Cell*,int> >& boundary_of_O = O->get_boundary();
	boundary_of_O.push_back( std::pair<Cell*,int>( K,1 ) );
	boundary_of_O.push_back( std::pair<Cell*,int>( C,1 ) );    

	//for cell P
	std::vector< std::pair<Cell*,int> >& boundary_of_P = P->get_boundary();
	boundary_of_P.push_back( std::pair<Cell*,int>( L,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( O,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( M,1 ) );
	boundary_of_P.push_back( std::pair<Cell*,int>( D,1 ) );   

	//std::vector< Cell* > vect_of_cells = {A,C,I,K,B,D,E,G,J,L,M,O,F,H,N,P};
	std::vector< Cell* > vect_of_cells = {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P};

	//Here is a construction of a standard Hasse diagram:
	Hasse_diag_persistence hd( vect_of_cells );  	
	
	Persistent_cohomology pcoh_hd(hd,true);  
	pcoh_hd.init_coefficients(2);    
	pcoh_hd.compute_persistent_cohomology(0);
	
	std::stringstream output_pers;
	pcoh_hd.output_diagram(output_pers);
	std::string obtained_output = output_pers.str();	
	
	std::string reference_output = "2  0 0 inf \n2  1 1 inf \n2  1 1 inf \n2  2 2 inf \n2  1 1 2 \n2  1 1 2 \n2  1 1 2 \n";
		
	BOOST_CHECK( obtained_output.compare( reference_output ) == 0 );
	
	delete A;
	delete B;
	delete C;
	delete D;
	delete E;
	delete F;
	delete G;
	delete H;
	delete I;
	delete J;
	delete K;
	delete L;
	delete M;
	delete N;
	delete O;
	delete P;
}





BOOST_AUTO_TEST_CASE( computations_of_persistence_test_cubical_complex_conversion ) 
{
	//in this case we create a random cubical complex, convert it to Hasse diagram and
	//compute its persistence.
	std::vector<unsigned> sizes = {5,5,5};
	std::vector< double > data = 
	{0.475925266,0.6909955496,0.2987723371,0.5941234897,0.6302306829,0.7566871834,
	0.0889347408,0.4145348414,0.4161299802,0.0158570393,0.8177982965,0.6363700503,
	0.5724720464,0.2612326392,0.7449052588,0.8596373471,0.9178112024,0.4578601108,
	0.4430730366,0.0967003377,0.9198298205,0.5071376241,0.6443051952,0.4634103079,
	0.4180808044,0.4797491042,0.7428110982,0.5080985199,0.8205517582,0.6370014616,
	0.044818389,0.798747987,0.0293595032,0.0851407051,0.1143798232,0.9264189743,
	0.2878147298,0.8654495201,0.6961601796,0.5064976541,0.285332279,0.8273268931,
	0.7466716764,0.1321406444,0.0380182192,0.1333596667,0.0697457688,0.2454088915,
	0.8949574926,0.4616186123,0.8748250473,0.3335150792,0.840947388,0.1003871823,
	0.2942723529,0.1268202693,0.6855756056,0.7634359002,0.3067176531,0.6026431809,
	0.4286906647,0.9293038943,0.701803914,0.9060462629,0.0807383519,0.311073167,
	0.303490595,0.4808853157,0.0987648794,0.7937564885,0.4426015958,0.9979149867,
	0.5212089573,0.7533499466,0.1688620416,0.1119819144,0.4341951646,0.1917495718,
	0.4823353051,0.4048604742,0.8277260433,0.1623254251,0.9911443421,0.8347346636,
	0.0845580238,0.1582094594,0.7839610698,0.0045175261,0.6061314014,0.3481750281,
	0.1866215509,0.9111697536,0.2863292079,0.3344595125,0.8805206581,0.9500726885,
	0.423807825,0.1526751057,0.2061331687,0.0951848407,0.0805850523,0.0612173504,
	0.1708136885,0.1219345015,0.3520677362,0.2274386305,0.5623495856,0.5673613751,
	0.0658727519,0.0521670314,0.8833213325,0.4828412561,0.9391740526,0.5193006489,
	0.9286763233,0.2625846171,0.6679543864,0.7168741848,0.4050888666,0.6722856015,
	0.1900163866,0.4200847822,0.0461202529,0.2983607783,0.1813348292};
	
	//Now contruct the cibical complex.
	typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    
	Bitmap_cubical_complex b(sizes,data);	
			
    Hasse_diag_persistence* hd = Gudhi::Hasse_diagram::convert_to_Hasse_diagram_persistence<Bitmap_cubical_complex,Cell>( b );       
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diag_persistence, Field_Zp> Persistent_cohomology;

	Persistent_cohomology pcoh(*hd,true);  
	unsigned field_characteristic = 11;
	double min_persistence = 0;
	
	pcoh.init_coefficients(field_characteristic);    
	pcoh.compute_persistent_cohomology(min_persistence);
	
	std::stringstream output_pers;
	pcoh.output_diagram( output_pers );
	std::string obtained_output = output_pers.str();	
  
	std::string reference_output = 
	"11  0 0.00451753 inf \n11  1 0.162325 0.701804 \n11  1 0.0987649 0.519301 \n11  1 0.170814 0.567361 \n11  1 0.287815 0.667954 \n11  1 0.333515 0.685576 \n11  1 0.306718 0.606131 \n11  1 0.19175 0.482335 \n11  1 0.261233 0.506498 \n11  2 0.763436 0.991144 \n11  2 0.701804 0.929304 \n11  1 0.100387 0.306718 \n11  2 0.701804 0.906046 \n11  1 0.303491 0.480885 \n11  0 0.0697458 0.245409 \n11  2 0.783961 0.91117 \n11  1 0.303491 0.428691 \n11  2 0.746672 0.86545 \n11  1 0.45786 0.572472 \n11  0 0.0461203 0.152675 \n11  0 0.015857 0.0851407 \n11  0 0.0612174 0.12682 \n11  2 0.742811 0.798748 \n11  0 0.0293595 0.0851407 \n11  0 0.0380182 0.084558 \n11  0 0.0448184 0.0889347 \n11  0 0.052167 0.0658728 \n11  1 0.420085 0.423808 \n11  0 0.0951848 0.0987649 \n";
	  
	BOOST_CHECK( obtained_output.compare( reference_output ) == 0 );
}
	
	

BOOST_AUTO_TEST_CASE( computations_of_persistence_test_simplex_tree_conversion ) 
{
	//in this case we create a random cubical complex, convert it to Hasse diagram and
	//compute its persistence.
	std::vector< std::vector< double > > points = 
{
{	0.244567377	,	0.4534274482	,	0.44555614	},{	0.0052390727	,	0.602575748	,	0.7414902009	},{	0.5654369846	,	0.582432837	,	0.2053795555	},
{	0.7793295668	,	0.8533353866	,	0.4397430255	},{	0.9557639691	,	0.9192109511	,	0.9259494832	},{	0.3840833744	,	0.082794898	,	0.8565300764	},
{	0.1428351549	,	0.1902121839	,	0.953872889	},{	0.9339357959	,	0.8610937127	,	0.8912499396	},{	0.5875723823	,	0.5118134613	,	0.486443941	},
{	0.0042094528	,	0.5490676572	,	0.3648811274	},{	0.6106395442	,	0.0479015107	,	0.7234964839	},{	0.9003575924	,	0.463473931	,	0.6906190405	},
{	0.8387085642	,	0.5958318673	,	0.946587868	},{	0.4139854845	,	0.7328906141	,	0.1533220368	},{	0.9693074175	,	0.6534106459	,	0.056046945	},
{	0.8216093997	,	0.3319417189	,	0.5898868984	},{	0.6106844896	,	0.8193073997	,	0.5189275641	},{	0.6388283474	,	0.7048335578	,	0.981393236	},
{	0.9783574799	,	0.2257944583	,	0.9884581352	},{	0.7560783883	,	0.2675755005	,	0.1937581748	},{	0.587111552	,	0.4347299838	,	0.3136795729	},
{	0.5738840688	,	0.7039033608	,	0.14652298	},{	0.5724811275	,	0.5281198355	,	0.2073067205	},{	0.8108193642	,	0.3070845357	,	0.5836749915	},
{	0.4393080305	,	0.7803897262	,	0.8373108173	},{	0.4532144398	,	0.5475729688	,	0.2293561632	},{	0.3355917046	,	0.9282942102	,	0.0465660496	},
{	0.8434037087	,	0.5648629228	,	0.9355533635	},{	0.7028795821	,	0.2255625632	,	0.8234450168	},{	0.9987782293	,	0.9843572681	,	0.1736919994	},
{	0.6658957372	,	0.6353726364	,	0.7423913458	},{	0.3802117028	,	0.4781723069	,	0.0636693444	},{	0.2421696743	,	0.837428093	,	0.8620829294	},
{	0.0010670582	,	0.980847826	,	0.6958072949	},{	0.42232705	,	0.226440059	,	0.6451251025	},{	0.8822138815	,	0.3885846001	,	0.3759150384	},
{	0.9030864094	,	0.340339537	,	0.8118740963	},{	0.1867483079	,	0.8130131601	,	0.0197051407	},{	0.252778834	,	0.2873789812	,	0.3754186316	},
{	0.0087783856	,	0.8993070263	,	0.0718508868	},{	0.3907823127	,	0.9211310176	,	0.5048685686	},{	0.5775422677	,	0.0697872201	,	0.4463922526	},
{	0.2269183863	,	0.3263768568	,	0.8177733568	},{	0.5879656055	,	0.9994577868	,	0.2880502767	},{	0.097509719	,	0.4626215256	,	0.381744294	},
{	0.2110773809	,	0.6381706239	,	0.2474126818	},{	0.4968352262	,	0.106845567	,	0.7738761613	},{	0.6734521776	,	0.3040526875	,	0.6202309472	},
{	0.4626756501	,	0.2990835512	,	0.4827392255	},{	0.1334659336	,	0.6301364193	,	0.5039012358	},{	0.8145505874	,	0.0362096264	,	0.9278269906	},
{	0.655572861	,	0.8479117453	,	0.1224969672	},{	0.5491489104	,	0.3026981994	,	0.8143578223	},{	0.1899994526	,	0.6801554773	,	0.6263206936	},
{	0.3362623113	,	0.5062377129	,	0.7144837207	},{	0.6701612703	,	0.1698187401	,	0.9116249306	},{	0.230611085	,	0.4262196971	,	0.8410191524	},
{	0.8942252414	,	0.1383566547	,	0.8868567133	},{	0.075280997	,	0.6843339154	,	0.8580178597	},{	0.4981495854	,	0.8566128612	,	0.4772330665	},
{	0.5481387693	,	0.0153311666	,	0.3579402522	},{	0.8718389694	,	0.2995902314	,	0.4648993181	},{	0.6801459836	,	0.2122102405	,	0.8892606667	},
{	0.2654320858	,	0.3941686389	,	0.8232600491	},{	0.7958464893	,	0.4984564262	,	0.3064013061	},{	0.6451589935	,	0.070217126	,	0.1188122472	},
{	0.8059301747	,	0.6109557769	,	0.0704263314	},{	0.2906281103	,	0.3204455904	,	0.9403931508	},{	0.4659039152	,	0.2630866987	,	0.8247926794	},
{	0.5212232072	,	0.4935525327	,	0.536168742	},{	0.9266397702	,	0.5609459437	,	0.4547862771	},{	0.2246078542	,	0.0940762695	,	0.3796125795	},
{	0.7074530728	,	0.721940404	,	0.365677546	},{	0.7244691269	,	0.0149161613	,	0.3685105206	},{	0.0448947994	,	0.0351540081	,	0.9516720418	},
{	0.8689939133	,	0.4313556731	,	0.4698950397	},{	0.6903580576	,	0.7219779345	,	0.8363126002	},{	0.9568317193	,	0.6625965836	,	0.8763575919	},
{	0.7797430824	,	0.845015499	,	0.4104402452	},{	0.0911103676	,	0.0236776166	,	0.660911154	},{	0.6340847414	,	0.5438042537	,	0.6510421513	},
{	0.8723665436	,	0.7206786426	,	0.7567551201	},{	0.1795819809	,	0.8183786662	,	0.731336494	},{	0.11894441	,	0.0561744536	,	0.5173398806	},
{	0.2410736179	,	0.5577847087	,	0.3582249796	},{	0.9805349875	,	0.0230989989	,	0.4587616941	},{	0.8640490128	,	0.2959632939	,	0.2840696587	},
{	0.421897894	,	0.0624545745	,	0.7313301964	},{	0.863567238	,	0.4127365707	,	0.9378369267	},{	0.0993859142	,	0.7942271791	,	0.2099778573	},
{	0.1999941028	,	0.4041564127	,	0.1558515604	},{	0.1592797104	,	0.4621188694	,	0.3141038383	},{	0.3937685955	,	0.7452657686	,	0.931068995	},
{	0.5295022607	,	0.5408681727	,	0.6619690931	},{	0.7703861897	,	0.8360621934	,	0.0098989184	},{	0.6660422327	,	0.5753107893	,	0.7834934308	},
{	0.7242438002	,	0.2786071568	,	0.2273587526	},{	0.6809390595	,	0.9278106436	,	0.017736956	},{	0.145308326	,	0.3814102353	,	0.0485017023	},
{	0.9698241742	,	0.287944451	,	0.1123728123	}		
};
	//and now let us create the Rips complex:
	double threshold = 0.3;
	unsigned dim_max = 2;//We are only intersed in see the 1 dimensional generator.
	unsigned coef_field = 11;//We will use Z_11 to compute persistence.
	double min_persistence = 0.1;//we are interested only in the intervals of persistence 0.1 or higher.	

	using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
	using Filtration_value = Simplex_tree::Filtration_value;
	using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

	// Construct the Rips complex in a Simplex Tree
	Rips_complex rips_complex_from_file(points, threshold, Gudhi::Euclidean_distance());

	// Construct the Rips complex in a Simplex Tree
	Simplex_tree simplex_tree;
	rips_complex_from_file.create_complex(simplex_tree, dim_max);
	simplex_tree.initialize_filtration();


    //using Persistent_cohomology_simplex_tree = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;	
	//Persistent_cohomology_simplex_tree pcoh_st(simplex_tree);	
	//pcoh_st.init_coefficients(coef_field);
	//pcoh_st.compute_persistent_cohomology(min_persistence);	
	//pcoh_st.output_diagram();
	
	std::string reference_output = 
"11  1 0.297153 inf \n11  1 0.294646 inf \n11  1 0.294634 inf \n11  1 0.263833 inf \n11  1 0.253499 inf \n11  1 0.25342 inf \n11  1 0.236033 inf \n11  1 0.217658 inf \n11  0 0 inf \n11  0 0 inf \n11  0 0 0.271628 \n11  0 0 0.259022 \n11  0 0 0.247532 \n11  0 0 0.24398 \n11  0 0 0.238475 \n11  0 0 0.235796 \n11  0 0 0.234403 \n11  0 0 0.215533 \n11  0 0 0.211884 \n11  0 0 0.210002 \n11  0 0 0.204546 \n11  0 0 0.203221 \n11  0 0 0.201823 \n11  0 0 0.200367 \n11  0 0 0.198641 \n11  0 0 0.197447 \n11  0 0 0.197388 \n11  0 0 0.195782 \n11  0 0 0.195525 \n11  0 0 0.19539 \n11  0 0 0.195028 \n11  0 0 0.193902 \n11  0 0 0.190355 \n11  0 0 0.190173 \n11  0 0 0.189392 \n11  0 0 0.189181 \n11  0 0 0.185264 \n11  0 0 0.183437 \n11  0 0 0.183413 \n11  0 0 0.182412 \n11  0 0 0.180441 \n11  0 0 0.17768 \n11  0 0 0.175111 \n11  0 0 0.173904 \n11  0 0 0.173381 \n11  0 0 0.172836 \n11  0 0 0.172237 \n11  0 0 0.17149 \n11  0 0 0.169415 \n11  0 0 0.167298 \n11  0 0 0.167214 \n11  0 0 0.162647 \n11  0 0 0.158647 \n11  0 0 0.158261 \n11  0 0 0.15752 \n11  0 0 0.156087 \n11  0 0 0.154912 \n11  0 0 0.153474 \n11  0 0 0.152784 \n11  0 0 0.150565 \n11  0 0 0.149812 \n11  0 0 0.14959 \n11  0 0 0.146201 \n11  0 0 0.143821 \n11  0 0 0.143596 \n11  0 0 0.142636 \n11  0 0 0.142306 \n11  0 0 0.14218 \n11  0 0 0.140146 \n11  0 0 0.13831 \n11  0 0 0.13771 \n11  0 0 0.136123 \n11  0 0 0.13587 \n11  0 0 0.135243 \n11  0 0 0.134659 \n11  0 0 0.134171 \n11  0 0 0.133743 \n11  0 0 0.133375 \n11  0 0 0.133197 \n11  0 0 0.132358 \n11  0 0 0.131698 \n11  0 0 0.130077 \n11  0 0 0.128374 \n11  0 0 0.128305 \n11  0 0 0.128273 \n11  0 0 0.126277 \n11  0 0 0.125675 \n11  0 0 0.122605 \n11  0 0 0.119933 \n11  0 0 0.109991 \n11  0 0 0.107953 \n11  0 0 0.105193 \n11  0 0 0.104098 \n";
			
	Hasse_diag_persistence* hd = Gudhi::Hasse_diagram::convert_to_Hasse_diagram_persistence<Simplex_tree,Cell>( simplex_tree );
	Persistent_cohomology pcoh_hd(*hd);  
	pcoh_hd.init_coefficients(coef_field);    
	pcoh_hd.compute_persistent_cohomology(min_persistence);
	std::stringstream ss;
	pcoh_hd.output_diagram(ss);
	std::string current_output = ss.str();
  
	//std::cout << current_output << "\n\n\n" << reference_output << "\n\n\n";
	
	BOOST_CHECK( current_output.compare( reference_output ) == 0 );
}

