#include <gudhi/reader_utils.h>
#include <gudhi/Hasse_diagram.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Persistent_cohomology.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main(int argc, char** argv) 
{ 
	typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;
	//in this example we will construct a CW decomposition of two dimensional 
	//torus:
	//  ______________________
	//  |V1|_E1__|V2|__E2_|V3|     
	//  |--|              |--| 
	//  |  |              |  |     
	//  |  |              |  |
	//  |E6|              |E3|
	//  |  |              |  |
	//  |  |              |  |
	//  |__|______________|__|
	//  |V6|_E5__|V5|__E4_|V4|
	//
	//Here is the corresponding filtration:
	//  ______________________
	//  |0 |_0___|0 |_1___| 1|     
	//  |--|              |--| 
	//  |  |              |  |     
	//  |  |              |  |
	//  |.1|              | 1|
	//  |  |              |  |
	//  |  |              |  |
	//  |__|______________|__|
	//  |.1|__.1_|0 |__0__| 0|
  
  // //cretign cells
  // Cell* V1 = new Cell(0,0.0);
  // Cell* V2 = new Cell(0,0.0);
  // Cell* V3 = new Cell(0,1.0);
  // Cell* V4 = new Cell(0,0.0);
  // Cell* V5 = new Cell(0,0.0);
  // Cell* V6 = new Cell(0,1e-1);
  // Cell* E1 = new Cell(1,0.0);
  // Cell* E2 = new Cell(1,1.0);
  // Cell* E3 = new Cell(1,1.0);
  // Cell* E4 = new Cell(1,0.0);
  // Cell* E5 = new Cell(1,1e-1);
  // Cell* E6 = new Cell(1,1e-1);

  // //setting up boundaries and coboundaries of cells:
  // std::vector< std::pair<Cell*,int> >& boundary_of_E1 = E1->get_boundary();
  // boundary_of_E1.push_back( std::pair<Cell*,int>( V1,1 ) );
  // boundary_of_E1.push_back( std::pair<Cell*,int>( V2,1 ) );

  // std::vector< std::pair<Cell*,int> >& boundary_of_E2 = E2->get_boundary();
  // boundary_of_E2.push_back( std::pair<Cell*,int>( V2,1 ) );
  // boundary_of_E2.push_back( std::pair<Cell*,int>( V3,1 ) );

  // std::vector< std::pair<Cell*,int> >& boundary_of_E3 = E3->get_boundary();
  // boundary_of_E3.push_back( std::pair<Cell*,int>( V3,1 ) );
  // boundary_of_E3.push_back( std::pair<Cell*,int>( V4,1 ) );

  // std::vector< std::pair<Cell*,int> >& boundary_of_E4 = E4->get_boundary();
  // boundary_of_E4.push_back( std::pair<Cell*,int>( V4,1 ) );
  // boundary_of_E4.push_back( std::pair<Cell*,int>( V5,1 ) );

  // std::vector< std::pair<Cell*,int> >& boundary_of_E5 = E5->get_boundary();
  // boundary_of_E5.push_back( std::pair<Cell*,int>( V5,1 ) );
  // boundary_of_E5.push_back( std::pair<Cell*,int>( V6,1 ) );

  // std::vector< std::pair<Cell*,int> >& boundary_of_E6 = E6->get_boundary();
  // boundary_of_E6.push_back( std::pair<Cell*,int>( V6,1 ) );
  // boundary_of_E6.push_back( std::pair<Cell*,int>( V1,1 ) );

  //cretign cells
  Cell* V000 = new Cell(0,0.0);
  Cell* V001 = new Cell(0,0.0);
  Cell* V010 = new Cell(0,1.0);
  Cell* V011 = new Cell(0,1.0);
  Cell* V100 = new Cell(0,1.0);
  Cell* V101 = new Cell(0,1.0);
  Cell* V110 = new Cell(0,0.0);
  Cell* V111 = new Cell(0,0.0);
  
  Cell* E_00 = new Cell(1,1.0);
  Cell* E_01 = new Cell(1,1.0);
  Cell* E_10 = new Cell(1,1.0);
  Cell* E_11 = new Cell(1,1.0);
  Cell* E0_0 = new Cell(1,1.0);
  Cell* E0_1 = new Cell(1,1.0);
  Cell* E1_0 = new Cell(1,1.0);
  Cell* E1_1 = new Cell(1,1.0);
  Cell* E00_ = new Cell(1,0.0);
  Cell* E01_ = new Cell(1,1.0);
  Cell* E10_ = new Cell(1,1.0);
  Cell* E11_ = new Cell(1,0.0);

  Cell* S0__ = new Cell(2,1.0);
  Cell* S1__ = new Cell(2,1.0);
  Cell* S_0_ = new Cell(2,1.0);
  Cell* S_1_ = new Cell(2,1.0);
  /* no upper and lower cell */
  
  //setting up boundaries and coboundaries of cells:
  std::vector< std::pair<Cell*,int> >& boundary_of_E_00 = E_00->get_boundary();
  boundary_of_E_00.push_back( std::pair<Cell*,int>( V000,1 ) );
  boundary_of_E_00.push_back( std::pair<Cell*,int>( V100,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E_01 = E_01->get_boundary();
  boundary_of_E_01.push_back( std::pair<Cell*,int>( V001,1 ) );
  boundary_of_E_01.push_back( std::pair<Cell*,int>( V101,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E_10 = E_10->get_boundary();
  boundary_of_E_10.push_back( std::pair<Cell*,int>( V010,1 ) );
  boundary_of_E_10.push_back( std::pair<Cell*,int>( V110,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E_11 = E_11->get_boundary();
  boundary_of_E_11.push_back( std::pair<Cell*,int>( V011,1 ) );
  boundary_of_E_11.push_back( std::pair<Cell*,int>( V111,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E0_0 = E0_0->get_boundary();
  boundary_of_E0_0.push_back( std::pair<Cell*,int>( V000,1 ) );
  boundary_of_E0_0.push_back( std::pair<Cell*,int>( V010,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E0_1 = E0_1->get_boundary();
  boundary_of_E0_1.push_back( std::pair<Cell*,int>( V001,1 ) );
  boundary_of_E0_1.push_back( std::pair<Cell*,int>( V011,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E1_0 = E1_0->get_boundary();
  boundary_of_E1_0.push_back( std::pair<Cell*,int>( V100,1 ) );
  boundary_of_E1_0.push_back( std::pair<Cell*,int>( V110,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E1_1 = E1_1->get_boundary();
  boundary_of_E1_1.push_back( std::pair<Cell*,int>( V101,1 ) );
  boundary_of_E1_1.push_back( std::pair<Cell*,int>( V111,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E00_ = E00_->get_boundary();
  boundary_of_E00_.push_back( std::pair<Cell*,int>( V000,1 ) );
  boundary_of_E00_.push_back( std::pair<Cell*,int>( V001,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E01_ = E01_->get_boundary();
  boundary_of_E01_.push_back( std::pair<Cell*,int>( V010,1 ) );
  boundary_of_E01_.push_back( std::pair<Cell*,int>( V011,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E10_ = E10_->get_boundary();
  boundary_of_E10_.push_back( std::pair<Cell*,int>( V100,1 ) );
  boundary_of_E10_.push_back( std::pair<Cell*,int>( V101,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_E11_ = E11_->get_boundary();
  boundary_of_E11_.push_back( std::pair<Cell*,int>( V110,1 ) );
  boundary_of_E11_.push_back( std::pair<Cell*,int>( V111,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_S0__ = S0__->get_boundary();
  boundary_of_S0__.push_back( std::pair<Cell*,int>( E00_,1 ) );
  boundary_of_S0__.push_back( std::pair<Cell*,int>( E01_,1 ) );
  boundary_of_S0__.push_back( std::pair<Cell*,int>( E0_0,1 ) );
  boundary_of_S0__.push_back( std::pair<Cell*,int>( E0_1,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_S1__ = S1__->get_boundary();
  boundary_of_S1__.push_back( std::pair<Cell*,int>( E10_,1 ) );
  boundary_of_S1__.push_back( std::pair<Cell*,int>( E11_,1 ) );
  boundary_of_S1__.push_back( std::pair<Cell*,int>( E1_0,1 ) );
  boundary_of_S1__.push_back( std::pair<Cell*,int>( E1_1,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_S_0_ = S_0_->get_boundary();
  boundary_of_S_0_.push_back( std::pair<Cell*,int>( E0_0,1 ) );
  boundary_of_S_0_.push_back( std::pair<Cell*,int>( E0_1,1 ) );
  boundary_of_S_0_.push_back( std::pair<Cell*,int>( E00_,1 ) );
  boundary_of_S_0_.push_back( std::pair<Cell*,int>( E01_,1 ) );

  std::vector< std::pair<Cell*,int> >& boundary_of_S_1_ = S_1_->get_boundary();
  boundary_of_S_1_.push_back( std::pair<Cell*,int>( E1_0,1 ) );
  boundary_of_S_1_.push_back( std::pair<Cell*,int>( E1_1,1 ) );
  boundary_of_S_1_.push_back( std::pair<Cell*,int>( E10_,1 ) );
  boundary_of_S_1_.push_back( std::pair<Cell*,int>( E11_,1 ) );

  
  std::vector< Cell* > vect_of_cells = {
    V000, V001, V010, V011, V100, V101, V110, V111,
    E_00, E_01, E_10, E_11, E0_0, E0_1, E1_0, E1_1, E00_, E01_, E10_, E11_,
    S0__, S1__, S_0_, S_1_
  };
  
  //Here is a construction of a standard Hasse diagram:
  Gudhi::Hasse_diagram::Hasse_diagram<Cell> hd( vect_of_cells );  
  std::cout << "Here is the Hasse diagam : " << std::endl << hd << std::endl;
  
  //Here is a construction of a Hasse_diagram_persistence and computations of 
  //persistent homology of the complex above. You should see both the fundamental
  //classes of torus and the generators of the three squares. 
  Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> hdp( vect_of_cells );  
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology 
  <Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell>, Field_Zp> Persistent_cohomology;

  Persistent_cohomology pcoh(hdp,true);  
  unsigned field_characteristic = 11;
  double min_persistence = 0;
	
  pcoh.init_coefficients(field_characteristic);    
  pcoh.compute_persistent_cohomology(min_persistence);
  
  std::cout << "Here is the persistent homology of the diagram." << std::endl;
  pcoh.output_diagram();
  
  
  return 0;
}
