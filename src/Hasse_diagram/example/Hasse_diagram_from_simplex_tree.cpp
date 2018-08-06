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

#include <gudhi/reader_utils.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main() 
{    
    //In this example we will first construct a Rips complex of points sampled
    //from a circle of a radius 1.
    
    //first the suitable typedefs:
    //for simplex tree:
    using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
	using Filtration_value = Simplex_tree::Filtration_value;
	using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
	using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
	using Persistent_cohomology_simplex_tree = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;	
	//for Hasse diagrams:
	typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;	
    typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> Hasse_diag;
	using Persistent_cohomology_Hasse_diagram = Gudhi::persistent_cohomology::Persistent_cohomology<Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell>, Field_Zp>;
	
	//Here is the point cloud of the circle:
	std::vector< std::vector<double> > circle_pt_cloud = 
	{
		{0,1},{0.0998334166,0.9950041653},{0.1986693308,0.9800665778},{0.2955202067,0.9553364891},
		{0.3894183423,0.921060994},{0.4794255386,0.8775825619},{0.5646424734,0.8253356149},{0.6442176872,0.7648421873},
		{0.7173560909,0.6967067093},{0.7833269096,0.6216099683},{0.8414709848,0.5403023059},{0.8912073601,0.4535961214},
		{0.932039086,0.3623577545},{0.9635581854,0.2674988286},{0.98544973,0.1699671429},{0.9974949866,0.0707372017},
		{0.999573603,-0.0291995223},{0.9916648105,-0.1288444943},{0.9738476309,-0.2272020947},{0.9463000877,-0.3232895669},
		{0.9092974268,-0.4161468365},{0.8632093666,-0.5048461046},{0.8084964038,-0.5885011173},{0.7457052122,-0.6662760213},
		{0.6754631806,-0.7373937155},{0.5984721441,-0.8011436155},{0.5155013718,-0.8568887534},{0.4273798802,-0.904072142},
		{0.3349881502,-0.9422223407},{0.2392493292,-0.9709581651},{0.1411200081,-0.9899924966},{0.0415806624,-0.9991351503},
		{-0.0583741434,-0.9982947758},{-0.1577456941,-0.9874797699},{-0.255541102,-0.9667981926},{-0.3507832277,-0.9364566873},
		{-0.4425204433,-0.8967584163},{-0.5298361409,-0.8481000317},{-0.6118578909,-0.7909677119},{-0.6877661592,-0.7259323042},
		{-0.7568024953,-0.6536436209},{-0.8182771111,-0.5748239465},{-0.8715757724,-0.4902608213},{-0.9161659367,-0.4007991721},
		{-0.9516020739,-0.30733287},{-0.9775301177,-0.2107957994},{-0.9936910036,-0.1121525269},{-0.9999232576,-0.0123886635},
		{-0.9961646088,0.0874989834},{-0.9824526126,0.1865123694},{-0.9589242747,0.2836621855},{-0.9258146823,0.3779777427},
		{-0.8834546557,0.4685166713},{-0.8322674422,0.5543743362},{-0.7727644876,0.6346928759},{-0.7055403256,0.7086697743},
		{-0.6312666379,0.7755658785},{-0.5506855426,0.8347127848},{-0.4646021794,0.8855195169},{-0.3738766648,0.9274784307},
		{-0.2794154982,0.9601702867},{-0.1821625043,0.9832684384},{-0.0830894028,0.996542097},{0.0168139005,0.9998586364}
	};
	
	//and now let us create the Rips complex:
	float threshold = 1;
	int dim_max = 2;//We are only intersed in see the 1 dimensional generator.
	int coef_field = 11;//We will use Z_11 to compute persistence.
	float min_persistence = 0.1;//we are interested only in the intervals of persistence 0.1 or higher.
	Rips_complex rips_complex_from_file(circle_pt_cloud, threshold, Gudhi::Euclidean_distance());

	// Construct the Rips complex in a Simplex Tree
	Simplex_tree simplex_tree;

	rips_complex_from_file.create_complex(simplex_tree, dim_max);

	// Sort the simplices in the order of the filtration
	simplex_tree.initialize_filtration();

	//Here is the computation of persistence of the simplex tree:
	Persistent_cohomology_simplex_tree pcoh_st(simplex_tree);	
	pcoh_st.init_coefficients(coef_field);
	pcoh_st.compute_persistent_cohomology(min_persistence);	
	std::cout << "Here is the persistence of the simplex tree: " << std::endl;
	pcoh_st.output_diagram();
	//The output you should see after executing those lines is:
	//The complex contains 3705 simplices 
    //and has dimension 2 
	//11  0 0 inf 
	//11  1 0.0999583 inf 
	
	//Now let us convert the simplex tree into the Hasse diagram and compute
	//persistence of the Hasse diagram:
    Hasse_diag* hd = Gudhi::Hasse_diagram::convert_to_Hasse_diagram_persistence<Simplex_tree,Cell>( simplex_tree );

	//And compute persistence of Hasse diagram with the same parametes as 
	//for the simplex tree above:	
	Persistent_cohomology_Hasse_diagram pcoh_hd(*hd);  
	pcoh_hd.init_coefficients(coef_field);    
	pcoh_hd.compute_persistent_cohomology(min_persistence);
	std::cout << "Here is the persistence of the Hasse diagram: " << std::endl;
	pcoh_hd.output_diagram();
  
    return 0;
}

