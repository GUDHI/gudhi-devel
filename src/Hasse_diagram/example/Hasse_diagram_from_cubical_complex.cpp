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
#include <gudhi/Bitmap_cubical_complex.h>

// standard stuff
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

int main() 
{    
	//This is an example of construction of Hasse diagram from a cubical complex. In the example
	//below we first define 3 by 3 by 3 cubical complex and assign a filtration on it. Later
	//we convert it to a standard Hasse_diagram and display it. At the end, we will convert it
	//to an object of a class Hasse_diagram_persistence and compute persitence of it.
	
    typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
    typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
    std::vector< unsigned > sizes(3);
    sizes[0] = sizes[1] = sizes[2] = 3;
    //this is 3 by 3 by 3 cubical complex representing a cubical sphere embedded on two dimensional torus. 
    //The lifespan of a sphere is from 1 to 10.
    
    std::vector<double> top_dimensional_cells_data = 
    {
		1,1,1,
		1,1,1,
		1,1,1,
//
		1,1,1,
		1,10,1,	
		1,1,1,
//
		1,1,1,
		1,1,1,
		1,1,1
	};
	Bitmap_cubical_complex b(sizes,top_dimensional_cells_data);
	
	typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int,double,double> Cell;	
    typedef Gudhi::Hasse_diagram::Hasse_diagram_persistence<Cell> Hasse_diag;
    Hasse_diag* hd = Gudhi::Hasse_diagram::convert_to_Hasse_diagram_persistence<Bitmap_cubical_complex,Cell>( b );
    
    std::cout << "Here is the Hasse diagram obtained from the cunical complex : " << *hd << std::endl;
    
   
	typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
	typedef Gudhi::persistent_cohomology::Persistent_cohomology<Hasse_diag, Field_Zp> Persistent_cohomology;

	Persistent_cohomology pcoh(*hd,true);  
	int field_characteristic = 11;
	double min_persistence = 0;
	
	pcoh.init_coefficients(field_characteristic);    
	pcoh.compute_persistent_cohomology(min_persistence);
	
	std::cout << "And here is the persistent homology of the Hasse diagram : " << std::endl;
	pcoh.output_diagram();
  
    return 0;
}

