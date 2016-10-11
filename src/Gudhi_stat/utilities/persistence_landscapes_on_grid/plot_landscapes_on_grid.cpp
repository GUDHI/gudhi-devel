/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA (France)
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



#include <gudhi/abstract_classes/Abs_Topological_data.h>
#include <gudhi/concretizations/Persistence_landscape_on_grid.h>



using namespace Gudhi;
using namespace Gudhi::Gudhi_stat;

#include <iostream>
#include <sstream>


int main( int argc , char** argv )
{
	std::cout << "This program plot persistence landscape on grid stored in a file (the file needs to be created beforehand). Please call the code with the name of a landsape on grid file \n";	
	if ( argc == 1 )
	{
		std::cout << "Wrong parameters of a program call, the program will now terminate \n";
		return 1;
	}
	Persistence_landscape_on_grid l;
	l.load_landscape_from_file( argv[1] );
	
	std::stringstream ss;
	ss << argv[1] << "_gnuplot_script";
	l.plot( ss.str().c_str() );
	
	std::cout << "Done \n";
		
	return 0;
}
