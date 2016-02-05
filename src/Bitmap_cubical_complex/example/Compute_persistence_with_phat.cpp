
 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       Pawel Dlotko
  *
  *    Copyright (C) 2015  INRIA Saclay (France)
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
#include <gudhi/Bitmap_cubical_complex_base.h>
#include <gudhi/Bitmap_cubical_complex_periodic_boundary_conditions_base.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Compute_persistence_with_phat.h>

using namespace Gudhi;
using namespace Gudhi::Cubical_complex;

//standard stuff
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;


int main( int argc , char** argv )
{
    if ( argc != 2 )
    {
        cout << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap at the input. The program will now terminate.\n";
        return 1;
    }

    Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > b( argv[1] );
    cerr << "Cubical complex created \n";


    Compute_persistence_with_phat< Bitmap_cubical_complex< Bitmap_cubical_complex_base<double> > , double > phat(&b);
    phat::persistence_pairs pairs = phat.compute_persistence_pairs_dualized_chunk_reduction();
    std::pair< std::vector< std::vector<double> > , std::vector< std::vector< std::pair<double,double> > > > persistence = phat.get_the_intervals( pairs );
    writeBettiNumbersAndPersistenceIntervalsToFile( "phat_persistence" , persistence );

    return 0;
}
