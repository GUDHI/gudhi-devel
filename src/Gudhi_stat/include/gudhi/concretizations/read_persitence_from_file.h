/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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


#ifndef Read_Persitence_From_File_H
#define Read_Persitence_From_File_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>


namespace Gudhi
{
namespace Gudhi_stat
{


std::vector< std::pair< double , double > > read_standard_file( const char* filename )
{
	bool dbg = false;
	
	std::ifstream in;
    in.open( filename );
    if ( !( access( filename, F_OK ) != -1 ) )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}

    std::string line;
    std::vector< std::pair<double,double> > barcode;

    while (!in.eof())
    {
        getline(in,line);
        if ( !(line.length() == 0 || line[0] == '#') )
        {
            std::stringstream lineSS;
            lineSS << line;
            double beginn, endd;
            lineSS >> beginn;
            lineSS >> endd;
            if ( beginn > endd )
            {
                double b = beginn;
                beginn = endd;
                endd = b;
            }
            barcode.push_back( std::make_pair( beginn , endd ) );
            if (dbg)
            {
				std::cerr << beginn << " , " << endd << std::endl;
			}
        }
	}
	in.close();
	return barcode;
}//read_standard_file

std::vector< std::pair< double , double > > read_gudhi_file( const char* filename , size_t dimension = 0 )
{
	bool dbg = false;
	std::ifstream in;
    in.open( filename );

    std::string line;
    std::vector< std::pair<double,double> > barcode;

    while (!in.eof())
    {
        getline(in,line);
        if ( !(line.length() == 0 || line[0] == '#') )
        {
			if ( line.find("inf") != std::string::npos )
			{
				if ( dbg )
				{
					std::cerr << "This line: " << line << " contains infinite interval. We will skip it. \n";
				}
				continue;
			}
            std::stringstream lineSS;
            lineSS << line;
            double beginn, endd, field, dim;
            lineSS >> field;
            lineSS >> dim;
            lineSS >> beginn;
            lineSS >> endd;
            if ( beginn > endd )
            {
                double b = beginn;
                beginn = endd;
                endd = b;
            }
            if ( dim == dimension )
            { 
				barcode.push_back( std::make_pair( beginn , endd ) );
				if (dbg)
				{
					std::cerr << beginn << " , " << endd << std::endl;
				}
			}
        }
	}
	in.close();
	return barcode;
}//read_gudhi_file

}//namespace Gudhi_stat
}//namespace Gudhi




#endif
