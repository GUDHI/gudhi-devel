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


#ifndef Read_Persitence_From_File_H
#define Read_Persitence_From_File_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unistd.h>


namespace Gudhi
{
namespace Gudhi_stat
{
			
	
/**
 * This procedure reads birth-death dagta from a file. We assume that in the file, there may be one type of string 'inf' or 'Inf'. If the second parameter of the program is set to -1,
 * then those vales are ignored. If the second parameter of this program is set to a positive value, then the infinite intervals will be substituted by that number.
**/ 
std::vector< std::pair< double,double > > read_persistence_file_that_may_contain_inf_string( char* filename , double what_to_substitute_for_infinite_bar = -1 )
{
	
	bool dbg = true;
	
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read do not exist. The program will now terminate \n";
	}
	
	std::string line;
    std::vector< std::pair<double,double> > barcode;

	std::ifstream in;
    in.open( filename );
    while (!in.eof())
    {
        getline(in,line);        
        if ( !(line.length() == 0 || line[0] == '#') )
        {
			std::stringstream lineSS(line);           
            double beginn, endd;
			if ( (line.find("inf") != std::string::npos) || (line.find("Inf") != std::string::npos) )
			{
				if ( dbg )
				{
					std::cerr << "This line: " << line << " contains infinite interval. We will skip it. \n";
				}
				if ( what_to_substitute_for_infinite_bar != -1 )
				{
					lineSS >> beginn;
					endd = what_to_substitute_for_infinite_bar;
				}
				else
				{
					continue;
				}
			}
			else
			{           				       
				lineSS >> beginn;
				lineSS >> endd;
			}
            if ( beginn > endd )
            {
                std::swap(beginn,endd);
            }  
            if ( dbg )
			{
				std::cerr << "Getting an interval : " << beginn << "," << endd << std::endl;
			}          
			barcode.push_back( std::make_pair( beginn , endd ) );			
        }
	}
	return barcode;    
}//readFileNames	
	
	
/**
 * This procedure reads names of files which are stored in a file. 
**/ 
std::vector< std::string > readFileNames( const char* filenameWithFilenames )
{
    bool dbg = false;
    
    if ( !( access( filenameWithFilenames, F_OK ) != -1 ) )
	{
		std::cerr << "The file : " << filenameWithFilenames << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read do not exist. The program will now terminate \n";
	}

    std::vector< std::string > result;
    std::ifstream in;
    in.open( filenameWithFilenames );
    std::string line;
    while (!in.eof())
    {
        getline(in,line);
        line.erase( std::remove_if( line.begin(), line.end(), ::isspace) , line.end() );

        if (dbg){std::cerr << "line : " << line << std::endl;}

        if ( (line.length() == 0) || (line[0] == '#') )
        {
            //in this case we have a file name. First we should remove all the white spaces.
            if ( dbg ){std::cerr << "This is a line with comment, it will be ignored n";}
        }
        else
        {
            result.push_back( line.c_str() );
            if (dbg){std::cerr << "Line after removing white spaces : " << line << std::endl;}
        }
	}
    in.close();

    return result;
}//readFileNames


/**
 * This method reads persistence from standalone file. The format of the file is as follows:
 * In every line there are two numbers which denotes birth and death of a persistence interval. The file is not supposed to contains any nonwhite characters excet from digits. In particular symbols 
 *  like 'inf' are not allowed.
 * If begin of the interval is greater than the end of the interval, those two numbers are swapped.
**/ 
std::vector< std::pair< double , double > > read_standard_persistence_file( const char* filename )
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
            std::stringstream lineSS(line);            
            double beginn, endd;
            lineSS >> beginn;
            lineSS >> endd;
            if ( beginn == endd )continue;
            if ( beginn > endd )
            {
                std::swap(beginn,endd);                
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



/**
 * This procedure reads Gudhi style file. The format of a Gudhi style file is as follows:
 * Each line consist of the following information:
 * A prime number p indicating that the class is nontrivial over Z_p field.
 * A positive numebr, being a dimension of a persistence class
 * A birth and a death time of a class.
 * Death time may be infitnity, in which case a string 'inf' is used.
 * If begin of the interval is greater than the end of the interval, those two numbers are swapped.
 * Note that this procedure reads persistence in a single dimension. The dimension of intervals that
 * are to be read are determined by the second parameter of the function.
**/ 
std::vector< std::pair< double , double > > read_gudhi_persistence_file_in_one_dimension( const char* filename , size_t dimension = 0 , double what_to_substitute_for_infinite_bar = -1)
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
				if ( what_to_substitute_for_infinite_bar != -1 )
				{
					double beginn, field, dim;
					std::stringstream lineSS(line);			
					lineSS >> field;
					lineSS >> dim;
					lineSS >> beginn;
					if ( dim == dimension )
					{
						if ( beginn > what_to_substitute_for_infinite_bar )
						{
							barcode.push_back( std::make_pair( what_to_substitute_for_infinite_bar , beginn ) );
						}
						else
						{
							barcode.push_back( std::make_pair( beginn , what_to_substitute_for_infinite_bar ) );
						}
						if (dbg)
						{
							std::cerr << beginn << " , " << what_to_substitute_for_infinite_bar << std::endl;
						}
					} 
				}
				continue;
			}
			std::stringstream lineSS(line);			
			double beginn, endd, field, dim;
			lineSS >> field;
			lineSS >> dim;
			lineSS >> beginn;
			lineSS >> endd;
			if ( beginn > endd )
			{
				std::swap(beginn,endd);				
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


std::vector< std::vector< double > > read_numbers_from_file_line_by_line( const char* filename )
{
	bool dbg = false;
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}

	std::vector< std::vector< double > > result;
	double number;

	std::ifstream in(filename);
	std::string line;
	while ( in.good() )
	{
		std::getline(in,line);
		std::stringstream ss(line);		

		if ( dbg )std::cerr << "\n Reading line : " << line << std::endl;

		std::vector< double > this_line;
		while ( ss.good() )
		{
			ss >> number;
			this_line.push_back( number );
			if ( dbg )std::cerr << number << " ";
		}
		if ( this_line.size() && in.good() ) result.push_back( this_line );
	}		
	in.close();

	return result;
}//read_numbers_from_file_line_by_line

}//namespace Gudhi_stat
}//namespace Gudhi




#endif
	
