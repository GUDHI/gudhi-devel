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



namespace Gudhi
{
namespace Persistence_representations
{
			
	
	
/**
 * This procedure reads names of files which are stored in a file. 
**/ 
std::vector< std::string > readFileNames( const char* filenameWithFilenames )
{
    bool dbg = false;
    
    std::ifstream in(filenameWithFilenames);
	if ( !in.good() )
	{
		std::cerr << "The file : " << filenameWithFilenames << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read do not exist. The program will now terminate \n";
	}

    std::vector< std::string > result;    
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
            result.push_back( line );
            if (dbg){std::cerr << "Line after removing white spaces : " << line << std::endl;}
        }
	}
    in.close();

    return result;
}//readFileNames



std::vector< std::vector< double > > read_numbers_from_file_line_by_line( const char* filename )
{
	bool dbg = false;
	std::ifstream in(filename);
	if ( !in.good() )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}

	std::vector< std::vector< double > > result;
	double number;

	
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


/**
 * Universal procedure to read files with persistence. It ignores the lines starting from # (treat them as comments). 
 * It reads the fist line which is not a comment and assume that there are some numerical entries over there. The program assume
 * that each other line in the file, which is not a comment, have the same number of numerical entries (2, 3 or 4). 
 * If there are two numerical entries per line, then the function assume that they are birth/death coordinates. 
 * If there are three numerical entries per line, then the function assume that they are: dimension and birth/death coordinates. 
 * If there are four numerical entries per line, then the function assume that they are: thc characteristic of a filed over which 
 * persistence was computed, dimension and birth/death coordinates. 
 * The 'inf' string can appear only as a last element of a line. 
 * The procedure returns vector of persistence pairs. 
**/ 
std::vector<std::pair<double,double> > read_persistence_intervals_in_one_dimension_from_file(std::string const& filename, int dimension=-1 , double what_to_substitute_for_infinite_bar = -1 )
{
	bool dbg = false;
	std::ifstream in;
	in.open( filename );
	//checking if the file exist:
	if ( !in.good() )
	{
		std::cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}
	
	
	

	std::string line;
	std::vector< std::pair<double,double> > barcode;
	
	int number_of_entries_per_line = -1;

	while (!in.eof())
	{
		getline(in,line);
		if ( dbg )std::cerr << "Reading line : " << line << std::endl;
		if ( !(line.length() == 0 || line[0] == '#') )
		{
			//If we do not know how many entries per line we have, we check it in below.
			if ( number_of_entries_per_line == -1 )
			{	
				number_of_entries_per_line = 0;
				std::string line_copy(line);			
				if ( line_copy.find("inf") != std::string::npos )
				{					
					size_t np = line_copy.find("inf");
					//replace symbols 'inf' in line_copy with whilespaces:					
					line_copy[np] = ' ';
					line_copy[np+1] = ' ';
					line_copy[np+2] = ' ';		
					number_of_entries_per_line = 1;			
				}				
				//check how many entries we have in the line.				
				std::stringstream ss( line_copy );
				double number;
				std::vector<double> this_line;
				while ( ss >> number )
				{
					this_line.push_back( number );
				}
				number_of_entries_per_line += (int)this_line.size();						
				if ( dbg )
				{
					std::cerr << "number_of_entries_per_line : " << number_of_entries_per_line << ". This number was obtained by analyzing this line : " << line << std::endl;
				}
				if ( (number_of_entries_per_line < 2) || ( number_of_entries_per_line > 4 ) )
				{
					std::cerr << "The input file you have provided have wrong number of numerical entries per line. The program will now terminate. \n";
					throw "The input file you have provided have wrong number of numerical entries per line. The program will now terminate. \n";
				}
			}
			//In case there is an 'inf' string in this line, we are dealing with this situation in below.
			if ( line.find("inf") != std::string::npos )
			{
				if ( dbg )
				{
					std::cerr << "This line: " << line << " contains infinite interval. \n";
				}
				//first we substitute inf by whitespaces:
				size_t np = line.find("inf");									
				line[np] = ' ';
				line[np+1] = ' ';
				line[np+2] = ' ';	
				if ( what_to_substitute_for_infinite_bar != -1 )
				{
					double beginn, field, dim;
					std::stringstream lineSS(line);			
					if ( number_of_entries_per_line == 4 )lineSS >> field;
					if ( number_of_entries_per_line >= 3 )
					{
						 lineSS >> dim;
					}
					else
					{
						dim = dimension;
					}
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
							std::cerr << "this is the line that is going to the output : " << beginn << " , " << what_to_substitute_for_infinite_bar << std::endl;
						}
					} 
				}
				else
				{
					//this is a line with infinity. Since the variable what_to_substitute_for_infinite_bar have not been set up, it means that this line will be skipped. 
					if ( dbg )
					{
						std::cerr << "We will skip it \n";
					}
				}
				continue;
			}
			else
			{
				//Then, we read the content of the line. We know that it do not contain 'inf' substring.
				std::stringstream lineSS(line);			
				double beginn, endd, field, dim;
				if ( number_of_entries_per_line == 4 )lineSS >> field;
				if ( number_of_entries_per_line >= 3 )
				{
					 lineSS >> dim;
				}
				else
				{
					dim = dimension;
				}
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
						std::cerr << "This is a line that is going to the output : " <<  beginn << " , " << endd << std::endl;
					}
				}
				else
				{
					if ( (number_of_entries_per_line==3) && (dimension == -1) )
					{
						barcode.push_back( std::make_pair( beginn , endd ) );
					}
				}
			}
		}
		else
		{
			if ( dbg )
			{
				std::cerr << "This is a comment line \n";
			}
		}
	}
	in.close();
	if ( dbg )std::cerr << "End of reading \n";  
	
	return barcode;
}//read_persistence_intervals_in_one_dimension_from_file

}//namespace Gudhi_stat
}//namespace Gudhi




#endif
	
