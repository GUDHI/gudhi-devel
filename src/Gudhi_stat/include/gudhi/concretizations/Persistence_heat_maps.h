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
 
 
 /*
I would not say that homology is rovust when you speak about invariance with respect to continuous deformation. 
* 
Definition of simplicial compelx is redundant. There is no need to include vertices separatell.
* 
term 'codimension' used without being defined.
* 
"because boundary of boundary is always empty' -- in an intro survey one should not assume it.
* 
*what you call 'simple homotopies' in the literature is called 'free face collapses'. They are based on a simple homotopy theory, but this is not how they are called.
* 
* "hovewer when 1 \neq -1 ... -- not clear for people who do not know what you are talking about.
* 
* It is a urban legend that Cech complexes are difficult. They are not, one use exactly the same algorithm as in Rips complex construction and use minibal to check the intersection. 
* 
* correlations or measures of dissimilarity are typically not distances.
* 
* topological inference means soemthing else than you say in section 4 ?????????????????????????????????
* 
* Definition of persistence is wrong. Authors in the next paragraph are trying to provide a 'refined and ultimately more useful' definition which takes into account that classes can merge.
* This is a very clumzy and extremally confusing and in my opinion useless. What is the aim of providing wrong definition and then trying to fix it in informal text? Now it is not clear what do 
* they mean by birth and death later in the paper.
* In the previous iteration the feedback on this was given by both reviewers. Still, the authors failed to copy-pasted the correct definition 
* from vast literature. The fact that despite critical comments in the previous review none of many authors of the paper have found a problem in definition clearly indicate that authors are not 
* experts in persistent homology. In my opinion such a type of survey should be written by scientists who are more experiences in a field. 
* 
* 
  */


//standard include
#include <vector>
#include <sstream>
#include <iostream>

//gudhi include
#include "persistence_intervals.h"
#include "Distances_of_points_in_diagram.h"
#include "read_files_names.h"


using namespace std;

template <typename T>
class Persistence_heat_maps
{
public:
	Persistence_heat_maps(){};
    Persistence_heat_maps( Persistence_intervals<T>& interval , int number_of_pixels = 1000 , int diameter_of_gaussian = 4 , T min_ = -1 , T max_ = -1  );
    Persistence_heat_maps( char* name_of_file_with_names_of_files_with_intervals , int number_of_pixels = 1000 , int diameter_of_gaussian = 4 , T min_ = -1 , T max_ = -1  );

	
	void load_mean( const std::vector<Persistence_heat_maps*>& maps );
	void load_median( const std::vector<Persistence_heat_maps*>& maps );	
	void load_percentage_of_active( const std::vector<Persistence_heat_maps*>& maps );
    
    //put to file subroutine
    void write_to_file( const char* filename );
    
    //read from file subroutine
    Persistence_heat_maps( const char* filename );
    
    void plot( const char* filename );
private:
	std::vector< std::vector<T> > check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps );
    void construct( const Persistence_intervals<T>& intervals_ , T min_ = -1 , T max_ = -1);
    const Persistence_intervals<T>& intervals;
    int number_of_pixels;
    //int diameter_of_gaussian;
    //we should keep here the Gaussian filter
    
    //here we should have the scalling function
    
    //here we should have the info if we cut out everything below diagonal or not

    T min_;
    T max_;
    std::vector< std::vector<T> > heat_map;    
};


template <typename T>
void Persistence_heat_maps<T>::construct( std::vector< Persistence_intervals<T>* > intervals_ , T min_ , T max_ )
{
    bool dbg = false;

    for ( size_t i = 0 ; i != intervals_.size() ; ++i )
    {
        this->intervals.push_back( intervals_[i] );
    }

    if ( min_ == max_ )
    {
        //in this case, we want the program to set up the min_ and max_ values by itself.
        min_ = INT_MAX;
        max_ = -INT_MAX;
        for ( size_t i = 0 ; i != intervals_.size() ; ++i )
        {
            std::pair<T,T> min_max = intervals_[i]->min_max();
            if ( min_max.first < min_ )min_ = min_max.first;
            if ( min_max.second > max_ )max_ = min_max.second;
        }
        //now we have the structure filled in, and moreover we know min_ and max_ values of the interval, so we know the range.

        //add some more space:
        min_ -= fabs(max_ - min_)/100;
        max_ += fabs(max_ - min_)/100;
    }


    cerr << "min_ : " << min_ << endl;
    cerr << "max_ : " << max_ << endl;

    this->min_ = min_;
    this->max_ = max_;


    //initialization of the structure heat_map
    std::vector< std::vector< std::vector<T> > > heat_map_;
    for ( size_t i = 0 ; i != this->number_of_pixels ; ++i )
    {
        std::vector< std::vector<T> > row;
        for ( size_t j = 0 ; j != this->number_of_pixels ; ++j )
        {
            std::vector<T> v( this->intervals.size() );
            std::fill( v.begin() , v.end() , 0 );
            row.push_back( v );
        }
        heat_map_.push_back( row );
    }
    this->heat_map = heat_map_;

    if (dbg)cerr << "Done creating of the heat map, now we will fill in the structure \n";

    //and filling-in the structure:
    //for every persistence diagram:
    for ( size_t diag_no = 0 ; diag_no != this->intervals.size() ; ++diag_no )
    {
        //for every point in this diagram:
        for ( size_t point_no = 0 ; point_no != this->intervals[diag_no]->intervals.size() ; ++point_no )
        {
            //take the coordinates of this->intervals[diag_no][point_no]
            //this->heat_map[diag_no]
            size_t x_position_matrix = ((this->intervals[diag_no]->intervals[point_no].first-min_)/( max_ - min_ )) * this->number_of_pixels;
            size_t y_position_matrix = ((this->intervals[diag_no]->intervals[point_no].second-min_)/( max_ - min_ )) * this->number_of_pixels;

            if ( (x_position_matrix > this->number_of_pixels) || (y_position_matrix > this->number_of_pixels) )
            {
                cerr << "Points coords out of range : (" << this->intervals[diag_no]->intervals[point_no].first << " , " << this->intervals[diag_no]->intervals[point_no].second << "). This is probably due to wrong values of min and max being setup. This point will be disregarded \n";
                continue;
            }

            if ( dbg )
            {
                cerr << "The coords of the point are : (" << this->intervals[diag_no]->intervals[point_no].first << " , " << this->intervals[diag_no]->intervals[point_no].second << ") \n";
                cerr << "The matrix coords are : (" << x_position_matrix << " , " << y_position_matrix << ") \n";
                getchar();
            }

            if ( x_position_matrix == y_position_matrix )continue;

            size_t lower_side = (size_t)std::max( (int)(y_position_matrix-this->diameter_of_gaussian) , (int)0 );
            size_t upper_side = (size_t)std::min( (int)(y_position_matrix+this->diameter_of_gaussian) , (int)(this->number_of_pixels-1) );
            size_t left_side = (size_t)std::max(  (int)(x_position_matrix - this->diameter_of_gaussian) , (int)0 );
            size_t right_side = (size_t)std::min( (int)(x_position_matrix + this->diameter_of_gaussian) , (int)(this->number_of_pixels-1) );

            if ( dbg )
            {
                cerr << "lower_side : " << lower_side << endl;
                cerr << "upper_side : " << upper_side << endl;
                cerr << "left_side : " << left_side << endl;
                cerr <<  "right_side : " << right_side << endl;
            }

            for ( size_t y = lower_side ; y <= upper_side ; ++y )
            {
                for ( size_t x = left_side ; x <= right_side ; ++x )
                {
                    //cerr << "(" << x << "," << y << ")";
                    //make it more ambitious:
                    this->heat_map[x][y][diag_no]++;
                }
                //cerr << endl;
            }
            //getchar();
        }
    }
}//construct

template <typename T>
Persistence_heat_maps<T>::Persistence_heat_maps( std::vector< Persistence_intervals<T>* > intervals_ , int number_of_pixels_ , int diameter_of_gaussian_ , T min_ , T max_ ):number_of_pixels(number_of_pixels_),diameter_of_gaussian(diameter_of_gaussian_)
{
    this->construct( intervals_ , min_ , max_ );
}

template <typename T>
Persistence_heat_maps<T>::Persistence_heat_maps( char* name_of_file_with_names_of_files_with_intervals , int number_of_pixels_ , int diameter_of_gaussian_ , T min_ , T max_):number_of_pixels(number_of_pixels_),diameter_of_gaussian(diameter_of_gaussian_)
{
    std::vector< Persistence_intervals<T>* > intervals;
    std::vector<std::string> names = readFileNames( name_of_file_with_names_of_files_with_intervals );
    for ( size_t file_no = 0 ; file_no != names.size() ; ++file_no )
    {
        cout << "Reading file : " << names[file_no] << endl;
        Persistence_intervals<T>* interval = new Persistence_intervals<T>( (char*)names[file_no].c_str() );
        intervals.push_back( interval );
    }
    this->construct( intervals , min_ , max_);
}


double vector_median(std::vector<double> vec)
{
    if(vec.empty()) return 0;
    else
    {
        std::sort(vec.begin(), vec.end());
        if(vec.size() % 2 == 0)
        {
            return (vec[vec.size()/2 - 1] + vec[vec.size()/2]) / 2;
        }
        else
        {
            return vec[vec.size()/2];
        }
    }
}

template <typename T>
std::vector< std::vector<T> > Persistence_heat_maps<T>::check_and_initialize_maps( const std::vector<Persistence_heat_maps*>& maps )
{
	//checking if all the heat maps are of the same size:
	for ( size_t i = 0 ; i != maps.heat_maps.size() ; ++i )
    {
		if ( maps.heat_maps[i].size() != maps.heat_maps[0].size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
		}
		if ( maps.heat_maps[i][0].size() != maps.heat_maps[0][0].size() )
		{
			std::cerr << "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";
			throw "Sizes of Persistence_heat_maps are not compatible. The program will terminate now \n";			
		}
	} 	
	std::vector< std::vector<T> > heat_maps( maps.heat_maps.size() );
	for ( size_t i = 0 ; i != maps.heat_maps.size() ; ++i )
    {
		 std::vector<T> v( maps.heat_maps[0].size() , 0 );
		 heat_maps[i] = v;
	}
	return heat_maps;
}

template <typename T>
void Persistence_heat_maps<T>::load_median( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<T> > heat_maps = this->check_and_initialize_maps( maps );
	
	std::vector<T> to_compute_median( maps.size() );
    for ( size_t i = 0 ; i != heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_map[i].size() ; ++j )
        {
			for ( size_t aa = 0 ; aa != maps.size() ; ++aa )
			{
				to_compute_median[aa] = maps[i][j];
			}
            heat_maps[i][j] =  std::nth_element(to_compute_median.begin(), to_compute_median.begin() + to_compute_median.size()/2, to_compute_median.end());
        }       
    }
	this->heat_map = heat_maps;
}


template <typename T>
void Persistence_heat_maps<T>::load_mean( const std::vector<Persistence_heat_maps*>& maps )
{
	std::vector< std::vector<T> > heat_maps = this->check_and_initialize_maps( maps );
	
	std::vector<T> to_compute_median( maps.size() );
    for ( size_t i = 0 ; i != heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != heat_map[i].size() ; ++j )
        {
			T mean = 0;
			for ( size_t aa = 0 ; aa != maps.size() ; ++aa )
			{
				mean += maps[i][j];
			}
            heat_maps[i][j] =  mean/(T)maps.size();
        }       
    }
	this->heat_map = heat_maps;
}




template <typename T>
void Persistence_heat_maps<T>::load_percentage_of_active( const std::vector<Persistence_heat_maps*> maps , size_t cutoff )
{
	std::vector< std::vector<T> > heat_maps = this->check_and_initialize_maps( maps );

    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            int number_of_active_levels = 0;
            for ( size_t k = 0 ; k != this->heat_map[j][i].size() ; ++k )
            {
                if ( this->heat_map[j][i][k] ) number_of_active_levels++;
            }
            if ( number_of_active_levels > cutoff )
            {
                heat_maps[i][j] =number_of_active_levels;
            }
            else
            {
                heat_maps[i][j] = 0;
            }
        }       
    }

    this->heat_map = heat_maps;
}


template <typename T>
void Persistence_heat_maps<T>::plot( const char* filename )
{
	 ofstream out;
    out.open( filename );
	out << "plot      '-' matrix with image" << std::endl;
    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << endl;
    }

    out.close();
}


template <typename T>
void Persistence_heat_maps<T>::write_to_file( const char* filename )
{
	ofstream out;
	out.open( filename );
	for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            out << this->heat_map[i][j] << " ";
        }
        out << endl;
    }
	out.close();
}

template <typename T>
Persistence_heat_maps<T>::Persistence_heat_maps( const char* filename )
{
	bool dbg = true;
	
	ifstream in;
	in.open( filename );
	
	//checking if the file exist / if it was open. 
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}
	
	//now we read the file one by one. 
	std::string line;
	 
	while (!in.eof())
	{
        getline(in,line);
        std::stringstream lineSS;
        lineSS << line;
        
        std::vector< std::vector<T> > line_of_heat_map;    
        while ( lineSS.good() )
        {
			T point;
			lineSS >> point;
			line_of_heat_map.push_back( point );
			if ( dbg )
			{
				std::cout << point << " ";
			}
		}
		if ( dbg )
		{
			std::cout << std::endl;
		}

		this->heat_map.push_back( line_of_heat_map );
	}	
	in.close();
	if ( dbg )std::cout << "Done \n";
}
