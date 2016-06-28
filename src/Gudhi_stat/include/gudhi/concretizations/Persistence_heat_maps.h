#include <vector>
#include <sstream>
#include <iostream>


#include "persistence_intervals.h"
#include "Distances_of_points_in_diagram.h"
#include "read_files_names.h"


using namespace std;

template <typename T>
class Persistence_heat_maps
{
public:
    Persistence_heat_maps( std::vector< Persistence_intervals<T>* > intervals , int number_of_pixels = 1000 , int diameter_of_gaussian = 4 , T min_ = -1 , T max_ = -1  );
    Persistence_heat_maps( char* name_of_file_with_names_of_files_with_intervals , int number_of_pixels = 1000 , int diameter_of_gaussian = 4 , T min_ = -1 , T max_ = -1  );

    void output_median( char* filename );
    void output_mean( char* filename );
    void output_percentage_of_active( char* filename, size_t cutoff = 0);
private:
    void construct( std::vector< Persistence_intervals<T>* > intervals_ , T min_ = -1 , T max_ = -1);
    std::vector< Persistence_intervals<T>* > intervals;
    int number_of_pixels;
    int diameter_of_gaussian;

    T min_;
    T max_;

    std::vector< std::vector< std::vector<T> > > heat_map;
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
void Persistence_heat_maps<T>::output_median( char* filename )
{
    ofstream out;
    out.open( filename );

    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {

            out << vector_median( this->heat_map[j][i] ) << " ";
        }
        out << endl;
    }

    out.close();
}


template <typename T>
void Persistence_heat_maps<T>::output_mean( char* filename )
{
    ofstream out;
    out.open( filename );

    for ( size_t i = 0 ; i != this->heat_map.size() ; ++i )
    {
        for ( size_t j = 0 ; j != this->heat_map[i].size() ; ++j )
        {
            T mean = 0;
            for ( size_t k = 0 ; k != this->heat_map[j][i].size() ; ++k )
            {
                mean += this->heat_map[j][i][k];
            }
            mean /= (T)this->heat_map[j].size();
            out << mean << " ";
        }
        out << endl;
    }

    out.close();
}



template <typename T>
void Persistence_heat_maps<T>::output_percentage_of_active( char* filename , size_t cutoff )
{
    ofstream out;
    out.open( filename );

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
                out << number_of_active_levels << " ";
            }
            else
            {
                out << "0 ";
            }
        }
        out << endl;
    }

    out.close();
}
