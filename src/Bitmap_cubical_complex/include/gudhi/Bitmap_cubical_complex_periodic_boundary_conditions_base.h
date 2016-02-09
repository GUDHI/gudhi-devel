#pragma once
#include <cmath>
#include "Bitmap_cubical_complex_base.h"

using namespace std;

namespace Gudhi
{

namespace Cubical_complex
{

//in this class, we are storing all the elements which are in normal bitmap (i.e. the bitmap without the periodic boundary conditions). But, we set up the iterators and the procedures
//to compute boundary and coboundary in the way that it is all right. We assume here that all the cells that are on the left / bottom and so on remains, while all the cells on the
//right / top are not in the Bitmap_cubical_complex_periodic_boundary_conditions_base

template <typename T>
class Bitmap_cubical_complex_periodic_boundary_conditions_base : public Bitmap_cubical_complex_base<T>
{
public:
    //constructors that take an extra parameter:
    Bitmap_cubical_complex_periodic_boundary_conditions_base(){};
    Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> sizes , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed );
    Bitmap_cubical_complex_periodic_boundary_conditions_base( const char* perseusStyleFile );
    Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed );

    //overwritten methods co compute boundary and coboundary
    virtual std::vector< size_t > get_boundary_of_a_cell( size_t cell )const;
    std::vector< size_t > get_coboundary_of_a_cell( size_t cell )const;
    //inline unsigned get_dimension_of_a_cell( size_t cell )const;

protected:
    std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed;
    void set_up_containers( const std::vector<unsigned>& sizes )
    {

        unsigned multiplier = 1;
        for ( size_t i = 0 ; i != sizes.size() ; ++i )
        {
            this->sizes.push_back(sizes[i]);
            this->multipliers.push_back(multiplier);

            if ( directions_in_which_periodic_b_cond_are_to_be_imposed[i] )
            {
                multiplier *= 2*sizes[i];
            }
            else
            {
                multiplier *= 2*sizes[i]+1;
            }
        }
        //std::reverse( this->sizes.begin() , this->sizes.end() );
        this->data = std::vector<T>(multiplier,std::numeric_limits<T>::max());
        this->total_number_of_cells = multiplier;
    }
    Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> sizes );
    Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells );
    void construct_complex_based_on_top_dimensional_cells( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed );
};

template <typename T>
void Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::construct_complex_based_on_top_dimensional_cells( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed )
{
    this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;
    this->set_up_containers( dimensions );

    size_t i = 0;
    for ( typename Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Top_dimensional_cells_iterator it = this->top_dimensional_cells_begin() ; it != this->top_dimensional_cells_end() ; ++it )
    {
        *it = topDimensionalCells[i];
        ++i;
    }
    this->impose_lower_star_filtration();
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> sizes , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed )
{
    this->directions_in_which_periodic_b_cond_are_to_be_imposed = directions_in_which_periodic_b_cond_are_to_be_imposed;
    this->set_up_containers( sizes );
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base( const char* perseus_style_file )
{
    
    //for Perseus style files:
    bool dbg = false;
    ifstream inFiltration;
    inFiltration.open( perseus_style_file );
    unsigned dimensionOfData;
    inFiltration >> dimensionOfData;

    this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>( dimensionOfData , false );

    std::vector<unsigned> sizes;
    sizes.reserve( dimensionOfData );
    for ( size_t i = 0 ; i != dimensionOfData ; ++i )
    {
        int size_in_this_dimension;
        inFiltration >> size_in_this_dimension;
        if ( size_in_this_dimension < 0 )
        {
            this->directions_in_which_periodic_b_cond_are_to_be_imposed[i] = true;
        }
        sizes.push_back( abs(size_in_this_dimension) );
    }
    this->set_up_containers( sizes );

    typename Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Top_dimensional_cells_iterator it(*this);
    it = this->top_dimensional_cells_begin();

    while ( !inFiltration.eof() )
    {
        double filtrationLevel;
        inFiltration >> filtrationLevel;
        if ( inFiltration.eof() )break;

        if ( dbg )
        {
            cerr << "Cell of an index : "
                 << it.compute_index_in_bitmap()
                 << " and dimension: "
                 << this->get_dimension_of_a_cell(it.compute_index_in_bitmap())
                 << " get the value : " << filtrationLevel << endl;
        }
        *it = filtrationLevel;
        ++it;
    }
    inFiltration.close();
    this->impose_lower_star_filtration();
    
/*
     char* filename = (char*)perseus_style_file;
     //char* filename = "combustionWithPeriodicBoundaryConditions/v0/tV0_000000.float";
     ifstream file( filename , ios::binary | ios::ate );
     unsigned realSizeOfFile = file.tellg();
     file.close();
     realSizeOfFile = realSizeOfFile/sizeof(T);

     unsigned w, h, d;

     w = h = d = ceil(pow( realSizeOfFile , (double)(1/(double)3) ));

     T* slice = new T[w*h*d];
     if (slice == NULL)
     {
        cerr << "Allocation error, cannot allocate " << w*h*d*sizeof(T) << " bytes to store the data from the file. The program will now terminate \n";
        exit(EXIT_FAILURE);
     }

     FILE* fp;
     if ((fp=fopen( filename, "rb" )) == NULL )
     {
         cerr << "Cannot open the file: " << filename << ". The program will now terminate \n";
         exit(1);
     }

     clock_t read_begin = clock();
     fread( slice,4,w*h*d,fp );
     fclose(fp);
     cerr << "Time of reading the file : " << double(clock() - read_begin) / CLOCKS_PER_SEC << endl;


     clock_t begin_creation_bitap = clock();
     std::vector<T> data(slice,slice+w*h*d);
     delete[] slice;
     std::vector< unsigned > sizes;
     sizes.push_back(w);
     sizes.push_back(w);
     sizes.push_back(w);

    this->directions_in_which_periodic_b_cond_are_to_be_imposed.push_back( true );
    this->directions_in_which_periodic_b_cond_are_to_be_imposed.push_back( true );
    this->directions_in_which_periodic_b_cond_are_to_be_imposed.push_back( true );
    this->set_up_containers( sizes );

    size_t i = 0;
    for ( typename Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Top_dimensional_cells_iterator it = this->top_dimensional_cells_begin() ; it != this->top_dimensional_cells_end() ; ++it )
    {
        *it = data[i];
        ++i;
    }
    this->impose_lower_star_filtration();
    cerr << "Time of creation of a bitmap : " << double(clock() - begin_creation_bitap ) / CLOCKS_PER_SEC << endl;
*/
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> sizes )
{
    this->directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>( sizes.size() , false );
    this->set_up_containers( sizes );
}

template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells )
{
    std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = std::vector<bool>( dimensions.size() , false );
    this->construct_complex_based_on_top_dimensional_cells( dimensions , topDimensionalCells , directions_in_which_periodic_b_cond_are_to_be_imposed );
}





template <typename T>
Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::Bitmap_cubical_complex_periodic_boundary_conditions_base( std::vector<unsigned> dimensions , std::vector<T> topDimensionalCells , std::vector< bool > directions_in_which_periodic_b_cond_are_to_be_imposed )
{
    this->construct_complex_based_on_top_dimensional_cells( dimensions , topDimensionalCells , directions_in_which_periodic_b_cond_are_to_be_imposed );
}

//***********************Methods************************//

template <typename T>
std::vector< size_t >  Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_boundary_of_a_cell( size_t cell )const
{
    bool dbg = false;
    if ( dbg ){cerr << "Computations of boundary of a cell : " << cell << endl;}

    std::vector< size_t > boundary_elements;
    size_t cell1 = cell;
    for ( size_t i = this->multipliers.size() ; i != 0 ; --i )
    {
        unsigned position = cell1/this->multipliers[i-1];
        //this cell have a nonzero length in this direction, therefore we can compute its boundary in this direction.

        if ( position%2 == 1 )
        {
            //if there are no periodic boundary conditions in this direction, we do not have to do anything.
            if ( !directions_in_which_periodic_b_cond_are_to_be_imposed[i-1] )
            {
                //cerr << "A\n";
                boundary_elements.push_back( cell - this->multipliers[ i-1 ] );
                boundary_elements.push_back( cell + this->multipliers[ i-1 ] );
                if (dbg){cerr << cell - this->multipliers[ i-1 ] << " " << cell + this->multipliers[ i-1 ] << " ";}
            }
            else
            {
                //in this direction we have to do boundary conditions. Therefore, we need to check if we are not at the end.
                if ( position != 2*this->sizes[ i-1 ]-1 )
                {
                    //cerr << "B\n";
                    boundary_elements.push_back( cell - this->multipliers[ i-1 ] );
                    boundary_elements.push_back( cell + this->multipliers[ i-1 ] );
                    if (dbg){cerr << cell - this->multipliers[ i-1 ] << " " << cell + this->multipliers[ i-1 ] << " ";}
                }
                else
                {
                    //cerr << "C\n";
                    boundary_elements.push_back( cell - this->multipliers[ i-1 ] );
                    boundary_elements.push_back( cell - (2*this->sizes[ i-1 ]-1)*this->multipliers[ i-1 ] );
                    if (dbg){cerr << cell - this->multipliers[ i-1 ] << " " << cell - (2*this->sizes[ i-1 ]-1)*this->multipliers[ i-1 ] << " ";}
                }
            }
        }
        cell1 = cell1%this->multipliers[i-1];
    }
    return boundary_elements;
}

template <typename T>
std::vector< size_t > Bitmap_cubical_complex_periodic_boundary_conditions_base<T>::get_coboundary_of_a_cell( size_t cell )const
{
    std::vector<unsigned> counter = this->compute_counter_for_given_cell( cell );
    std::vector< size_t > coboundary_elements;
    size_t cell1 = cell;
    for ( size_t i = this->multipliers.size() ; i != 0 ; --i )
    {
        unsigned position = cell1/this->multipliers[i-1];
        //if the cell has zero length in this direction, then it will have cbd in this direction.
        if ( position%2 == 0 )
        {
            if ( !this->directions_in_which_periodic_b_cond_are_to_be_imposed[i-1] )
            {
                //no periodic boundary conditions in this direction
                if ( (counter[i-1] != 0) && (cell > this->multipliers[i-1]) )
                {
                    coboundary_elements.push_back( cell - this->multipliers[i-1] );
                }
                if ( (counter[i-1] != 2*this->sizes[i-1]) && (cell + this->multipliers[i-1] < this->data.size()) )
                {
                    coboundary_elements.push_back( cell + this->multipliers[i-1] );
                }
            }
            else
            {
                //we want to have periodic boundary conditions in this direction
                if ( counter[i-1] != 0 )
                {
                    coboundary_elements.push_back( cell - this->multipliers[i-1] );
                    coboundary_elements.push_back( cell + this->multipliers[i-1] );
                }
                else
                {
                    //in this case counter[i-1] == 0.
                    coboundary_elements.push_back( cell + this->multipliers[i-1] );
                    coboundary_elements.push_back( cell + (2*this->sizes[ i-1 ]-1)*this->multipliers[i-1] );
                }
            }
        }

        cell1 = cell1%this->multipliers[i-1];
    }
    return coboundary_elements;
}



}//Cubical_complex
}//namespace Gudhi
