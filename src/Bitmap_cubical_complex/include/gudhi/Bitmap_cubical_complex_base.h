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

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <limits>
#include "counter.h"


using namespace std;

namespace Gudhi
{

namespace Cubical_complex
{



/**
 * This is a class implementing a basic bitmap data structure to store cubical complexes. 
 * It implements only the most basic subroutines.
 * The idea of the bitmap is the following. Our aim is to have a memory efficient 
 * data structure to store d-dimensional cubical complex 
 * C being a cubical decomposition
 * of a rectangular region of a space. This is achieved by storing C as a 
 * vector of bits (this is where the name 'bitmap' came from). 
 * Each cell is represented by a single
 * bit (in case of black and white bitmaps, or by a single element of a type T 
 * (here T is a filtration type of a bitmap, typically a double). 
 * All the informations needed for homology and
 * persistent homology computations (like dimension of a cell, boundary and 
 * coboundary elements of a cell, are then obtained from the 
 * position of the element in C.
 * The default filtration used in this implementation is the lower star filtration.
 */
template <typename T>
class Bitmap_cubical_complex_base
{
public:
    /**
    * There are a few constructors of a Bitmap_cubical_complex_base class. 
    * First one, that takes vector<unsigned>, creates an empty bitmap of a dimension equal 
    * the number of elements in the
    * input vector and size in the i-th dimension equal the number in the position i-of the input vector.
    */
    Bitmap_cubical_complex_base( std::vector<unsigned>& sizes );
    /**
    * The second constructor takes as a input a Perseus style file. For more details, 
    * please consult the documentations of 
    * Perseus software as well as examples attached to this
    * implementation.
    **/
    Bitmap_cubical_complex_base( const char* perseus_style_file );
    /**
    * The last constructor of a Bitmap_cubical_complex_base class accepts vector of dimensions (as the first one) 
    * together with vector of filtration values of top dimensional cells.
    **/
    Bitmap_cubical_complex_base( std::vector<unsigned>& dimensions , const std::vector<T>& top_dimensional_cells );

    /**
    * The functions get_boundary_of_a_cell, get_coboundary_of_a_cell, get_dimension_of_a_cell 
    * and get_cell_data are the basic 
     * functions that compute boundary / coboundary / dimension and the filtration
     * value form a position of a cell in the structure of a bitmap. The input parameter of all of those function is a 
     * non-negative integer, indicating a position of a cube in the data structure.
     * In the case of functions that compute (co)boundary, the output is a vector if non-negative integers pointing to 
     * the positions of (co)boundary element of the input cell.
    */
    inline std::vector< size_t > get_boundary_of_a_cell( size_t cell )const;
    /**
    * The functions get_coboundary_of_a_cell, get_coboundary_of_a_cell, 
    * get_dimension_of_a_cell and get_cell_data are the basic 
    * functions that compute boundary / coboundary / dimension and the filtration
    * value form a position of a cell in the structure of a bitmap. 
    * The input parameter of all of those function is a non-negative integer, 
    * indicating a position of a cube in the data structure.
    * In the case of functions that compute (co)boundary, the output is a vector if 
    * non-negative integers pointing to the 
    * positions of (co)boundary element of the input cell.
    **/
    inline std::vector< size_t > get_coboundary_of_a_cell( size_t cell )const;
    /**
    * In the case of get_dimension_of_a_cell function, the output is a non-negative integer 
    * indicating the dimension of a cell.
    **/
    inline unsigned get_dimension_of_a_cell( size_t cell )const;
    /**
    * In the case of get_cell_data, the output parameter is a reference to the value of a cube in a given position.
    **/
    inline T& get_cell_data( size_t cell );


    /**
    * Typical input used to construct a baseBitmap class is a filtration given at the top dimensional cells. 
    * Then, there are a few ways one can pick the filtration of lower dimensional
    * cells. The most typical one is by so called lower star filtration. This function is always called by any 
    * constructor which takes the top dimensional cells. If you use such a constructor,
    * then there is no need to call this function. Call it only if you are putting the filtration 
    * of the cells by your own (for instance by using Top_dimensional_cells_iterator).
    **/
    void impose_lower_star_filtration();//assume that top dimensional cells are already set.

    /**
    * Returns dimension of a complex.
    **/
    inline unsigned dimension()const{ return sizes.size(); }

    /**
    * Returns number of all cubes in the data structure.
    **/
    inline unsigned size_of_bitmap()const
    {
        return this->data.size();
    }

    /**
    * Writing to stream operator.
    **/
    template <typename K>
    friend ostream& operator << ( ostream & os , const Bitmap_cubical_complex_base<K>& b );

    //ITERATORS

    /**
    * Iterator through all cells in the complex (in order they appear in the structure -- i.e. 
    * in lexicographical order).
    **/
    typedef typename std::vector< T >::iterator all_cells_iterator;
    all_cells_iterator all_cells_begin()const
    {
        return this->data.begin();
    }
    all_cells_iterator all_cells_end()const
    {
        return this->data.end();
    }


    typedef typename std::vector< T >::const_iterator all_cells_const_iterator;
    all_cells_const_iterator all_cells_const_begin()const
    {
        return this->data.begin();
    }
    all_cells_const_iterator all_cells_const_end()const
    {
        return this->data.end();
    }

    /**
    * Iterator through top dimensional cells of the complex. The cells appear in order they are stored 
    * in the structure (i.e. in lexicographical order)
    **/
    class Top_dimensional_cells_iterator : std::iterator< std::input_iterator_tag, double >
    {
        public:
            Top_dimensional_cells_iterator( Bitmap_cubical_complex_base& b ):b(b)
            {
                for ( size_t i = 0 ; i != b.dimension() ; ++i )
                {
                    this->counter.push_back(0);
                }
            }
            Top_dimensional_cells_iterator operator++()
            {
                //first find first element of the counter that can be increased:
                size_t dim = 0;
                while ( ( dim != this->b.dimension() ) && ( this->counter[dim] == this->b.sizes[dim]-1 ) )++dim;

                if ( dim != this->b.dimension() )
                {
                    ++this->counter[dim];
                    for ( size_t i = 0 ; i != dim ; ++i )
                    {
                        this->counter[i] = 0;
                    }
                }
                else
                {
                    ++this->counter[0];
                }
                return *this;
            }
            Top_dimensional_cells_iterator operator++(int)
            {
                Top_dimensional_cells_iterator result = *this;
                ++(*this);
                return result;
            }
            Top_dimensional_cells_iterator operator =( const Top_dimensional_cells_iterator& rhs )
            {
                this->counter = rhs.counter;
                this->b = rhs.b;
                return *this;
            }
            bool operator == ( const Top_dimensional_cells_iterator& rhs )const
            {
                if ( &this->b != &rhs.b )return false;
                if ( this->counter.size() != rhs.counter.size() )return false;
                for ( size_t i = 0 ; i != this->counter.size() ; ++i )
                {
                    if ( this->counter[i] != rhs.counter[i] )return false;
                }
                return true;
            }
            bool operator != ( const Top_dimensional_cells_iterator& rhs )const
            {
                return !(*this == rhs);
            }

            T& operator*()
            {
                //given the counter, compute the index in the array and return this element.
                unsigned index = 0;
                for ( size_t i = 0 ; i != this->counter.size() ; ++i )
                {
                    index += (2*this->counter[i]+1)*this->b.multipliers[i];
                }
                return this->b.data[index];
            }

            size_t compute_index_in_bitmap()const
            {
                size_t index = 0;
                for ( size_t i = 0 ; i != this->counter.size() ; ++i )
                {
                    index += (2*this->counter[i]+1)*this->b.multipliers[i];
                }
                return index;
            }

            void print_counter()const
            {
                for ( size_t i = 0 ; i != this->counter.size() ; ++i )
                {
                    cout << this->counter[i] << " ";
                }
            }
            friend class Bitmap_cubical_complex_base;
        protected:
            std::vector< unsigned > counter;
            Bitmap_cubical_complex_base& b;
    };
    Top_dimensional_cells_iterator top_dimensional_cells_begin()
    {
        Top_dimensional_cells_iterator a(*this);
        return a;
    }
    Top_dimensional_cells_iterator top_dimensional_cells_end()
    {
        Top_dimensional_cells_iterator a(*this);
        for ( size_t i = 0 ; i != this->dimension() ; ++i )
        {
            a.counter[i] = this->sizes[i]-1;
        }
        a.counter[0]++;
        return a;
    }


//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//


//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//
//****************************************************************************************************************//

protected:
    std::vector<unsigned> sizes;
    std::vector<unsigned> multipliers;
    std::vector<T> data;
    size_t total_number_of_cells;
    void set_up_containers( std::vector<unsigned>& sizes )
    {
        unsigned multiplier = 1;
        for ( size_t i = 0 ; i != sizes.size() ; ++i )
        {
            this->sizes.push_back(sizes[i]);
            this->multipliers.push_back(multiplier);
            //multiplier *= 2*(sizes[i]+1)+1;
            multiplier *= 2*sizes[i]+1;
        }
        //std::reverse( this->sizes.begin() , this->sizes.end() );
        std::vector<T> data(multiplier);
        std::fill( data.begin() , data.end() , std::numeric_limits<int>::max() );
        this->total_number_of_cells = multiplier;
        this->data = data;
    }

    size_t compute_position_in_bitmap( std::vector< unsigned >& counter )
    {
        size_t position = 0;
        for ( size_t i = 0 ; i != this->multipliers.size() ; ++i )
        {
            position += this->multipliers[i]*counter[i];
        }
        return position;
    }

     std::vector<unsigned> compute_counter_for_given_cell( size_t cell )const
    {
        std::vector<unsigned> counter;
        for ( size_t dim = this->sizes.size() ; dim != 0 ; --dim )
        {
            counter.push_back(cell/this->multipliers[dim-1]);
            cell = cell%this->multipliers[dim-1];
        }
        std::reverse( counter.begin() , counter.end() );
        return counter;
    }

    std::vector< size_t > 
    generate_vector_of_shifts_for_bitmaps_with_periodic_boundary_conditions
    ( std::vector< bool >& directions_for_periodic_b_cond );
};




template <typename K>
ostream& operator << ( ostream & out , const Bitmap_cubical_complex_base<K>& b )
{
    for ( typename Bitmap_cubical_complex_base<K>::all_cells_const_iterator 
          it = b.all_cells_const_begin() ; it != b.all_cells_const_end() ; ++it )
    {
        out << *it << " ";
    }
    return out;
}


template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base
( std::vector<unsigned>& sizes )
{
    this->set_up_containers( sizes );
}

template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base
( std::vector<unsigned>& sizes_in_following_directions , const std::vector<T>& top_dimensional_cells )
{
    this->set_up_containers( sizes_in_following_directions );

    size_t number_of_top_dimensional_elements = 1;
    for ( size_t i = 0 ; i != sizes_in_following_directions.size() ; ++i )
    {
        number_of_top_dimensional_elements *= sizes_in_following_directions[i];
    }
    if ( number_of_top_dimensional_elements != top_dimensional_cells.size() )
    {
        cerr << 
       "Error in constructor\
        Bitmap_cubical_complex_base\
       ( std::vector<size_t> sizes_in_following_directions , std::vector<float> top_dimensional_cells ).\
       Number of top dimensional elements that follow from sizes_in_following_directions vector is different\
       than the size of top_dimensional_cells vector." << endl;
        throw("Error in constructor Bitmap_cubical_complex_base( std::vector<size_t> sizes_in_following_directions,\
        std::vector<float> top_dimensional_cells )\
      . Number of top dimensional elements that follow from sizes_in_following_directions vector is different than the\
       size of top_dimensional_cells vector.");
    }

    Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
    size_t index = 0;
    for ( it = this->top_dimensional_cells_begin() ; it != this->top_dimensional_cells_end() ; ++it )
    {
        (*it) = top_dimensional_cells[index];
        ++index;
    }
    this->impose_lower_star_filtration();
}


template <typename T>
Bitmap_cubical_complex_base<T>::Bitmap_cubical_complex_base( const char* perseus_style_file )
{
    bool dbg = false;
    ifstream inFiltration, inIds;
    inFiltration.open( perseus_style_file );
    unsigned dimensionOfData;
    inFiltration >> dimensionOfData;

    if (dbg){cerr << "dimensionOfData : " << dimensionOfData << endl;}

    std::vector<unsigned> sizes;
    for ( size_t i = 0 ; i != dimensionOfData ; ++i )
    {
        unsigned size_in_this_dimension;
        inFiltration >> size_in_this_dimension;
        size_in_this_dimension = abs(size_in_this_dimension);
        sizes.push_back( size_in_this_dimension );
        if (dbg){cerr << "size_in_this_dimension : " << size_in_this_dimension << endl;}
    }
    this->set_up_containers( sizes );

    Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
    it = this->top_dimensional_cells_begin();

    //TODO -- over here we also need to read id's of cell and put them to bitmapElement structure!
    while ( !inFiltration.eof() )
    {
        double filtrationLevel;
        inFiltration >> filtrationLevel;
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
}


template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::get_boundary_of_a_cell( size_t cell )const
{
    bool bdg = false;
    //first of all, we need to take the list of coordinates in which the cell has nonzero length. 
    //We do it by using modified version to compute dimension of a cell:
    std::vector< unsigned > dimensions_in_which_cell_has_nonzero_length;
    unsigned dimension = 0;
    size_t cell1 = cell;
    for ( size_t i = this->multipliers.size() ; i != 0 ; --i )
    {
        unsigned position = cell1/multipliers[i-1];
        if ( position%2 == 1 )
        {
            dimensions_in_which_cell_has_nonzero_length.push_back(i-1);
            dimension++;
        }
        cell1 = cell1%multipliers[i-1];
    }

    if (bdg)
    {
        cerr << "dimensions_in_which_cell_has_nonzero_length : \n";
        for ( size_t i = 0 ; i != dimensions_in_which_cell_has_nonzero_length.size() ; ++i )
        {
            cerr << dimensions_in_which_cell_has_nonzero_length[i] << endl;
        }
        getchar();
    }

    std::vector< size_t > boundary_elements;
    if ( dimensions_in_which_cell_has_nonzero_length.size() == 0 )return boundary_elements;
    for ( size_t i = 0 ; i != dimensions_in_which_cell_has_nonzero_length.size() ; ++i )
    {
        boundary_elements.push_back( cell - multipliers[ dimensions_in_which_cell_has_nonzero_length[i] ] );
        boundary_elements.push_back( cell + multipliers[ dimensions_in_which_cell_has_nonzero_length[i] ] );

        if (bdg) cerr << "multipliers[dimensions_in_which_cell_has_nonzero_length[i]] : " 
        << multipliers[dimensions_in_which_cell_has_nonzero_length[i]] << endl;
        if (bdg) cerr << "cell - multipliers[dimensions_in_which_cell_has_nonzero_length[i]] : " 
        << cell - multipliers[dimensions_in_which_cell_has_nonzero_length[i]] << endl;
        if (bdg) cerr << "cell + multipliers[dimensions_in_which_cell_has_nonzero_length[i]] : " 
        << cell + multipliers[dimensions_in_which_cell_has_nonzero_length[i]] << endl;
    }
    return boundary_elements;
}




template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::get_coboundary_of_a_cell( size_t cell )const
{
    bool bdg = false;
    //first of all, we need to take the list of coordinates in which the cell has nonzero length. 
    //We do it by using modified version to compute dimension of a cell:
    std::vector< unsigned > dimensions_in_which_cell_has_zero_length;
    unsigned dimension = 0;
    size_t cell1 = cell;
    for ( size_t i = this->multipliers.size() ; i != 0 ; --i )
    {
        unsigned position = cell1/multipliers[i-1];
        if ( position%2 == 0 )
        {
            dimensions_in_which_cell_has_zero_length.push_back(i-1);
            dimension++;
        }
        cell1 = cell1%multipliers[i-1];
    }

    std::vector<unsigned> counter = this->compute_counter_for_given_cell( cell );
    //reverse(counter.begin() , counter.end());

    if (bdg)
    {
        cerr << "dimensions_in_which_cell_has_zero_length : \n";
        for ( size_t i = 0 ; i != dimensions_in_which_cell_has_zero_length.size() ; ++i )
        {
            cerr << dimensions_in_which_cell_has_zero_length[i] << endl;
        }
        cerr << "\n counter : " << endl;
        for ( size_t i = 0 ; i != counter.size() ; ++i )
        {
            cerr << counter[i] << endl;
        }
        getchar();
    }

    std::vector< size_t > coboundary_elements;
    if ( dimensions_in_which_cell_has_zero_length.size() == 0 )return coboundary_elements;
    for ( size_t i = 0 ; i != dimensions_in_which_cell_has_zero_length.size() ; ++i )
    {
       if ( bdg )
        {
            cerr << "Dimension : " << i << endl;
            if (counter[dimensions_in_which_cell_has_zero_length[i]] == 0)
            {
                cerr << "In dimension : " << i 
                << " we cannot substract, since we will jump out of a Bitmap_cubical_complex_base \n";
            }
            if ( counter[dimensions_in_which_cell_has_zero_length[i]] 
                 == 
                 2*this->sizes[dimensions_in_which_cell_has_zero_length[i]] )
            {
                cerr << "In dimension : " << i 
                     << " we cannot substract, since we will jump out of a Bitmap_cubical_complex_base \n";
            }
        }


        if ( (cell > multipliers[dimensions_in_which_cell_has_zero_length[i]]) 
              && (counter[dimensions_in_which_cell_has_zero_length[i]] != 0) )
        //if ( counter[dimensions_in_which_cell_has_zero_length[i]] != 0 )
        {
            if ( bdg )
            {
                 cerr << "Subtracting : " << cell - multipliers[dimensions_in_which_cell_has_zero_length[i]] << endl;
            }
            coboundary_elements.push_back( cell - multipliers[dimensions_in_which_cell_has_zero_length[i]] );
        }
        if ( 
           (cell + multipliers[dimensions_in_which_cell_has_zero_length[i]] < this->data.size()) && 
           (counter[dimensions_in_which_cell_has_zero_length[i]] 
           != 
           2*this->sizes[dimensions_in_which_cell_has_zero_length[i]]) 
           )
        //if ( counter[dimensions_in_which_cell_has_zero_length[i]] != 
        //2*this->sizes[dimensions_in_which_cell_has_zero_length[i]] )
        {
            coboundary_elements.push_back( cell + multipliers[dimensions_in_which_cell_has_zero_length[i]] );
            if ( bdg )cerr << "Adding : " << cell + multipliers[dimensions_in_which_cell_has_zero_length[i]] << endl;
        }
    }
    return coboundary_elements;
}






template <typename T>
unsigned Bitmap_cubical_complex_base<T>::get_dimension_of_a_cell( size_t cell )const
{
    bool dbg = false;
    if (dbg)cerr << "\n\n\n Computing position o a cell of an index : " << cell << endl;
    unsigned dimension = 0;
    for ( size_t i = this->multipliers.size() ; i != 0 ; --i )
    {
        unsigned position = cell/multipliers[i-1];

        if (dbg)cerr << "i-1 :" << i-1 << endl;
        if (dbg)cerr << "cell : " << cell << endl;
        if (dbg)cerr << "position : " << position << endl;
        if (dbg)cerr << "multipliers["<<i-1<<"] = " << multipliers[i-1] << endl;
        if (dbg)getchar();

        if ( position%2 == 1 )
        {
            if (dbg)cerr << "Nonzero length in this direction \n";
            dimension++;
        }
        cell = cell%multipliers[i-1];
    }
    return dimension;
}

template <typename T>
T& Bitmap_cubical_complex_base<T>::get_cell_data( size_t cell )
{
    return this->data[cell];
}


template <typename T>
void Bitmap_cubical_complex_base<T>::impose_lower_star_filtration()
{
    bool dbg = false;

    //this vector will be used to check which elements have already been taken care of 
    //in imposing lower star filtration:
    std::vector<bool> is_this_cell_considered( this->data.size() , false );

    std::vector<size_t> indices_to_consider;
    //we assume here that we already have a filtration on the top dimensional cells and 
    //we have to extend it to lower ones.
    typename Bitmap_cubical_complex_base<T>::Top_dimensional_cells_iterator it(*this);
    for ( it = this->top_dimensional_cells_begin() ; it != this->top_dimensional_cells_end() ; ++it )
    {
        indices_to_consider.push_back( it.compute_index_in_bitmap() );
    }

    while ( indices_to_consider.size() )
    {
        if ( dbg )
        {
            cerr << "indices_to_consider in this iteration \n";
            for ( size_t i = 0 ; i != indices_to_consider.size() ; ++i )
            {
                cout << indices_to_consider[i] << "  ";
            }
            getchar();
        }
        std::vector<size_t> new_indices_to_consider;
        for ( size_t i = 0 ; i != indices_to_consider.size() ; ++i )
        {
            std::vector<size_t> bd = this->get_boundary_of_a_cell( indices_to_consider[i] );
            for ( size_t boundaryIt = 0 ; boundaryIt != bd.size() ; ++boundaryIt )
            {
                if ( this->data[ bd[boundaryIt] ] > this->data[ indices_to_consider[i] ] )
                {
                    this->data[ bd[boundaryIt] ] = this->data[ indices_to_consider[i] ];
                }
                if ( is_this_cell_considered[ bd[boundaryIt] ] == false )
                {
                    new_indices_to_consider.push_back( bd[boundaryIt] );
                    is_this_cell_considered[ bd[boundaryIt] ] = true;
                }
            }
        }
        indices_to_consider.swap(new_indices_to_consider);
    }
}


template <typename T>
bool compareFirstElementsOfTuples( const std::pair< std::pair< T , size_t > , char >& first ,
                                   const std::pair< std::pair< T , size_t > , char >& second )
{
    if ( first.first.first < second.first.first )
    {
        return true;
    }
    else
    {
        if ( first.first.first > second.first.first )
        {
            return false;
        }
        //in this case first.first.first == second.first.first, so we need to compare dimensions
        return first.second < second.second;
    }
}



template <typename T>
std::vector< size_t > Bitmap_cubical_complex_base<T>::
generate_vector_of_shifts_for_bitmaps_with_periodic_boundary_conditions
( std::vector< bool >& directions_for_periodic_b_cond )
{
    bool dbg = false;
    if ( this->sizes.size() != directions_for_periodic_b_cond.size() )
    throw "directions_for_periodic_b_cond vector size is different from the size of the bitmap. Program terminate \n";

    std::vector<unsigned> sizes( this->sizes.size() );
    for ( size_t i = 0 ; i != this->sizes.size() ; ++i )sizes[i] = 2*this->sizes[i];

    counter c( sizes );

    std::vector< size_t > result;

    for ( size_t i = 0 ; i != this->data.size() ; ++i )
    {
        size_t position;
        if ( !c.isFinal() )
        {
            position = i;
            //result.push_back( i );
        }
        else
        {
            std::vector< bool > finals = c.directions_of_finals();
            bool jump_in_position = false;
            for ( size_t dir = 0 ; dir != finals.size() ; ++dir )
            {
                if ( finals[dir] == false )continue;
                if ( directions_for_periodic_b_cond[dir] )
                {
                    jump_in_position = true;
                }
            }
            if ( jump_in_position == true )
            {
                //in this case this guy is final, so we need to find 'the opposite one'
                position = compute_position_in_bitmap(  c.find_opposite( directions_for_periodic_b_cond ) );
            }
            else
            {
                position = i;
            }
        }
        result.push_back( position );
        if ( dbg )
        {
            cerr << " position : " << position << endl;
            cerr << c << endl;
            getchar();
        }

        c.increment();
    }

    return result;
}

}

}