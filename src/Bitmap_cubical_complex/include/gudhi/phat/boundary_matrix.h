/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include "helpers/misc.h"
#include "representations/bit_tree_pivot_column.h"

// interface class for the main data structure -- implementations of the interface can be found in ./representations
namespace phat {
    template< class Representation = bit_tree_pivot_column >
    class boundary_matrix
    {

    protected:
        Representation rep;

    // interface functions -- actual implementation and complexity depends on chosen @Representation template
    public:
        // get overall number of columns in boundary_matrix
        index get_num_cols() const { return rep._get_num_cols(); }

        // set overall number of columns in boundary_matrix
        void set_num_cols( index nr_of_columns ) { rep._set_num_cols( nr_of_columns ); }

        // get dimension of given index
        dimension get_dim( index idx ) const { return rep._get_dim( idx ); }

        // set dimension of given index
        void set_dim( index idx, dimension dim ) { rep._set_dim( idx, dim ); }

        // replaces content of @col with boundary of given index
        void get_col( index idx, column& col  ) const { rep._get_col( idx, col ); }

        // set column @idx to the values contained in @col
        void set_col( index idx, const column& col  ) { rep._set_col( idx, col ); }

        // true iff boundary of given column is empty
        bool is_empty( index idx ) const { return rep._is_empty( idx ); }

        // largest index of given column (new name for lowestOne())
        index get_max_index( index idx ) const { return rep._get_max_index( idx ); }

        // removes maximal index from given column
        void remove_max( index idx ) { rep._remove_max( idx ); }

        // adds column @source to column @target'
        void add_to( index source, index target ) { rep._add_to( source, target ); }

        // clears given column
        void clear( index idx ) { rep._clear( idx ); }

        // syncronizes all internal data structures -- has to be called before and after any multithreaded access!
        void sync() { rep._sync(); }

    // info functions -- independent of chosen 'Representation'
    public:
        // maximal dimension
        dimension get_max_dim() const {
            dimension cur_max_dim = 0;
            for( index idx = 0; idx < get_num_cols(); idx++ )
                cur_max_dim = get_dim( idx ) > cur_max_dim ? get_dim( idx ) : cur_max_dim;
            return cur_max_dim;
        }

        // number of nonzero rows for given column @idx
        index get_num_rows( index idx ) const {
            column cur_col;
            get_col( idx, cur_col );
            return cur_col.size();
        }

        // maximal number of nonzero rows of all columns
        index get_max_col_entries() const {
            index max_col_entries = -1;
            const index nr_of_columns = get_num_cols();
            for( index idx = 0; idx < nr_of_columns; idx++ )
                max_col_entries = get_num_rows( idx ) > max_col_entries ? get_num_rows( idx ) : max_col_entries;
            return max_col_entries;
        }

        // maximal number of nonzero cols of all rows
        index get_max_row_entries() const {
            size_t max_row_entries = 0;
            const index nr_of_columns = get_num_cols();
            std::vector< std::vector< index > > transposed_matrix( nr_of_columns );
            column temp_col;
            for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ ) {
                get_col( cur_col, temp_col );
                for( index idx = 0; idx < (index)temp_col.size(); idx++)
                    transposed_matrix[ temp_col[ idx ]  ].push_back( cur_col );
            }
            for( index idx = 0; idx < nr_of_columns; idx++ )
                max_row_entries = transposed_matrix[ idx ].size() > max_row_entries ? transposed_matrix[ idx ].size() : max_row_entries;
            return max_row_entries;
        }

        // overall number of entries in the matrix
        index get_num_entries() const {
            index number_of_nonzero_entries = 0;
            const index nr_of_columns = get_num_cols();
            for( index idx = 0; idx < nr_of_columns; idx++ )
                number_of_nonzero_entries += get_num_rows( idx );
            return number_of_nonzero_entries;
        }

    // operators / constructors
    public:
        boundary_matrix() {};

        template< class OtherRepresentation >
        boundary_matrix( const boundary_matrix< OtherRepresentation >& other ) {
            *this = other;
        }

        template< typename OtherRepresentation >
        bool operator==( const boundary_matrix< OtherRepresentation >& other_boundary_matrix ) const {
            const index number_of_columns = this->get_num_cols();

            if( number_of_columns != other_boundary_matrix.get_num_cols() )
                return false;

            column temp_col;
            column other_temp_col;
            for( index idx = 0; idx < number_of_columns; idx++ ) {
                this->get_col( idx, temp_col );
                other_boundary_matrix.get_col( idx, other_temp_col );
                if( temp_col != other_temp_col || this->get_dim( idx ) != other_boundary_matrix.get_dim( idx ) )
                    return false;
            }
            return true;
        }

        template< typename OtherRepresentation >
        bool operator!=( const boundary_matrix< OtherRepresentation >& other_boundary_matrix ) const {
            return !( *this == other_boundary_matrix );
        }

        template< typename OtherRepresentation >
        boundary_matrix< Representation >& operator=( const boundary_matrix< OtherRepresentation >& other )
        {
            const index nr_of_columns = other.get_num_cols();
            this->set_num_cols( nr_of_columns );
            column temp_col;
            for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                this->set_dim( cur_col, other.get_dim( cur_col ) );
                other.get_col( cur_col, temp_col );
                this->set_col( cur_col, temp_col );
            }

            // by convention, always return *this
            return *this;
        }

    // I/O -- independent of chosen 'Representation'
    public:

        // initializes boundary_matrix from (vector<vector>, vector) pair -- untested
        template< typename index_type, typename dimemsion_type >
        void load_vector_vector( const std::vector< std::vector< index_type > >& input_matrix, const std::vector< dimemsion_type >& input_dims ) {
            const index nr_of_columns = (index)input_matrix.size();
            this->set_num_cols( nr_of_columns );
            column temp_col;
            for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                this->set_dim( cur_col, (dimension)input_dims[ cur_col ] );

                index num_rows = input_matrix[ cur_col ].size();
                temp_col.resize( num_rows );
                for( index cur_row = 0; cur_row <  num_rows; cur_row++ )
                    temp_col[ cur_row ] = (index)input_matrix[ cur_col ][ cur_row ];
                this->set_col( cur_col, temp_col );
            }
        }

        template< typename index_type, typename dimemsion_type >
        void save_vector_vector( std::vector< std::vector< index_type > >& output_matrix, std::vector< dimemsion_type >& output_dims ) {
            const index nr_of_columns = get_num_cols();
            output_matrix.resize( nr_of_columns );
            output_dims.resize( nr_of_columns );
            column temp_col;
            for( index cur_col = 0; cur_col <  nr_of_columns; cur_col++ ) {
                output_dims[ cur_col ] = (dimemsion_type)get_dim( cur_col );
                get_col( cur_col, temp_col );
                index num_rows = temp_col.size();
                output_matrix[ cur_col ].clear();
                output_matrix[ cur_col ].resize( num_rows );
                for( index cur_row = 0; cur_row <  num_rows; cur_row++ )
                    output_matrix[ cur_col ][ cur_row ] = (index_type)temp_col[ cur_row ];
            }
        }

        // Loads the boundary_matrix from given file in ascii format
        // Format: each line represents a column, first number is dimension, other numbers are the content of the column.
        // Ignores empty lines and lines starting with a '#'.
        bool load_ascii( std::string filename ) {
            // first count number of columns:
            std::string cur_line;
            std::ifstream dummy( filename .c_str() );
            if( dummy.fail() )
                return false;

            index number_of_columns = 0;
            while( getline( dummy, cur_line ) ) {
                cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
                if( cur_line != "" && cur_line[ 0 ] != '#' )
                    number_of_columns++;

            }
            this->set_num_cols( number_of_columns );
            dummy.close();

            std::ifstream input_stream( filename.c_str() );
            if( input_stream.fail() )
                return false;

            column temp_col;
            index cur_col = -1;
            while( getline( input_stream, cur_line ) ) {
                cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
                if( cur_line != "" && cur_line[ 0 ] != '#' ) {
                    cur_col++;
                    std::stringstream ss( cur_line );

                    int64_t temp_dim;
                    ss >> temp_dim;
                    this->set_dim( cur_col, (dimension) temp_dim );

                    int64_t temp_index;
                    temp_col.clear();
                    while( ss.good() ) {
                        ss >> temp_index;
                        temp_col.push_back( (index)temp_index );
                    }
                    std::sort( temp_col.begin(), temp_col.end() );
                    this->set_col( cur_col, temp_col );
                }
            }

            input_stream.close();
            return true;
        }

        // Saves the boundary_matrix to given file in ascii format
        // Format: each line represents a column, first number is dimension, other numbers are the content of the column
        bool save_ascii( std::string filename ) {
            std::ofstream output_stream( filename.c_str() );
            if( output_stream.fail() )
                return false;

            const index nr_columns = this->get_num_cols();
            column tempCol;
            for( index cur_col = 0; cur_col < nr_columns; cur_col++ ) {
                output_stream << (int64_t)this->get_dim( cur_col );
                this->get_col( cur_col, tempCol );
                for( index cur_row_idx = 0; cur_row_idx < (index)tempCol.size(); cur_row_idx++ )
                    output_stream << " " << tempCol[ cur_row_idx ];
                output_stream << std::endl;
            }

            output_stream.close();
            return true;
        }

        // Loads boundary_matrix from given file
        // Format: nr_columns % dim1 % N1 % row1 row2 % ...% rowN1 % dim2 % N2 % ...
        bool load_binary( std::string filename ) {
            std::ifstream input_stream( filename.c_str(), std::ios_base::binary | std::ios_base::in );
            if( input_stream.fail() )
                return false;

            int64_t nr_columns;
            input_stream.read( (char*)&nr_columns, sizeof( int64_t ) );
            this->set_num_cols( (index)nr_columns );

            column temp_col;
            for( index cur_col = 0; cur_col < nr_columns; cur_col++ ) {
                int64_t cur_dim;
                input_stream.read( (char*)&cur_dim, sizeof( int64_t ) );
                this->set_dim( cur_col, (dimension) cur_dim );
                int64_t nr_rows;
                input_stream.read( (char*)&nr_rows, sizeof( int64_t ) );
                temp_col.resize( (std::size_t)nr_rows );
                for( index idx = 0; idx < nr_rows; idx++ ) {
                    int64_t cur_row;
                    input_stream.read( (char*)&cur_row, sizeof( int64_t ) );
                    temp_col[ idx ] = (index)cur_row;
                }
                this->set_col( cur_col, temp_col );
            }

            input_stream.close();
            return true;
        }

        // Saves the boundary_matrix to given file in binary format
        // Format: nr_columns % dim1 % N1 % row1 row2 % ...% rowN1 % dim2 % N2 % ...
        bool save_binary( std::string filename ) {
            std::ofstream output_stream( filename.c_str(), std::ios_base::binary | std::ios_base::out );
            if( output_stream.fail() )
                return false;

            const int64_t nr_columns = this->get_num_cols();
            output_stream.write( (char*)&nr_columns, sizeof( int64_t ) );
            column tempCol;
            for( index cur_col = 0; cur_col < nr_columns; cur_col++ ) {
                int64_t cur_dim = this->get_dim( cur_col );
                output_stream.write( (char*)&cur_dim, sizeof( int64_t ) );
                this->get_col( cur_col, tempCol );
                int64_t cur_nr_rows = tempCol.size();
                output_stream.write( (char*)&cur_nr_rows, sizeof( int64_t ) );
                for( index cur_row_idx = 0; cur_row_idx < (index)tempCol.size(); cur_row_idx++ ) {
                    int64_t cur_row = tempCol[ cur_row_idx ];
                    output_stream.write( (char*)&cur_row, sizeof( int64_t ) );
                }
            }

            output_stream.close();
            return true;
        }
    };
}
