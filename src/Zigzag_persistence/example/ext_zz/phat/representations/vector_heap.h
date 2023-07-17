/*  Copyright 2013 IST Austria
Contributed by: Jan Reininghaus

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

#include <phat/helpers/misc.h>

namespace phat {
    class vector_heap {

    protected:
        std::vector< dimension > dims;
        std::vector< column > matrix;

        std::vector< index > inserts_since_last_prune;

        mutable thread_local_storage< column > temp_column_buffer;

    protected:
        void _prune( index idx )
        {
            column& col = matrix[ idx ];
            column& temp_col = temp_column_buffer();
            temp_col.clear();
            index max_index = _pop_max_index( col );
            while( max_index != -1 ) {
                temp_col.push_back( max_index );
                max_index = _pop_max_index( col );
            }
            col = temp_col;
            std::reverse( col.begin( ), col.end( ) );
            std::make_heap( col.begin( ), col.end( ) );
            inserts_since_last_prune[ idx ] = 0;
        }

        index _pop_max_index( index idx )
        {
            return _pop_max_index( matrix[ idx ] );
        }

        index _pop_max_index( column& col ) const
        {
            if( col.empty( ) )
                return -1;
            else {
                index max_element = col.front( );
                std::pop_heap( col.begin( ), col.end( ) );
                col.pop_back( );
                while( !col.empty( ) && col.front( ) == max_element ) {
                    std::pop_heap( col.begin( ), col.end( ) );
                    col.pop_back( );
                    if( col.empty( ) )
                        return -1;
                    else {
                        max_element = col.front( );
                        std::pop_heap( col.begin( ), col.end( ) );
                        col.pop_back( );
                    }
                }
                return max_element;
            }
        }

    public:
        // overall number of cells in boundary_matrix
        index _get_num_cols( ) const
        {
            return (index)matrix.size( );
        }
        void _set_num_cols( index nr_of_columns )
        {
            dims.resize( nr_of_columns );
            matrix.resize( nr_of_columns );
            inserts_since_last_prune.assign( nr_of_columns, 0 );
        }

        // dimension of given index
        dimension _get_dim( index idx ) const
        {
            return dims[ idx ];
        }
        void _set_dim( index idx, dimension dim )
        {
            dims[ idx ] = dim;
        }

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( index idx, column& col ) const
        {
            temp_column_buffer( ) = matrix[ idx ];
            
            index max_index = _pop_max_index( temp_column_buffer() );
            while( max_index != -1 ) {
                col.push_back( max_index );
                max_index = _pop_max_index( temp_column_buffer( ) );
            }
            std::reverse( col.begin( ), col.end( ) );
        }
        void _set_col( index idx, const column& col )
        {
            matrix[ idx ] = col;
            std::make_heap( matrix[ idx ].begin( ), matrix[ idx ].end( ) );
        }

        // true iff boundary of given idx is empty
        bool _is_empty( index idx ) const
        {
            return _get_max_index( idx ) == -1;
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index( index idx ) const
        {
            column& col = const_cast< column& >( matrix[ idx ] );
            index max_element = _pop_max_index( col );
            col.push_back( max_element );
            std::push_heap( col.begin( ), col.end( ) );
            return max_element;
        }

        // removes the maximal index of a column
        void _remove_max( index idx )
        {
            _pop_max_index( idx );
        }

        // clears given column
        void _clear( index idx )
        {
            matrix[ idx ].clear( );
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync( ) {}

        // adds column 'source' to column 'target'
        void _add_to( index source, index target )
        {              
            for( index idx = 0; idx < (index)matrix[ source ].size( ); idx++ ) {
                matrix[ target ].push_back( matrix[ source ][ idx ] );
                std::push_heap( matrix[ target ].begin(), matrix[ target ].end() );
            }
            inserts_since_last_prune[ target ] += matrix[ source ].size();

            if( 2 * inserts_since_last_prune[ target ] > ( index )matrix[ target ].size() )
                _prune( target );
        }
        
        // finalizes given column
        void _finalize( index idx ) {
            _prune( idx );
        }

    };
}
