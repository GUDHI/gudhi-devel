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
    class vector_list {

    protected:
        std::vector< dimension > dims;
        std::vector< std::list< index > > matrix;

    public:
        // overall number of cells in boundary_matrix
        index _get_num_cols() const {
            return (index)matrix.size(); 
        }
        void _set_num_cols( index nr_of_columns ) {
            dims.resize( nr_of_columns );
            matrix.resize( nr_of_columns );
        }

         // dimension of given index
        dimension _get_dim( index idx ) const { 
            return dims[ idx ]; 
        }
        void _set_dim( index idx, dimension dim ) { 
            dims[ idx ] = dim; 
        }

        // replaces(!) content of 'col' with boundary of given index
        void _get_col( index idx, column& col  ) const {
            col.clear();
            col.reserve( matrix[idx].size() );
            std::copy (matrix[idx].begin(), matrix[idx].end(), std::back_inserter(col) );
        }

        void _set_col( index idx, const column& col  ) {
            matrix[ idx ].clear();
            matrix[ idx ].resize( col.size() );
            std::copy (col.begin(), col.end(), matrix[ idx ].begin() );
        }

        // true iff boundary of given idx is empty
        bool _is_empty( index idx ) const {
            return matrix[ idx ].empty();    
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index( index idx ) const {
            return matrix[ idx ].empty() ? -1 : *matrix[ idx ].rbegin();
        }

        // removes the maximal index of a column
        void _remove_max( index idx ) {
            std::list< index >::iterator it = matrix[ idx ].end();
            it--;
            matrix[ idx ].erase( it );
        }
        
        // clears given column
        void _clear( index idx ) {
            matrix[ idx ].clear();    
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync() {}

        // adds column 'source' to column 'target'
        void _add_to( index source, index target ) {
            std::list< index >& source_col = matrix[ source ];
            std::list< index >& target_col = matrix[ target ];
            std::list< index > temp_col;
            target_col.swap( temp_col );
            std::set_symmetric_difference( temp_col.begin(), temp_col.end(),
                                           source_col.begin(), source_col.end(),
                                           std::back_inserter( target_col ) );
        }
        
        // finalizes given column
        void _finalize( index idx ) {
        }
    };
}
