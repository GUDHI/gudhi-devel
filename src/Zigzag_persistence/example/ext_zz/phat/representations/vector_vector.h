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

#include <phat/helpers/misc.h>

namespace phat {
    class vector_vector {

    protected:
        std::vector< dimension > dims;
        std::vector< column > matrix;

        thread_local_storage< column > temp_column_buffer;

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
            col = matrix[ idx ]; 
        }
        void _set_col( index idx, const column& col  ) { 
            matrix[ idx ] = col; 
        }

        // true iff boundary of given idx is empty
        bool _is_empty( index idx ) const { 
            return matrix[ idx ].empty(); 
        }

        // largest row index of given column idx (new name for lowestOne())
        index _get_max_index( index idx ) const { 
            return matrix[ idx ].empty() ? -1 : matrix[ idx ].back(); 
        }

        // removes the maximal index of a column
        void _remove_max( index idx ) {
            matrix[ idx ].pop_back();
        }

        // clears given column
        void _clear( index idx ) { 
            matrix[ idx ].clear();
        }

        // syncronizes all data structures (essential for openmp stuff)
        void _sync() {}

        // adds column 'source' to column 'target'
        void _add_to( index source, index target ) {
            column& source_col = matrix[ source ];
            column& target_col = matrix[ target ];
            column& temp_col = temp_column_buffer();
            
            
            size_t new_size = source_col.size() + target_col.size();
            
            if (new_size > temp_col.size()) temp_col.resize(new_size);
            
            std::vector<index>::iterator col_end = std::set_symmetric_difference( target_col.begin(), target_col.end(),
                                           source_col.begin(), source_col.end(),
                                                                                 temp_col.begin() );
            temp_col.erase(col_end, temp_col.end());

            
            target_col.swap(temp_col);
        }
        
        // finalizes given column
        void _finalize( index idx ) {
            column& col = matrix[ idx ];
            column(col.begin(), col.end()).swap(col);
        }
    };
}
