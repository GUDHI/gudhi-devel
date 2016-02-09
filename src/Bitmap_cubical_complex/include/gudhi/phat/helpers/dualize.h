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

#include "misc.h"
#include "../boundary_matrix.h"

namespace phat {
    template< typename Representation >
    void dualize( boundary_matrix< Representation >& boundary_matrix ) {

        std::vector< dimension > dual_dims;
        std::vector< std::vector< index > > dual_matrix;

        index nr_of_columns = boundary_matrix.get_num_cols();
        dual_matrix.resize( nr_of_columns );
        dual_dims.resize( nr_of_columns );

        std::vector< index > dual_sizes( nr_of_columns, 0 );

        column temp_col;
        for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ ) {
            boundary_matrix.get_col( cur_col, temp_col );
            for( index idx = 0; idx < (index)temp_col.size(); idx++)
                dual_sizes[ nr_of_columns - 1 - temp_col[ idx ]  ]++;
        }

        for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ ) {
            dual_matrix[cur_col].reserve(dual_sizes[cur_col]);
        }

        for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ ) {
            boundary_matrix.get_col( cur_col, temp_col );
            for( index idx = 0; idx < (index)temp_col.size(); idx++)
                dual_matrix[ nr_of_columns - 1 - temp_col[ idx ]  ].push_back( nr_of_columns - 1 - cur_col );
        }

        const dimension max_dim = boundary_matrix.get_max_dim();
        for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ )
            dual_dims[ nr_of_columns - 1 - cur_col ] = max_dim - boundary_matrix.get_dim( cur_col );

        for( index cur_col = 0; cur_col < nr_of_columns; cur_col++ )
            std::reverse( dual_matrix[ cur_col ].begin(), dual_matrix[ cur_col ].end() );

        boundary_matrix.load_vector_vector( dual_matrix, dual_dims );
    }
}
