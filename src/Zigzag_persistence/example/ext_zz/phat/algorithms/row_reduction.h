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
#include <phat/boundary_matrix.h>

namespace phat {
    class row_reduction {
    public:
        template< typename Representation >
        void operator() ( boundary_matrix< Representation >& boundary_matrix ) {
            
            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< std::vector< index > > lowest_one_lookup( nr_columns );
            
            for( index cur_col = nr_columns - 1; cur_col >= 0; cur_col-- ) {
                if( !boundary_matrix.is_empty( cur_col ) )
                    lowest_one_lookup[ boundary_matrix.get_max_index( cur_col ) ].push_back( cur_col );
                
                if( !lowest_one_lookup[ cur_col ].empty() ) {
                    boundary_matrix.clear( cur_col );
                    boundary_matrix.finalize( cur_col );
                    std::vector< index >& cols_with_cur_lowest = lowest_one_lookup[ cur_col ];
                    index source = *min_element( cols_with_cur_lowest.begin(), cols_with_cur_lowest.end() );
                    for( index idx = 0; idx < (index)cols_with_cur_lowest.size(); idx++ ) {
                        index target = cols_with_cur_lowest[ idx ];
                        if( target != source && !boundary_matrix.is_empty( target ) ) {
                            boundary_matrix.add_to( source, target );     
                            if( !boundary_matrix.is_empty( target ) ) {
                                index lowest_one_of_target = boundary_matrix.get_max_index( target );
                                lowest_one_lookup[ lowest_one_of_target ].push_back( target );
                            }
                        } 
                    }
                }
            }
        }
    };
}