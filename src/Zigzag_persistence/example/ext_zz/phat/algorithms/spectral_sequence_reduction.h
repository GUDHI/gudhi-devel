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
#include <phat/boundary_matrix.h>

namespace phat {
    class spectral_sequence_reduction {
    public:
        template< typename Representation >
        void operator () ( boundary_matrix< Representation >& boundary_matrix ) {

            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< index > lowest_one_lookup( nr_columns, -1 );

            //const index num_stripes = (index) sqrt( (double)nr_columns );
            const index num_stripes = omp_get_max_threads();

            index block_size = ( nr_columns % num_stripes == 0 ) ? nr_columns / num_stripes : block_size = nr_columns / num_stripes + 1;
            
            std::vector< std::vector< index > > unreduced_cols_cur_pass( num_stripes );
            std::vector< std::vector< index > > unreduced_cols_next_pass( num_stripes );
            
            for( index cur_dim = boundary_matrix.get_max_dim(); cur_dim >= 1 ; cur_dim-- ) {
                #pragma omp parallel for schedule( guided, 1 )
                for( index cur_stripe = 0; cur_stripe < num_stripes; cur_stripe++ ) {
                    index col_begin = cur_stripe * block_size;
                    index col_end = std::min( (cur_stripe+1) * block_size, nr_columns );
                    for( index cur_col = col_begin; cur_col < col_end; cur_col++ )
                        if( boundary_matrix.get_dim( cur_col ) == cur_dim && boundary_matrix.get_max_index( cur_col ) != -1 )
                            unreduced_cols_cur_pass[ cur_stripe ].push_back( cur_col );
                }
                for( index cur_pass = 0; cur_pass < num_stripes; cur_pass++ ) {
                    boundary_matrix.sync();
                    #pragma omp parallel for schedule( guided, 1 )
                    for( int cur_stripe = 0; cur_stripe < num_stripes; cur_stripe++ ) {
                        index row_begin = (cur_stripe - cur_pass) * block_size;
                        index row_end = row_begin + block_size;
                        unreduced_cols_next_pass[ cur_stripe ].clear();
                        for( index idx = 0; idx < (index)unreduced_cols_cur_pass[ cur_stripe ].size(); idx++ ) {
                            index cur_col = unreduced_cols_cur_pass[ cur_stripe ][ idx ];
                            index lowest_one = boundary_matrix.get_max_index( cur_col );
                            while( lowest_one != -1 && lowest_one >= row_begin && lowest_one < row_end  && lowest_one_lookup[ lowest_one ] != -1 ) {
                                boundary_matrix.add_to( lowest_one_lookup[ lowest_one ], cur_col );
                                lowest_one = boundary_matrix.get_max_index( cur_col );
                            }
                            if( lowest_one != -1 ) {
                                if( lowest_one >= row_begin && lowest_one < row_end ) {
                                    lowest_one_lookup[ lowest_one ] = cur_col;
                                    boundary_matrix.clear( lowest_one );
                                    boundary_matrix.finalize( cur_col );
                                } else {
                                    unreduced_cols_next_pass[ cur_stripe ].push_back( cur_col );
                                }
                            }
                        }
                        unreduced_cols_next_pass[ cur_stripe ].swap( unreduced_cols_cur_pass[ cur_stripe ] );
                    }
                }
            }
        }
    };
}
