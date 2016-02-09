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

#include "../helpers/misc.h"
#include "../boundary_matrix.h"

namespace phat {
    class chunk_reduction {
    public:
        enum column_type { GLOBAL
                         , LOCAL_POSITIVE
                         , LOCAL_NEGATIVE };

    public:
        template< typename Representation >
        void operator() ( boundary_matrix< Representation >& boundary_matrix ) {

            const index nr_columns = boundary_matrix.get_num_cols();
            const dimension max_dim = boundary_matrix.get_max_dim();

            std::vector< index > lowest_one_lookup( nr_columns, -1 );
            std::vector < column_type > column_type( nr_columns, GLOBAL );
            std::vector< char > is_active( nr_columns, false );

            std::vector<index> chunk_boundaries;
            _get_chunks( boundary_matrix, chunk_boundaries );

            // Phase 1: Reduce chunks locally -- 1st pass
            #pragma omp parallel for schedule( guided, 1 )
            for( index chunk_id = 0; chunk_id < (index) chunk_boundaries.size() - 2; chunk_id += 2 )
                _local_chunk_reduction( boundary_matrix, lowest_one_lookup, column_type, max_dim,
                                        chunk_boundaries[chunk_id], chunk_boundaries[chunk_id+2] - 1 );
            boundary_matrix.sync();

            // Phase 1: Reduce chunks locally -- 2nd pass
            #pragma omp parallel for schedule( guided, 1 )
            for( index chunk_id = 1; chunk_id < (index) chunk_boundaries.size() - 2; chunk_id += 2 )
                _local_chunk_reduction( boundary_matrix, lowest_one_lookup, column_type, max_dim,
                                        chunk_boundaries[chunk_id], chunk_boundaries[chunk_id+2] - 1 );
            boundary_matrix.sync();

            // get global columns
            std::vector< index > global_columns;
            for( index cur_col_idx = 0; cur_col_idx < nr_columns; cur_col_idx++ )
                if( column_type[ cur_col_idx ] == GLOBAL )
                    global_columns.push_back( cur_col_idx );

            // get active columns
            #pragma omp parallel for
            for( index idx = 0; idx < (index)global_columns.size(); idx++ )
                is_active[ global_columns[ idx ] ] = true;
            _get_active_columns( boundary_matrix, lowest_one_lookup, column_type, global_columns, is_active );

            // Phase 2+3: Simplify columns and reduce them
            for( dimension cur_dim = max_dim; cur_dim >= 1; cur_dim-- ) {
                // Phase 2: Simplify columns
                std::vector< index > temp_col;
                #pragma omp parallel for schedule( guided, 1 ), private( temp_col )
                for( index idx = 0; idx < (index)global_columns.size(); idx++ )
                    if( boundary_matrix.get_dim( global_columns[ idx ] ) == cur_dim )
                        _global_column_simplification( global_columns[ idx ], boundary_matrix, lowest_one_lookup, column_type, is_active, temp_col );
                boundary_matrix.sync();

                // Phase 3: Reduce columns
                for( index idx = 0; idx < (index)global_columns.size(); idx++ ) {
                    index cur_col = global_columns[ idx ];
                    if( boundary_matrix.get_dim( cur_col ) == cur_dim && column_type[ cur_col ] == GLOBAL ) {
                        index lowest_one = boundary_matrix.get_max_index( cur_col );
                        while( lowest_one != -1 && lowest_one_lookup[ lowest_one ] != -1 ) {
                            boundary_matrix.add_to( lowest_one_lookup[ lowest_one ], cur_col );
                            lowest_one = boundary_matrix.get_max_index( cur_col );
                        }
                        if( lowest_one != -1 ) {
                            lowest_one_lookup[ lowest_one ] = cur_col;
                            boundary_matrix.clear( lowest_one );
                        }
                    }
                }
            }

            boundary_matrix.sync();
        }

    protected:
        template< typename Representation >
        void _get_chunks( const boundary_matrix< Representation >& boundary_matrix
                        , std::vector< index >& chunk_boundaries)
        {
            chunk_boundaries.clear();
            std::vector<index> temp_chunk_boundaries;
            const index nr_columns = boundary_matrix.get_num_cols();

            // size of chuks = sqrt(N)
            const index chunk_size = (index) sqrt( (float)nr_columns );

            // size of chunks = N / num_threads
            //const index chunk_size = nr_columns / omp_get_max_threads();

            for ( index cur_col = 0; cur_col < nr_columns; cur_col++ )
                if( cur_col % chunk_size == 0 )
                    temp_chunk_boundaries.push_back( cur_col );
            temp_chunk_boundaries.push_back( nr_columns );

            // subdivide chunks for interleaved 2 pass appraoch
            for( index chunk_id = 0; chunk_id < (index) temp_chunk_boundaries.size(); chunk_id ++ ) {
                chunk_boundaries.push_back( temp_chunk_boundaries[ chunk_id ] );
                if( chunk_id < (index) temp_chunk_boundaries.size() - 1 ) {
                    index midPoint = ( temp_chunk_boundaries[ chunk_id ] + temp_chunk_boundaries[ chunk_id + 1 ] ) / 2;
                    chunk_boundaries.push_back( midPoint );
                }
            }
        }

        template< typename Representation >
        void _local_chunk_reduction( boundary_matrix< Representation >& boundary_matrix
                                   , std::vector<index>& lowest_one_lookup
                                   , std::vector< column_type >& column_type
                                   , const dimension max_dim
                                   , const index chunk_begin
                                   , const index chunk_end ) {
            for( dimension cur_dim = max_dim; cur_dim >= 1; cur_dim-- ) {
                for( index cur_col = chunk_begin; cur_col <= chunk_end; cur_col++ ) {
                    if( column_type[ cur_col ] == GLOBAL && boundary_matrix.get_dim( cur_col ) == cur_dim ) {
                        index lowest_one = boundary_matrix.get_max_index( cur_col );
                        while( lowest_one != -1 && lowest_one >= chunk_begin && lowest_one_lookup[ lowest_one ] != -1  ) {
                            boundary_matrix.add_to( lowest_one_lookup[ lowest_one ], cur_col );
                            lowest_one = boundary_matrix.get_max_index( cur_col );
                        }
                        if( lowest_one >= chunk_begin ) {
                            lowest_one_lookup[ lowest_one ] = cur_col;
                            column_type[ cur_col ] = LOCAL_NEGATIVE;
                            column_type[ lowest_one ] = LOCAL_POSITIVE;
                            boundary_matrix.clear( lowest_one );
                        }
                    }
                }
            }
        }

        template< typename Representation >
        void _get_active_columns( const boundary_matrix< Representation >& boundary_matrix
                                , const std::vector< index >& lowest_one_lookup
                                , const std::vector< column_type >& column_type
                                , const std::vector< index >& global_columns
                                , std::vector< char >& is_active ) {

            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< char > finished( nr_columns, false );

            std::vector< std::pair < index, index > > stack;
            std::vector< index > cur_col_values;
            #pragma omp parallel for schedule( guided, 1 ), private( stack, cur_col_values )
            for( index idx = 0; idx < (index)global_columns.size(); idx++ ) {
                bool pop_next = false;
                index start_col = global_columns[ idx ];
                stack.push_back( std::pair< index, index >( start_col, -1 ) );
                while( !stack.empty() ) {
                    index cur_col = stack.back().first;
                    index prev_col = stack.back().second;
                    if( pop_next ) {
                        stack.pop_back();
                        pop_next = false;
                        if( prev_col != -1 ) {
                            if( is_active[ cur_col ] ) {
                                is_active[ prev_col ] = true;
                            }
                            if( prev_col == stack.back().first ) {
                                finished[ prev_col ] = true;
                                pop_next = true;
                            }
                        }
                    } else {
                        pop_next = true;
                        boundary_matrix.get_col( cur_col, cur_col_values );
                        for( index idx = 0; idx < (index) cur_col_values.size(); idx++ ) {
                            index cur_row = cur_col_values[ idx ];
                            if( ( column_type[ cur_row ] == GLOBAL ) ) {
                                is_active[ cur_col ] = true;
                            } else if( column_type[ cur_row ] == LOCAL_POSITIVE ) {
                                index next_col = lowest_one_lookup[ cur_row ];
                                if( next_col != cur_col && !finished[ cur_col ] ) {
                                    stack.push_back( std::make_pair( next_col, cur_col ) );
                                    pop_next = false;
                                }
                            }
                        }
                    }
                }
            }
        }

        template< typename Representation >
        void _global_column_simplification( const index col_idx
                                          , boundary_matrix< Representation >& boundary_matrix
                                          , const std::vector< index >& lowest_one_lookup
                                          , const std::vector< column_type >& column_type
                                          , const std::vector< char >& is_active
                                          , std::vector< index >& temp_col )
        {
            temp_col.clear();
            while( !boundary_matrix.is_empty( col_idx ) ) {
                index cur_row = boundary_matrix.get_max_index( col_idx );
                switch( column_type[ cur_row ] ) {
                case GLOBAL:
                    temp_col.push_back( cur_row );
                    boundary_matrix.remove_max( col_idx );
                    break;
                case LOCAL_NEGATIVE:
                    boundary_matrix.remove_max( col_idx );
                    break;
                case LOCAL_POSITIVE:
                    if( is_active[ lowest_one_lookup[ cur_row ] ] )
                        boundary_matrix.add_to( lowest_one_lookup[ cur_row ], col_idx );
                    else
                        boundary_matrix.remove_max( col_idx );
                    break;
                }
            }
            std::reverse( temp_col.begin(), temp_col.end() );
            boundary_matrix.set_col( col_idx, temp_col );
        }
    };
}
