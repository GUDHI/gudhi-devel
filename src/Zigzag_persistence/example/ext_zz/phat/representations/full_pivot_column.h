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
#include <phat/representations/abstract_pivot_column.h>

namespace phat {
    class full_column {

    protected:
        std::priority_queue< index > history;
        std::vector< char > is_in_history;
        std::vector< char > col_bit_field;

    public:
        void init( const index total_size ) {
            col_bit_field.resize( total_size, false );
            is_in_history.resize( total_size, false );
        }

        void add_col( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ ) {
                add_index( col[ idx ] );
            }
        }

        void add_index( const index idx ) {
            if( !is_in_history[ idx ] ) {
                history.push( idx );
                is_in_history[ idx ] = true; 
            }

            col_bit_field[ idx ] = !col_bit_field[ idx ];
        }

        index get_max_index() {
            while( history.size() > 0 ) {
                index topIndex = history.top();
                if( col_bit_field[ topIndex ] ) {
                    return topIndex;
                } else {
                    history.pop();
                    is_in_history[ topIndex ] = false;
                }
            }
    
            return -1;
        }

        void get_col_and_clear( column& col ) {
            while( !is_empty() ) {
                col.push_back( get_max_index() );
                add_index( get_max_index() );
            }
            std::reverse( col.begin(), col.end() );
        }

        bool is_empty() {
            return (get_max_index() == -1);   
        }

		void clear() {
			while( !is_empty() )
				add_index( get_max_index() );
		}

		void remove_max() {
            add_index( get_max_index() );
        }

        void set_col( const column& col  ) {
            clear();
            add_col( col );
        }

        void get_col( column& col  ) {
            get_col_and_clear( col );
            add_col( col );
        }
    };

    typedef abstract_pivot_column< full_column > full_pivot_column;
}
