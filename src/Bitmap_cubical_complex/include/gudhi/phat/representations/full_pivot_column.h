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
        std::priority_queue< index > m_history;
        std::vector< char > m_isInHistory;
        std::vector< char > m_data;

    public:
        void init( const index total_size ) {
            m_data.resize( total_size, false );
            m_isInHistory.resize( total_size, false );
        }

        void add_column( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ ) {
                add_index( col[ idx ] );
            }
        }
        void add_index( const index idx ) {
            if( !m_isInHistory[ idx ] ) {
                m_history.push( idx );
                m_isInHistory[ idx ] = true; 
            }

            m_data[ idx ] = !m_data[ idx ];
        }

        index max_index() {
            while( m_history.size() > 0 ) {
                index topIndex = m_history.top();
                if( m_data[ topIndex ] ) {
                    return topIndex;
                } else {
                    m_history.pop();
                    m_isInHistory[ topIndex ] = false;
                }
            }
    
            return -1;
        }

        void get_column_and_clear( column& col ) {
            col.clear();
            while( !empty() ) {
                col.push_back( max_index() );
                add_index( max_index() );
            }
            std::reverse( col.begin(), col.end() );
        }

        bool empty() {
            return (max_index() == -1);   
        }
    };

    typedef abstract_pivot_column< full_column > full_pivot_column;
}
