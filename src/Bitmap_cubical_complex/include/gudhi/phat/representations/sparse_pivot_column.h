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
    class sparse_column {

    protected:
        std::set< index > m_data;

    public:
        void init( const index total_size ) {
            m_data.clear(); 
        }

        void add_column( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ )
                add_index( col[ idx ] );
        }

        void add_index( const index idx ) {
            std::pair< std::set< index >::iterator, bool > result = m_data.insert( idx );
            if( result.second == false )
                m_data.erase( result.first );
        }

        index max_index() {
            return m_data.empty() ? -1 : *m_data.rbegin();
        }

        void get_column_and_clear( column& col ) {
            col.clear();
            col.assign( m_data.begin(), m_data.end() );
            m_data.clear();
        }

        bool empty() {
            return m_data.empty();
        }
    };

    typedef abstract_pivot_column< sparse_column > sparse_pivot_column;
}
