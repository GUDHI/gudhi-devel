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

// should ideally be equal to the cache line size of the CPU
#define PHAT_TLS_SPACING_FACTOR 64

// ThreadLocalStorage with some spacing to avoid "false sharing" (see wikipedia)
template< typename T > 
class thread_local_storage
{
public:

    thread_local_storage() : per_thread_storage( omp_get_max_threads() * PHAT_TLS_SPACING_FACTOR ) {};

    T& operator()() {
        return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator()() const {
        return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    T& operator[]( int tid ) {
        return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator[]( int tid ) const {
        return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

protected:
    std::vector< T > per_thread_storage; 
};
