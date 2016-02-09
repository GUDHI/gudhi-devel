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

#include "persistence_pairs.h"
#include "boundary_matrix.h"
#include "helpers/dualize.h"
#include "algorithms/twist_reduction.h"

namespace phat {

    template< typename ReductionAlgorithm, typename Representation >
    void compute_persistence_pairs( persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        ReductionAlgorithm reduce;
        reduce( boundary_matrix );
        pairs.clear();
        for( index idx = 0; idx < boundary_matrix.get_num_cols(); idx++ ) {
            if( !boundary_matrix.is_empty( idx ) ) {
                index birth = boundary_matrix.get_max_index( idx );
                index death = idx;
                pairs.append_pair( birth, death );
            }
        }
    }

    template< typename ReductionAlgorithm, typename Representation >
    void compute_persistence_pairs_dualized( persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        ReductionAlgorithm reduce;
        const index nr_columns = boundary_matrix.get_num_cols();
        dualize( boundary_matrix );
        reduce( boundary_matrix );
        pairs.clear();
        for( index idx = 0; idx < nr_columns; idx++ ) {
            if( !boundary_matrix.is_empty( idx ) ) {
                index death = nr_columns - 1 - boundary_matrix.get_max_index( idx );
                index birth = nr_columns - 1 - idx;
                pairs.append_pair( birth, death );
            }
        }
    }

    template< typename Representation >
    void compute_persistence_pairs( persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        phat::compute_persistence_pairs< twist_reduction >( pairs, boundary_matrix );
    }


    template< typename Representation >
    void compute_persistence_pairs_dualized( persistence_pairs& pairs, boundary_matrix< Representation >& boundary_matrix ) {
        compute_persistence_pairs_dualized< twist_reduction >( pairs, boundary_matrix );
    }

}
