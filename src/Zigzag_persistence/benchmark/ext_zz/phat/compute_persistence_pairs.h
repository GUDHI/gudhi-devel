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

#include <phat/persistence_pairs.h>
#include <phat/boundary_matrix.h>
#include <phat/helpers/dualize.h>
#include <phat/algorithms/twist_reduction.h>

namespace phat {
    // Extracts persistence pairs in separate dimensions from a reduced
    // boundary matrix representing ``double`` filtration. The pairs
    // give persistent relative homology of the pair of filtrations.
    // TODO: Use it with standard reduction algorithm (no template option).
    template< typename ReductionAlgorithm, typename Representation >
    void compute_relative_persistence_pairs(std::vector<persistence_pairs>& pairs, boundary_matrix<Representation>& boundary_matrix, const std::map<int, int>& L) {
        ReductionAlgorithm reduce;
        reduce(boundary_matrix);
        std::map<int, bool> free;
        std::map<int, int> invL;
        for (std::map<int, int>::const_iterator it = L.begin(); it != L.end(); ++it) { invL[it->second] = it->first; }
        for (std::vector<persistence_pairs>::iterator it = pairs.begin(); it != pairs.end(); ++it) { it->clear(); }
        for (index idx = 0; idx < boundary_matrix.get_num_cols(); ++idx) {
            int dimension = boundary_matrix.get_dim(idx);
            if (L.find(idx) != L.end()) { ++dimension; }
            free[idx] = true;
            if (!boundary_matrix.is_empty(idx)) {
                index birth = boundary_matrix.get_max_index(idx);
                index death = idx;
                pairs[dimension-1].append_pair(birth, death);
                free[birth] = false;
                free[death] = false;
            } else {
                // This is an L-simplex and a (dimension+1)-dimensional cycle
                if (L.find(idx) != L.end()) {
                    assert(dimension < pairs.size());
                    pairs[dimension].append_pair(idx, -1);
                }
            }
        }
        for (std::map<int, bool>::iterator it = free.begin(); it != free.end(); ++it) {
            if (it->second) {
                int dimension = boundary_matrix.get_dim(it->first);
                if (invL.find(it->first) == invL.end() && L.find(it->first) == L.end()) {
                    assert(dimension < pairs.size());
                    pairs[dimension].append_pair(it->first, -1);
                }
            }
        }
    }

    // Extracts persistence pairs in separate dimensions; expects a d-dimensional vector of persistent_pairs
    template< typename ReductionAlgorithm, typename Representation >
    void compute_persistence_pairs(std::vector<persistence_pairs>& pairs, boundary_matrix<Representation>& boundary_matrix) {
        ReductionAlgorithm reduce;
        reduce(boundary_matrix);
        std::map<int, bool> free;
        for (std::vector<persistence_pairs>::iterator it = pairs.begin(); it != pairs.end(); ++it) { it->clear(); }
        for (index idx = 0; idx < boundary_matrix.get_num_cols(); ++idx) {
            int dimension = boundary_matrix.get_dim(idx);
            free[idx] = true;
            if (!boundary_matrix.is_empty(idx)) {
                index birth = boundary_matrix.get_max_index(idx);
                index death = idx;
                pairs[dimension-1].append_pair(birth, death);
                // Cannot be of the form (a, infinity)
                free[birth] = false;
                free[death] = false;
            }
        }
        for (std::map<int, bool>::iterator it = free.begin(); it != free.end(); ++it) {
            if (it->second) {
                int dimension = boundary_matrix.get_dim(it->first);
                pairs[dimension].append_pair(it->first, -1);
            }
        }
    }

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

        dualize( boundary_matrix );
        compute_persistence_pairs< ReductionAlgorithm >( pairs, boundary_matrix );
        dualize_persistence_pairs( pairs, boundary_matrix.get_num_cols() );
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
