/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WITNESS_COMPLEX_ALL_FACES_IN_H_
#define WITNESS_COMPLEX_ALL_FACES_IN_H_

/* \brief Check if the facets of the k-dimensional simplex witnessed 
 *  by witness witness_id are already in the complex.
 *  inserted_vertex is the handle of the (k+1)-th vertex witnessed by witness_id
 */
template < typename SimplicialComplexForWitness,
           typename Simplex >
  bool all_faces_in(Simplex& simplex,
                    double* filtration_value,
                    SimplicialComplexForWitness& sc) {
    typedef typename SimplicialComplexForWitness::Simplex_handle Simplex_handle;

    if (simplex.size() == 1)
      return true;                  /* Add vertices unconditionally */

    Simplex facet;
    for (typename Simplex::iterator not_it = simplex.begin(); not_it != simplex.end(); ++not_it) {
        facet.clear();
        for (typename Simplex::iterator it = simplex.begin(); it != simplex.end(); ++it)
          if (it != not_it)
            facet.push_back(*it);
        Simplex_handle facet_sh = sc.find(facet);
        if (facet_sh == sc.null_simplex())
          return false;
        else if (sc.filtration(facet_sh) > *filtration_value)
          *filtration_value = sc.filtration(facet_sh);
      }
    return true;
  }

#endif  // WITNESS_COMPLEX_ALL_FACES_IN_H_
