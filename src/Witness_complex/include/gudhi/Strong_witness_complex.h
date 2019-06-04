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

#ifndef STRONG_WITNESS_COMPLEX_H_
#define STRONG_WITNESS_COMPLEX_H_

#include <gudhi/Active_witness/Active_witness.h>

#include <utility>
#include <vector>
#include <list>
#include <limits>

namespace Gudhi {

namespace witness_complex {

  /**
   * \private
 * \class Strong_witness_complex
 * \brief Constructs strong witness complex for a given table of nearest landmarks with respect to witnesses.
 * \ingroup witness_complex
 *
 * \tparam Nearest_landmark_table_ needs to be a range (one entry per witness)
 *         of sorted ranges of pairs of nearest landmarks and distances.
 *         The class Nearest_landmark_table_::value_type must be a copiable range.
 *         The range of pairs must admit a member type 'iterator'. The dereference type 
 *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'
 *         where the first element is the index of the landmark, and the second its
 *         (squared) distance to the witness.
 */
template< class Nearest_landmark_table_ >
class Strong_witness_complex {
 private:
  typedef typename Nearest_landmark_table_::value_type               Nearest_landmark_range;
  typedef std::size_t                                                Witness_id;
  typedef std::size_t                                                Landmark_id;
  typedef std::pair<Landmark_id, double>                             Id_distance_pair;
  typedef Active_witness<Id_distance_pair, Nearest_landmark_range>   ActiveWitness;
  typedef std::list< ActiveWitness >                                 ActiveWitnessList;
  typedef std::vector< Landmark_id >                                 typeVectorVertex;
  typedef std::vector<Nearest_landmark_range>                        Nearest_landmark_table_internal;
  typedef Landmark_id Vertex_handle;

 protected:
  Nearest_landmark_table_internal              nearest_landmark_table_;

 public:
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* @name Constructor
   */

  //@{

  Strong_witness_complex() {
  }

  /**
   *  \brief Initializes member variables before constructing simplicial complex.
   *  \details Records nearest landmark table.
   *  @param[in] nearest_landmark_table needs to be a range of a range of pairs of nearest landmarks and distances.
   *         The class Nearest_landmark_table_::value_type must be a copiable range.
   *         The range of pairs must admit a member type 'iterator'. The dereference type 
   *         of the pair range iterator needs to be 'std::pair<std::size_t, double>'.
   */
  Strong_witness_complex(Nearest_landmark_table_ const & nearest_landmark_table)
    : nearest_landmark_table_(std::begin(nearest_landmark_table), std::end(nearest_landmark_table)) {
  }

  /** \brief Outputs the strong witness complex of relaxation 'max_alpha_square' 
   *         in a simplicial complex data structure.
   *  \details The function returns true if the construction is successful and false otherwise.
   *  @param[out] complex Simplicial complex data structure, which is a model of
   *              SimplicialComplexForWitness concept.
   *  @param[in] max_alpha_square Maximal squared relaxation parameter.
   *  @param[in] limit_dimension Represents the maximal dimension of the simplicial complex
   *         (default value = no limit).
   */
  template < typename SimplicialComplexForWitness >
  bool create_complex(SimplicialComplexForWitness& complex,
                      double  max_alpha_square,
                      Landmark_id limit_dimension = std::numeric_limits<Landmark_id>::max()) const {
    Landmark_id complex_dim = 0;
    if (complex.num_vertices() > 0) {
      std::cerr << "Strong witness complex cannot create complex - complex is not empty.\n";
      return false;
    }
    if (max_alpha_square < 0) {
      std::cerr << "Strong witness complex cannot create complex - squared relaxation parameter must be "
                << "non-negative.\n";
      return false;
    }
    for (auto w : nearest_landmark_table_) {
      ActiveWitness aw(w);
      typeVectorVertex simplex;
      typename ActiveWitness::iterator aw_it = aw.begin();
      float lim_dist2 = aw.begin()->second + max_alpha_square;
      while ((Landmark_id)simplex.size() <= limit_dimension && aw_it != aw.end() && aw_it->second < lim_dist2) {
        simplex.push_back(aw_it->first);
        complex.insert_simplex_and_subfaces(simplex, aw_it->second - aw.begin()->second);
        aw_it++;
      }
      // continue inserting limD-faces of the following simplices
      typeVectorVertex& vertices = simplex;  // 'simplex' now will be called vertices
      while (aw_it != aw.end() && aw_it->second < lim_dist2) {
        typeVectorVertex facet = {};
        add_all_faces_of_dimension(limit_dimension, vertices, vertices.begin(), aw_it,
                                   aw_it->second - aw.begin()->second, facet, complex);
        vertices.push_back(aw_it->first);
        aw_it++;
      }
      if ((Landmark_id)simplex.size() - 1 > complex_dim)
        complex_dim = simplex.size() - 1;
    }
    return true;
  }

  //@}

 private:
    /* \brief Adds recursively all the faces of a certain dimension dim-1 witnessed by the same witness.
     * Iterator is needed to know until how far we can take landmarks to form simplexes.
     * simplex is the prefix of the simplexes to insert.
     * The landmark pointed by aw_it is added to all formed simplices.
     */
  template < typename SimplicialComplexForWitness >
  void add_all_faces_of_dimension(Landmark_id dim,
                                  typeVectorVertex& vertices,
                                  typename typeVectorVertex::iterator curr_it,
                                  typename ActiveWitness::iterator aw_it,
                                  double filtration_value,
                                  typeVectorVertex& simplex,
                                  SimplicialComplexForWitness& sc) const {
    if (dim > 0) {
      while (curr_it != vertices.end()) {
        simplex.push_back(*curr_it);
        ++curr_it;
        add_all_faces_of_dimension(dim-1,
                                   vertices,
                                   curr_it,
                                   aw_it,
                                   filtration_value,
                                   simplex,
                                   sc);
        simplex.pop_back();
        add_all_faces_of_dimension(dim,
                                   vertices,
                                   curr_it,
                                   aw_it,
                                   filtration_value,
                                   simplex,
                                   sc);
      }
    } else if (dim == 0) {
      simplex.push_back(aw_it->first);
      sc.insert_simplex_and_subfaces(simplex, filtration_value);
      simplex.pop_back();
    }
  }
};

}  // namespace witness_complex

}  // namespace Gudhi

#endif  // STRONG_WITNESS_COMPLEX_H_
