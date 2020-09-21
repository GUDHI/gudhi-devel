/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONCEPT_COXETER_TRIANGULATION_SIMPLEX_IN_COXETER_TRIANGULATION_H_
#define CONCEPT_COXETER_TRIANGULATION_SIMPLEX_IN_COXETER_TRIANGULATION_H_

namespace Gudhi {

namespace coxeter_triangulation {

/** \brief The concept SimplexInCoxeterTriangulation describes the requirements 
 * for a type to implement a representation of simplices in Freudenthal_triangulation
 * or in Coxeter_triangulation.
 */
struct SimplexInCoxeterTriangulation {

  /** \brief Type of the vertex. */
  typedef Vertex_ Vertex;
  
  /** \brief Type of the ordered partition. */
  typedef Ordered_set_partition_ OrderedSetPartition;

  /** \brief Dimension of the simplex. */
  unsigned dimension() const;

  /** \brief Type of a range of vertices, each of type Vertex. */
  typedef Vertex_range;
  
  /** \brief Returns a range of vertices of the simplex.
   */
  Vertex_range vertex_range() const;

  /** \brief Type of a range of faces, each of type that 
   *  is a model of the concept SimplexInCoxeterTriangulation. 
   */
  typedef Face_range;

  /** \brief Returns a range of permutahedral representations of k-dimensional faces
   *   of the simplex for some given integer parameter 'k'.
   */
  Face_range face_range(std::size_t k) const;

  /** \brief Returns a range of permutahedral representations of facets of the simplex.
   * The dimension of the simplex must be strictly positive.
   */
  Face_range facet_range() const;

  /** \brief Type of a range of cofaces, each of type that 
   *  is a model of the concept SimplexInCoxeterTriangulation. 
   */
  typedef Coface_range;

  /** \brief Returns a range of permutahedral representations of k-dimensional cofaces
   *   of the simplex for some given integer parameter 'k'.
   */
  Coface_range coface_range(std::size_t k) const;

  /** \brief Returns a range of permutahedral representations of cofacets of the simplex.
   * The dimension of the simplex must be strictly different from the ambient dimension.
   */
  Coface_range cofacet_range() const;

  /** \brief Returns true, if the simplex is a face of other simplex. */
  bool is_face_of(const Permutahedral_representation& other) const;
  
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi


#endif
