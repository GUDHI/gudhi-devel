/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM author: description of the modification.
 *         
 */
/** \brief The concept DynamicFilteredFlagComplex is a refinement of 
  * <CODE>FilteredComplex</CODE>, which describes the requirements 
  * for a type to implement a filtered flag complex, that can be dynamically 
  * updated by the additions and removals of edges and their cofaces.
  */
struct DynamicFilteredFlagComplex : public FilteredComplex
{
/** \brief Store a number for a simplex, which can later be retrieved with `key(sh)`.
  *
  * This is never called on null_simplex(). */
  void                     assign_key(Simplex_handle sh, Simplex_key n);

  /** \brief Add an edge to the flag complex, together with all of its cofaces up 
   * to an input dimension.
   * 
   * \details The edge is given by the handles of its two endpoints ; in case 
   * <CODE>u == v</CODE>, this corresponds to the insertion of the vertex of 
   * <CODE>Vertex_handle u</CODE>.
   *
   * All newly inserted faces, including the input edge/vertex, are recorded in an 
   * <CODE>std::vector</CODE> sorted in increasing filtration order (subfaces 
   * before cofaces).
   * 
   * All newly inserted faces are given filtration value <CODE>fil</CODE>.
   * 
   * The maximal dimension of expansion if given by <CODE>dim_max</CODE>.
   */
  void flag_add_edge(Vertex_handle u, Vertex_handle v, Filtration_value fil, 
    int dim_max, std::vector< Simplex_handle > & zz_filtration); 

/** \brief Removes a maximal cell from the complex.*/ 
  void remove_maximal_simplex(Simplex_handle sh);
/** \brief Records, but does not remove, all cofaces of an input edge \f$(u,v)\f$ 
 * from the complex.
 *
 * \details In preparation of the removal of an edge. All to be removed faces, 
 * including the input edge/vertex, are recorded in an 
 * <CODE>std::vector</CODE>. The cells are not sorted in <CODE>zz_filtration<\CODE> 
 * after computation.
 *
 * @param[in] u,v Vertex handles representing the edge (or vertex if u=v) to 
 *                be removed.
 * @param[in] zz_filtration Cofaces of the edge are pushed at the back of this 
 *                          vector. It may not be empty at the begining of 
 *                          computation, and may already contain cofaces of the edge.
 */ 
  void flag_lazy_remove_edge(Vertex_handle u, Vertex_handle v, 
    std::vector< Simplex_handle > & zz_filtration );
/** \brief Records, but does not remove, all all cells in the complex.
 *
 * \details In preparation of the removal of all faces. They are recorded in an 
 * <CODE>std::vector</CODE> sorted in decreasing filtration order (cofaces 
 * before subfaces).
 */ 
  void flag_lazy_empty_complex(std::vector< Simplex_handle > & zz_filtration);
};
