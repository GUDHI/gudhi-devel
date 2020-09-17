/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief The concept FlagZigzagFilteredComplex describes the requirements 
  * for a type to implement a filtered cell complex, from which 
  * one can compute persistent homology via a model of the concept 
  * PersistentHomology. 
  */
struct FilteredComplex
{
/** \brief Handle to specify a simplex. */
  typedef unspecified      Simplex_handle;
/** \brief Handle to specify a vertex. */
  typedef unspecified      Vertex_handle;
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;

/** \brief Specifies the nature of the indexing scheme. 
  * 
  * is model of IndexingTag. */
  // typedef unspecified      Indexing_tag;

/** \brief Data stored for each simplex. 
  *
  * Must be an integer type. */
  typedef unspecified      Simplex_key;
/** \brief Returns a constant dummy number that is either negative,
  * or at least as large as `num_simplices()`.  Suggested value: -1.  */
  // Simplex_key              null_key ();
/** \brief Returns the number stored for a simplex by `assign_key`.
  *
  * This is never called on null_simplex(). */
  // Simplex_key              key      ( Simplex_handle sh );
/** \brief Store a number for a simplex, which can later be retrieved with `key(sh)`.
  *
  * This is never called on null_simplex(). */
  void                     assign_key(Simplex_handle sh, Simplex_key n);

  void remove_maximal_simplex(Simplex_handle sh);
  void flag_add_edge(Vertex_handle u, Vertex_handle v, Filtration_value fil, 
    int dim_max, std::vector< Simplex_handle > & zz_filtration); 
  void flag_lazy_remove_edge(Vertex_handle u, Vertex_handle v, 
    std::vector< Simplex_handle > & zz_filtration );
  void flag_lazy_empty_complex(std::vector< Simplex_handle > & zz_filtration);

/** @} */


};
