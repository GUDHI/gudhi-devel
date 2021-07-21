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

/** \brief The concept ZigzagFiltrationSimplexIterator describes the requirements 
  * for a type to implement an iterator over the insertions and deletions of simplices 
  * in a a zigzag filtration.
  * 
  * The dynamic nature of computations requires the iterator to encode methods, 
  * that are usually encoded by the data structure storing the simplicial complex.
  */
struct ZigzagFiltrationSimplexIterator {
/** Handle to specify a simplex in a simplicial complex data structure. */
  typedef unspecified      Simplex_handle;  
/** Star operator of the iterator, of value_type Simplex_handle. */
  Simplex_handle & operator*();
/** Comparison of iterators. */
  bool operator!=();
/** Increment the iterator to point to the next insertion or deletion of a simplex in a zigzag filtration. */
  void operator++();
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns the filtration value of the simplex pointed at, at the time of insertion of the simplex.
  */
  Filtration_value filtration();
/** \brief Returns an upper bound on the dimension of the simplices iterated over.
  *
  * In particular used in in order to not record maximal dimension intervals 
  * in the zigzag persistence barcode.
  */
  int dim_max();
/** \brief Returns the direction (insertion or deletion) of the corresponding map.
  *
  * true for insertion, false for deletion.
  */
  bool arrow_direction();
};


/** \brief The concept ZigzagFilteredComplex describes the requirements 
  * for a type to implement a filtered simplicial complex, from which 
  * one can compute zigzag persistent homology via the class 
  * Zigzag_persistence< ZigzagFilteredComplex >
  */
struct ZigzagFilteredComplex
{
/** \brief Data stored for each simplex. 
  *
  * Must be an integer type, signed. */
  typedef unspecified      Simplex_key;
/** Handle to specify a simplex. */
  typedef unspecified      Simplex_handle;
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns a number assigned to the simplex pointed to by sh.
  * The key of a simplex is unique, and is equal to the index of insertion 
  * of the simplex in a zigzag filtration (hence at least 0).
  *
  * For example, in a filtration: insert a, insert b, remove a, insert c,
  * we can call key(simplex) when the simplex is present in the complex (i.e., there 
  * exists a valid Simplex_handle pointing to the simplex), and, when valid, we have 
  * key(a)==0, key(b)==1, key(c)==3.
  *
  * This is never called on null_simplex(). */
  Simplex_key              key      ( Simplex_handle sh );
/** \brief Iterator over all simplices of the complex 
  * in the order of the indexing scheme. Must be a model of
  * ZigzagFiltrationSimplexIterator
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Zigzag_filtration_simplex_iterator;
/** \brief Range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Zigzag_filtration_simplex_iterator.*/
typedef unspecified Zigzag_filtration_simplex_range;
/** \brief Returns a range over the simplices of the complex
  * in the order of the zigzag filtration. The iterator encodes whether this 
  * is an insertion of a deletion (see concept ZigzagFiltrationSimplexIterator).
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
Zigzag_filtration_simplex_range filtration_simplex_range();
/** \brief Returns the dimension of a simplex. */
  int dimension(Simplex_handle sh);
/** \brief Iterator over the simplices belonging to the
  * boundary of a simplex.
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Boundary_simplex_iterator;
/** \brief Range giving access to the simplices in the boundary of 
  * a simplex.
  *
  * .begin() and .end() return type Boundary_simplex_iterator.
  */
typedef unspecified Boundary_simplex_range;

/** \brief Returns a range giving access to all simplices of the 
  * boundary of a simplex, i.e.
  * the set of codimension 1 subsimplices of the Simplex.
  *
  * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
  * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
  * simplices of the boundary in the order: 
  * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
  *
  * We note that the alternate sum of the simplices given by the iterator
  * gives the chains corresponding to the boundary of the simplex.*/
Boundary_simplex_range boundary_simplex_range(Simplex_handle sh);

};









/** \brief The concept ZigzagFilteredMorseComplex is a refinement of 
  * ZigzagFiltered Complex with a discrete Morse theoretical layer. 
  * It describes the requirements 
  * for a type to implement a filtered simplicial complex, from which 
  * one can compute zigzag persistent homology via the class 
  * Zigzag_persistence< ZigzagFilteredComplex >, and a Morse matching.
  */
struct ZigzagFilteredMorseComplex
{
/** \brief Data stored for each simplex. 
  *
  * Must be an integer type, signed. */
  typedef unspecified      Simplex_key;
/** Handle to specify a simplex. */
  typedef unspecified      Simplex_handle;
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns a number assigned to the simplex pointed to by sh.
  * The key of a simplex is unique, and is equal to the index of insertion 
  * of the simplex in a zigzag filtration (hence at least 0).
  *
  * For example, in a filtration: insert a, insert b, remove a, insert c,
  * we can call key(simplex) when the simplex is present in the complex (i.e., there 
  * exists a valid Simplex_handle pointing to the simplex), and, when valid, we have 
  * key(a)==0, key(b)==1, key(c)==3.
  *
  * This is never called on null_simplex(). */
  Simplex_key              key      ( Simplex_handle sh );
/** \brief Iterator over all simplices of the complex 
  * in the order of the indexing scheme. Must be a model of
  * ZigzagFiltrationSimplexIterator
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Zigzag_filtration_simplex_iterator;
/** \brief Range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Zigzag_filtration_simplex_iterator.*/
typedef unspecified Zigzag_filtration_simplex_range;
/** \brief Returns a range over the simplices of the complex
  * in the order of the zigzag filtration. The iterator encodes whether this 
  * is an insertion of a deletion (see concept ZigzagFiltrationSimplexIterator).
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
Zigzag_filtration_simplex_range filtration_simplex_range();
/** \brief Returns the dimension of a simplex. */
  int dimension(Simplex_handle sh);
/** \brief Iterator over the simplices belonging to the
  * boundary of a simplex.
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Boundary_simplex_iterator;
/** \brief Range giving access to the simplices in the boundary of 
  * a simplex.
  *
  * .begin() and .end() return type Boundary_simplex_iterator.
  */
typedef unspecified Boundary_simplex_range;

/** \brief Returns a range giving access to all simplices of the 
  * boundary of a simplex, i.e.
  * the set of codimension 1 subsimplices of the Simplex.
  *
  * If the simplex is \f$[v_0, \cdots ,v_d]\f$, with canonical orientation
  * induced by \f$ v_0 < \cdots < v_d \f$, the iterator enumerates the 
  * simplices of the boundary in the order: 
  * \f$[v_0,\cdots,\widehat{v_i},\cdots,v_d]\f$ for \f$i\f$ from 0 to d
  *
  * We note that the alternate sum of the simplices given by the iterator
  * gives the chains corresponding to the boundary of the simplex.*/
Boundary_simplex_range boundary_simplex_range(Simplex_handle sh);
/** \brief In a Morse theoretical framework, returns true if the simplex is critical
  * with regards to a Morse matching, false if it is paired.
  */
  bool critical(Simplex_handle sh);
/** \brief If the simplex is paired in a Morse matching, returns a Simplex_handle 
  * to the simplex it is paired with. If sh is not paired (i.e., critical), 
  * morse_pair returns sh itself.
  */
  Simplex_handle morse_pair(Simplex_handle sh);
/** \brief Iterator over the simplices belonging to the
  * coboundary of a simplex (cofaces of codimension 1).
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
typedef unspecified Coboundary_simplex_iterator;
/** \brief Range giving access to the simplices in the coboundary of 
  * a simplex.
  *
  * .begin() and .end() return type Coboundary_simplex_iterator.
  */
typedef unspecified Coboundary_simplex_range;

/** \brief Returns a range giving access to all simplices of the 
  * coboundary of a simplex, i.e.
  * the set of codimension 1 cofaces of the Simplex.
  */
Coboundary_simplex_range coboundary_simplex_range(Simplex_handle sh);
/** Returns a Simplex_handle that is different from all simplex handles 
  * of the simplices. */
  // Simplex_handle           null_simplex();
/** \brief Returns a constant dummy number that is either negative,
  * or at least as large as `num_simplices()`.  Suggested value: -1.  */
  // Simplex_key              null_key ();
};
