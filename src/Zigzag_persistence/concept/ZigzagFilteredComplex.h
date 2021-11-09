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
/** \brief The concept ZigzagFilteredComplex describes the requirements 
  * for a type to implement a filtered cell complex, from which 
  * one can compute zigzag persistent homology via the class 
  * <CODE>Zigzag_persistence< ZigzagFilteredComplex ></CODE>.
  * 
  * \details It is a refinement of the concept <CODE>FilteredComplex</CODE>
  */
struct ZigzagFilteredComplex : public FilteredComplex 
{
/** \brief Data stored for each simplex. 
  *
  * Must be an integer type. 
  * 
  * Additional requirements for <CODE>ZigzagFilteredComplex</CODE>: Must be signed.*/
  typedef unspecified      Simplex_key;
/** \brief Returns the number stored for a simplex by `assign_key`.
  *
  * \details This is never called on null_simplex(). 
  * 
  * Additional requirements for <CODE>ZigzagFilteredComplex</CODE>: 
  * The key of a simplex is unique, and is equal to the index of insertion 
  * of the simplex in a zigzag filtration (hence at least 0).
  *
  * For example, in a filtration: insert a, insert b, remove a, insert c,
  * we can call <CODE>key(sh)</CODE> when the cell is in the complex currently 
  * maintained (i.e., there exists a valid <CODE>Simplex_handle</CODE> pointing 
  * to the cell), and, when valid, we have 
  * <CODE>key(a)==i</CODE>, <CODE>key(b)==i+1</CODE>, </CODE>key(c)==i+3<CODE> for 
  * some positive integer \f$i\f$.
  */
  Simplex_key              key      ( Simplex_handle sh );
/** \brief Iterator over all cells of the complex 
  * in the order of the indexing scheme. 
  * 
  * Additional requirements for <CODE>ZigzagFilteredComplex</CODE>:
  * Must be a model of
  * <CODE>ZigzagFiltrationSimplexIterator</CODE>
  *
  * <CODE>value_type</CODE> must be 'Simplex_handle'.
  */
  typedef unspecified Filtration_simplex_iterator;
  /** \brief Range over the simplices of the complex
    * in the order of the filtration.
    *
    * \details <CODE>.begin()</CODE> and <CODE>.end()</CODE> return type 
    * <CODE>Filtration_simplex_iterator</CODE> as defined above.*/
  typedef unspecified Filtration_simplex_range;
  /** \brief Returns a range over the simplices of the complex
    * in the order of the zigzag filtration. The iterator encodes whether this 
    * is an insertion of a deletion 
    * (see concept <CODE>ZigzagFiltrationSimplexIterator</CODE>).
    */
    Filtration_simplex_range filtration_simplex_range();
};
