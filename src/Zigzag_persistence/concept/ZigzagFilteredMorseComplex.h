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

/** \brief The concept ZigzagFilteredMorseComplex is the union of requirements for 
 * a model of the concepts <CODE>ZigzagFilteredComplex</CODE> and <CODE>FilteredMorseComplex</CODE>.
  */
struct ZigzagFilteredMorseComplex : public ZigzagFilteredComplex, public FilteredMorseComplex
{
/** \brief Iterator over all cells of the complex 
  * in the order of the indexing scheme. 
  * 
  * Additional requirements for <CODE>ZigzagFilteredMorseComplex</CODE>:
  * Must be a model of
  * <CODE>ZigzagMorseFiltrationSimplexIterator</CODE>
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
    * (see concept <CODE>ZigzagMorseFiltrationSimplexIterator</CODE>).
    */
    Filtration_simplex_range filtration_simplex_range();    
};
