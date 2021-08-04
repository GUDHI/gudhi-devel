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
/** \brief The concept ZigzagMorseFiltrationSimplexIterator is a refinement of 
  * ZigzagFiltrationSimplexIterator, which describes the requirements 
  * for a type to implement an iterator over the insertions and deletions of 
  * cells in a zigzag filtration, where the complex stores a Morse matching.
  */
struct ZigzagMorseFiltrationSimplexIterator : public ZigzagFiltrationSimplexIterator
{
/** \brief Indicates if a Morse pair \f$(\tau,\sigma)\f$ of non-critical cells is 
 * removed at once in a zigzag filtration.
 *
 * \detail Returns <CODE>true</CODE> iff the iterator is pointing to the removal of 
 * a cell \f$\sigma\f$, that is: 
 * - non-critical, and paired with a cell \f$\tau\f$ of the complex, \f$\tau \prec 
 * \sigma\f$, forming the Morse pair \f$(\tau,\sigma)\f$, and 
 * - the next element, pointed by the iterator incremented once, is not \f$\tau\f$.
 * 
 * Otherwise, returns <CODE>false</CODE>.
 */
  bool break_morse_pair();
};
