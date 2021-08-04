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
  * for a type to implement an iterator over the insertions and deletions of 
  * cells in a zigzag filtration.
  * 
  * The dynamic nature of computations requires the iterator to encode methods, 
  * that are usually encoded by the data structure storing the cell complex.
  */
struct ZigzagFiltrationSimplexIterator {
/** A model of the concept ZigzagFilteredComplex. */
  typedef unspecified Complex;
/** Handle to specify a cell in a cell complex data structure. Must match 
 * ZigzagFilteredComplex::Simplex_handle.*/
  typedef unspecified      Simplex_handle;  
/** Star operator of the iterator, of value_type Simplex_handle. */
  Simplex_handle & operator*();
/** Comparison of iterators. */
  bool operator!=();
/** Increment the iterator to point to the next insertion or deletion of a cell in 
 * a zigzag filtration. */
  void operator++();
/** \brief Type for the value of the filtration function.
  *
  * Must be convertible from ZigzagFilteredComplex::Filtration_value. Must be 
  * comparable with <. */
  typedef unspecified      Filtration_value;
/** \brief Returns the filtration value of the event (insertion or removal of a 
 * cell in the zigzag filtraiton) currently pointed at by the iterator.
 *
 * \details Note that this may not be equal to the filtration value stored in a cell 
 * represented by a data structure for a cell complex, as, in a cell complex, the 
 * filtration value is usually the one of the insertion of the cell. The deletion of 
 * this same cell may in general happen at a different filtration value.
 */
  Filtration_value filtration();
/** \brief Returns an upper bound on the dimension of the cells iterated over.
  *
  * \details In particular used in order to not record persistence interval of maximal dimension.
  */
  int dim_max();
/** \brief Returns the direction (insertion or deletion) of the corresponding map.
  *
  * \details <CODE>true</CODE> for insertion, <CODE>false</CODE> for deletion.
  */
  bool arrow_direction();
};
