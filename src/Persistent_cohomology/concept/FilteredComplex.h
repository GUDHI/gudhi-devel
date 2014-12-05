 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
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

/** \brief The concept FilteredComplex describes the requirements 
  * for a type to implement a filtered cell complex, from which 
  * one can compute persistent homology via a model of the concept 
  * PersistentHomology. 
  */
struct FilteredComplex
{
/** Handle to specify a simplex. */
  typedef unspecified      Simplex_handle;
/** \brief Key associated to each simplex. 
  *
  * Must be a signed integer type. */
  typedef unspecified      Simplex_key;
/** \brief Type for the value of the filtration function.
  *
  * Must be comparable with <. */
  typedef unspecified      Filtration_value;

/** \brief Specifies the nature of the indexing scheme. 
  * 
  * is model of IndexingTag. */
  typedef unspecified      Indexing_tag;

/** Returns a Simplex_hanlde that is different from all simplex handles 
  * of the simplices. */
  Simplex_handle           null_simplex();
/** \brief Returns the number of simplices in the complex.
  *
  * Does not count the empty simplex. */
  size_t                   num_simplices();
/** \brief Returns the dimension of a simplex. */
  int                      dimension(Simplex_handle sh);
/** \brief Returns the filtration value of a simplex. 
  * 
  * If sh is null_simplex(), returns the maximal value of the
  * filtration function on the complex. */
  Filtration_value         filtration(Simplex_handle sh);

/** \brief Returns a key that is different from the keys associated 
  * to the simplices. */
  Simplex_key              null_key ();
/** \brief Returns the key associated to a simplex. */  
  Simplex_key              key      ( Simplex_handle sh );
/** \brief Returns the simplex associated to a key.
  *
  * If key is different from null_key(), there must be a unique 
  * simplex having this key. */  
  Simplex_handle           simplex  ( Simplex_key key );
/** \brief Assign a key to a simplex. */
  void                     assign_key(Simplex_handle sh, Simplex_key key);
 
/** \brief Iterator on the simplices belonging to the
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

/** \brief Iterator over all simplices of the complex 
  * in the order of the indexing scheme.
  *
  * 'value_type' must be 'Simplex_handle'.
  */
typedef unspecified Filtration_simplex_iterator;
/** \brief Range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
typedef unspecified Filtration_simplex_range;
/** \brief Returns a range over the simplices of the complex
  * in the order of the filtration.
  *
  * .begin() and .end() return type Filtration_simplex_iterator.*/
Filtration_simplex_range filtration_simplex_range();


/* \brief Iterator over the simplices of the complex,
  * in an arbitrary order.
  *
  * 'value_type' must be 'Simplex_handle'.*/
//typedef unspecified Complex_simplex_iterator;
//typedef unspecified Complex_simplex_range;

/*
* Returns a range over all the simplices of a
* complex.
*/
//Complex_simplex_range complex_simplex_range();

/*************************************************/     /**
* @details Compares the filtration values of simplices s and t
*
* @return -1 if s comes before t in the filtration, +1
* if it comes later, and 0 if they come at the same time
*
* @note OPTIONAL
* @todo use an enum? Just a bool?
*/
//int is_before_in_filtration(Simplex_handle s, Simplex_handle t);
/*************************************************/     

};
