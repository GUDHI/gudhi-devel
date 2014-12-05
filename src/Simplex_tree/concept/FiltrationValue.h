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

/** \brief Value type for a filtration function on a cell complex.
  *
  * A <EM>filtration</EM> of a cell complex (see FilteredComplex) is
  * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq 
  * f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. Ordering the simplices 
  * by increasing filtration values (breaking ties so as a simplex appears after 
  * its subsimplices of same filtration value) provides an indexing scheme 
  * (see IndexingTag).
  */
  struct FiltrationValue {
    /** \brief Operator < is a StrictWeakOrdering. */
    bool operator<(FiltrationValue f1, FiltrationValue f2);
  };
