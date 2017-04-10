/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2017  INRIA
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

#ifndef NULL_OUTPUT_ITERATOR_H_
#define NULL_OUTPUT_ITERATOR_H_

#include <iterator>

namespace Gudhi {

/** An output iterator that ignores whatever it is given. */
struct Null_output_iterator {
  typedef std::output_iterator_tag iterator_category;
  typedef void                     value_type;
  typedef void                     difference_type;
  typedef void                     pointer;
  typedef void                     reference;

  Null_output_iterator& operator++(){return *this;}
  Null_output_iterator operator++(int){return *this;}
  struct proxy {
    template<class T>
    proxy& operator=(T&&){return *this;}
  };
  proxy operator*()const{return {};}
};
}  // namespace Gudhi

#endif  // NULL_OUTPUT_ITERATOR_H_
