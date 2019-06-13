/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Marc Glisse
 *
 *    Copyright (C) 2017 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
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

  Null_output_iterator& operator++() {return *this;}
  Null_output_iterator operator++(int) {return *this;}
  struct proxy {
    template<class T>
    proxy& operator=(T&&){return *this;}
  };
  proxy operator*()const{return {};}
};
}  // namespace Gudhi

#endif  // NULL_OUTPUT_ITERATOR_H_
