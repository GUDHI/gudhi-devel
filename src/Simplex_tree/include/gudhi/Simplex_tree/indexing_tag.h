/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SIMPLEX_TREE_INDEXING_TAG_H_
#define SIMPLEX_TREE_INDEXING_TAG_H_

namespace Gudhi {

/** \brief Tag for a linear ordering of simplices. 
 *
 * \implements IndexingTag
 */
struct linear_indexing_tag {
};

/* \brief Tag for a zigzag ordering of simplices. */
//  struct zigzag_indexing_tag {};
}  // namespace Gudhi

#endif  // SIMPLEX_TREE_INDEXING_TAG_H_
