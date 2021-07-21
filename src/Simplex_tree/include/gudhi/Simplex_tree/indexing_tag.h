/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2021/07 Clément Maria: Add zigzag_indexing_tag
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

/** \brief Tag for a zigzag ordering of simplices. 
 *
 * \implements IndexingTag
 */ 
 struct zigzag_indexing_tag {
};
 
} // namespace Gudhi

#endif  // SIMPLEX_TREE_INDEXING_TAG_H_
