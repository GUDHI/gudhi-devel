/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef PERSISTENT_COHOMOLOGY_PERSISTENT_COHOMOLOGY_COLUMN_H_
#define PERSISTENT_COHOMOLOGY_PERSISTENT_COHOMOLOGY_COLUMN_H_

#include <boost/intrusive/set.hpp>
#include <boost/intrusive/list.hpp>

#include <list>

namespace Gudhi {

namespace persistent_cohomology {

template<typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_column;

struct cam_h_tag;
// for horizontal traversal in the CAM
struct cam_v_tag;
// for vertical traversal in the CAM

typedef boost::intrusive::list_base_hook<boost::intrusive::tag<cam_h_tag>,
    boost::intrusive::link_mode<boost::intrusive::auto_unlink>  // allows .unlink()
> base_hook_cam_h;

typedef boost::intrusive::list_base_hook<boost::intrusive::tag<cam_v_tag>,
    boost::intrusive::link_mode<boost::intrusive::normal_link>  // faster hook, less safe
> base_hook_cam_v;

/** \internal
 * \brief
 *
 */
template<typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_cell : public base_hook_cam_h,
    public base_hook_cam_v {
 public:
  template<class T1, class T2> friend class Persistent_cohomology;
  friend class Persistent_cohomology_column<SimplexKey, ArithmeticElement>;

  typedef Persistent_cohomology_column<SimplexKey, ArithmeticElement> Column;

  Persistent_cohomology_cell(SimplexKey key, ArithmeticElement x,
                             Column * self_col)
      : key_(key),
        coefficient_(x),
        self_col_(self_col) {
  }

  SimplexKey key_;
  ArithmeticElement coefficient_;
  Column * self_col_;
};

/* 
 * \brief Sparse column for the Compressed Annotation Matrix.
 *
 * The non-zero coefficients of the column are stored in a
 * boost::intrusive::list. Contains a hook to be stored in a
 * boost::intrusive::set.
 *
 * Movable but not Copyable.
 */
template<typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_column : public boost::intrusive::set_base_hook<
    boost::intrusive::link_mode<boost::intrusive::normal_link> > {
  template<class T1, class T2> friend class Persistent_cohomology;

 public:
  typedef Persistent_cohomology_cell<SimplexKey, ArithmeticElement> Cell;
  typedef boost::intrusive::list<Cell,
      boost::intrusive::constant_time_size<false>,
      boost::intrusive::base_hook<base_hook_cam_v> > Col_type;

  /** \brief Creates an empty column.*/
  explicit Persistent_cohomology_column(SimplexKey key)
      : col_(),
        class_key_(key) {}

  /** \brief Returns true iff the column is null.*/
  bool is_null() const {
    return col_.empty();
  }
  /** \brief Returns the key of the representative simplex of
   * the set of simplices having this column as annotation vector
   * in the compressed annotation matrix.*/
  SimplexKey class_key() const {
    return class_key_;
  }

  /** \brief Lexicographic comparison of two columns.*/
  friend bool operator<(const Persistent_cohomology_column& c1,
                        const Persistent_cohomology_column& c2) {
    typename Col_type::const_iterator it1 = c1.col_.begin();
    typename Col_type::const_iterator it2 = c2.col_.begin();
    while (it1 != c1.col_.end() && it2 != c2.col_.end()) {
      if (it1->key_ == it2->key_) {
        if (it1->coefficient_ == it2->coefficient_) {
          ++it1;
          ++it2;
        } else {
          return it1->coefficient_ < it2->coefficient_;
        }
      } else {
        return it1->key_ < it2->key_;
      }
    }
    return (it2 != c2.col_.end());
  }

  Col_type col_;
  SimplexKey class_key_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_PERSISTENT_COHOMOLOGY_COLUMN_H_
