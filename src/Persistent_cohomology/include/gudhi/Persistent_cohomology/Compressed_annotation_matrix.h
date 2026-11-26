/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria, Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_PERSISTENT_COHOMOLOGY_CAM_H_
#define GUDHI_PERSISTENT_COHOMOLOGY_CAM_H_

#include <stdexcept>  // for std::logic_error
#include <utility>    // std::move, std::pair
#include <vector>

#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <boost/unordered_map.hpp>
#endif
#include <boost/intrusive/set.hpp>
#include <boost/intrusive/list.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/Simple_object_pool.h>
#include <gudhi/Persistent_cohomology/Simple_union_find.h>

namespace Gudhi {

namespace persistent_cohomology {

template <typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_column;

struct cam_h_tag;  // for horizontal traversal in the CAM
struct cam_v_tag;  // for vertical traversal in the CAM

using base_hook_cam_h =
    boost::intrusive::list_base_hook<boost::intrusive::tag<cam_h_tag>,
                                     boost::intrusive::link_mode<boost::intrusive::auto_unlink>  // allows .unlink()
                                     >;

using base_hook_cam_v =
    boost::intrusive::list_base_hook<boost::intrusive::tag<cam_v_tag>,
                                     boost::intrusive::link_mode<boost::intrusive::normal_link>  // faster hook,
                                                                                                 // less safe
                                     >;

/** \internal
 * \brief
 *
 */
template <typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_entry : public base_hook_cam_h, public base_hook_cam_v {
 public:
  template <class T1, class T2>
  friend class Persistent_cohomology;
  friend class Persistent_cohomology_column<SimplexKey, ArithmeticElement>;

  using Column = Persistent_cohomology_column<SimplexKey, ArithmeticElement>;

  Persistent_cohomology_entry(SimplexKey key, ArithmeticElement x, Column* self_col)
      : key_(key), coefficient_(x), self_col_(self_col) {}

  ArithmeticElement get_element() const { return coefficient_; }
  Column const* get_column_index() const { return self_col_; }
  SimplexKey get_row_index() const { return key_; }

  SimplexKey key_;
  ArithmeticElement coefficient_;
  Column* self_col_;
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
template <typename SimplexKey, typename ArithmeticElement>
class Persistent_cohomology_column
    : public boost::intrusive::set_base_hook<boost::intrusive::link_mode<boost::intrusive::normal_link>> {
  template <class T1, class T2>
  friend class Persistent_cohomology;

 public:
  using Entry = Persistent_cohomology_entry<SimplexKey, ArithmeticElement>;
  using Col_type = boost::intrusive::list<Entry, boost::intrusive::constant_time_size<false>,
                                          boost::intrusive::base_hook<base_hook_cam_v>>;
  using iterator = typename Col_type::iterator;
  using const_iterator = typename Col_type::const_iterator;

  /** \brief Creates an empty column.*/
  explicit Persistent_cohomology_column(SimplexKey key = -1) : col_(), class_key_(key) {}

  /** \brief Returns true iff the column is null.*/
  [[nodiscard]] bool is_empty() const { return col_.empty(); }

  /** \brief Returns the key of the representative simplex of
   * the set of simplices having this column as annotation vector
   * in the compressed annotation matrix.*/
  SimplexKey class_key() const { return class_key_; }

  /** \brief Lexicographic comparison of two columns.*/
  friend bool operator<(const Persistent_cohomology_column& c1, const Persistent_cohomology_column& c2) {
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

  iterator begin() noexcept { return col_.begin(); }

  const_iterator begin() const noexcept { return col_.begin(); }

  iterator end() noexcept { return col_.end(); }

  const_iterator end() const noexcept { return col_.end(); }

  Col_type col_;
  SimplexKey class_key_;
};

template <typename Cell_key, class Field_operators>
class Compressed_annotation_matrix {
 public:
  using Element = typename Field_operators::Element;
  using Characteristic = typename Field_operators::Element;
  using ID_index = Cell_key;
  using Column = Persistent_cohomology_column<ID_index, Characteristic>;
  using Index = Column const*;
  // Entry type
  using Entry = typename Column::Entry;  // contains 2 list_hooks
  // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
  using Row = boost::intrusive::list<Entry, boost::intrusive::constant_time_size<false>,
                                     boost::intrusive::base_hook<base_hook_cam_h>>;

  using Cam = boost::intrusive::set<Column, boost::intrusive::constant_time_size<false>>;

  Compressed_annotation_matrix(unsigned int numberOfColumns)
      : coeff_field_(),
        ds_repr_(numberOfColumns, nullptr),
        transverse_idx_(),
        cam_(),
        cell_sets_(numberOfColumns),
        column_pool_(),
        entry_pool_() {}

  Compressed_annotation_matrix(const Compressed_annotation_matrix& matrixToCopy) = delete;

  Compressed_annotation_matrix(Compressed_annotation_matrix&& other) noexcept = default;

  ~Compressed_annotation_matrix() {
    // Clean the transversal lists
    for (auto& row : transverse_idx_) {
      // Destruct all the entries
      row.second.clear_and_dispose([&](Entry* p) { p->~Entry(); });
    }
  }

  void set_characteristic(Characteristic characteristic) { coeff_field_.init(characteristic); }

  void set_characteristic(int min, int max) { coeff_field_.init(min, max); }

  const Field_operators& get_field_operators() const { return coeff_field_; }

  void insert_column(ID_index idx, Element e) {
    // Create a column containing only one entry,
    Column* new_col = column_pool_.construct(idx);
    Entry* new_entry = entry_pool_.construct(idx, e, new_col);
    new_col->col_.push_back(*new_entry);
    // and insert it in the matrix, in constant time thanks to the hint cam_.end().
    // Indeed *new_col has the biggest lexicographic value because key is the
    // biggest key used so far.
    cam_.insert(cam_.end(), *new_col);
    // Update the disjoint sets data structure.
    transverse_idx_[idx] = Row();  // insert the new row
    transverse_idx_[idx].push_back(*new_entry);
    cell_sets_.make_set(idx);
    ds_repr_[idx] = new_col;
  }

  const Column& get_column(ID_index columnIndex) {
    auto* col = ds_repr_[cell_sets_.find_set(columnIndex)];
    if (col == nullptr) return empty_column_;
    return *col;
  }

  const Row& get_row(ID_index rowIndex) const { return transverse_idx_.find(rowIndex)->second; }

  void erase_empty_row(ID_index rowIndex) {
    // row should be empty
    GUDHI_CHECK(transverse_idx_.find(rowIndex)->second.size() == 0,
                std::logic_error("Internal problem: row should be empty."));
    transverse_idx_.erase(rowIndex);
  }

  template <class Entry_range>
  void multiply_source_and_add_to(Element coefficient, const Entry_range& sourceColumn, Index targetColumnIndex) {
    Column* targetColumn = const_cast<Column*>(targetColumnIndex);

    // Disconnect the column from the rows in the CAM.
    for (auto& col_entry : targetColumn->col_) {
      col_entry.base_hook_cam_h::unlink();
    }

    // Remove the column from the CAM before modifying its value
    cam_.erase(cam_.iterator_to(*targetColumn));
    // Proceed to the reduction of the column
    plus_equal_column(*targetColumn, sourceColumn, coefficient);

    if (targetColumn->col_.empty()) {  // If the column is null
      ds_repr_[targetColumn->class_key_] = nullptr;
      column_pool_.destroy(targetColumn);  // delete targetColumn;
    } else {
      std::pair<typename Cam::iterator, bool> result_insert_cam;

      // Find whether the column obtained is already in the CAM
      result_insert_cam = cam_.insert(*targetColumn);
      if (result_insert_cam.second) {  // If it was not in the CAM before: insertion has succeeded
        for (auto& col_entry : targetColumn->col_) {
          // re-establish the row links
          transverse_idx_[col_entry.key_].push_front(col_entry);
        }
      } else {  // There is already an identical column in the CAM:
        // merge two disjoint sets.
        cell_sets_.merge_sets(targetColumn->class_key_, result_insert_cam.first->class_key_);

        Cell_key key_tmp = cell_sets_.find_set(targetColumn->class_key_);
        ds_repr_[key_tmp] = &(*(result_insert_cam.first));
        result_insert_cam.first->class_key_ = key_tmp;
        // intrusive containers don't own their elements, we have to release them manually
        targetColumn->col_.clear_and_dispose([&](Entry* p) { entry_pool_.destroy(p); });
        column_pool_.destroy(targetColumn);  // delete targetColumn;
      }
    }
  }

  Compressed_annotation_matrix& operator=(Compressed_annotation_matrix other) = delete;

  Compressed_annotation_matrix& operator=(Compressed_annotation_matrix&& other) noexcept {
    if (&ds_repr_ == &(other.ds_repr_)) return *this;

    for (auto& row : transverse_idx_) {
      // Destruct all potential entries
      row.second.clear_and_dispose([&](Entry* p) { p->~Entry(); });
    }

    coeff_field_ = std::move(other.coeff_field_);
    ds_repr_ = std::move(other.ds_repr_);
    transverse_idx_ = std::move(other.transverse_idx_);
    cam_ = std::move(other.cam_);
    cell_sets_ = std::move(other.cell_sets_);
    column_pool_ = std::move(other.column_pool_);
    entry_pool_ = std::move(other.entry_pool_);
  }

  void print();  // for debug

 private:
  Field_operators coeff_field_;
  std::vector<Column*> ds_repr_;
#if BOOST_VERSION >= 108100
  boost::unordered_flat_map<Cell_key, Row> transverse_idx_;
#else
  boost::unordered_map<Cell_key, Row> transverse_idx_;
#endif
  Cam cam_;
  Simple_union_find<Cell_key> cell_sets_;

  Simple_object_pool<Column> column_pool_;
  Simple_object_pool<Entry> entry_pool_;

  inline static const Column empty_column_; /**< Representative for empty columns. */

  /*
   * Assign:    target <- target + w * other.
   */
  template <class Entry_range>
  void plus_equal_column(Column& target,
                         Entry_range const& other,  // value_type is pair<Cell_key,Arith_element>
                         Characteristic w) {
    auto target_it = target.col_.begin();
    auto other_it = other.begin();
    while (target_it != target.col_.end() && other_it != other.end()) {
      if (target_it->key_ < other_it->first) {
        ++target_it;
      } else {
        if (target_it->key_ > other_it->first) {
          Entry* entry_tmp = entry_pool_.construct(other_it->first, coeff_field_.additive_identity(), &target);

          entry_tmp->coefficient_ = coeff_field_.plus_times_equal(entry_tmp->coefficient_, other_it->second, w);

          target.col_.insert(target_it, *entry_tmp);

          ++other_it;
        } else {  // it1->key == it2->key
          // target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
          target_it->coefficient_ = coeff_field_.plus_times_equal(target_it->coefficient_, other_it->second, w);
          if (target_it->coefficient_ == coeff_field_.additive_identity()) {
            auto tmp_it = target_it;
            ++target_it;
            ++other_it;  // iterators remain valid
            Entry* tmp_entry_ptr = &(*tmp_it);
            target.col_.erase(tmp_it);  // removed from column

            entry_pool_.destroy(tmp_entry_ptr);  // delete from memory
          } else {
            ++target_it;
            ++other_it;
          }
        }
      }
    }
    while (other_it != other.end()) {
      Entry* entry_tmp = entry_pool_.construct(other_it->first, coeff_field_.additive_identity(), &target);
      entry_tmp->coefficient_ = coeff_field_.plus_times_equal(entry_tmp->coefficient_, other_it->second, w);
      target.col_.insert(target.col_.end(), *entry_tmp);

      ++other_it;
    }
  }
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // GUDHI_PERSISTENT_COHOMOLOGY_CAM_H_
