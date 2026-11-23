/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria, Hannah Schreiber
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_PERSISTENT_COHOMOLOGY_BASE_H_
#define GUDHI_PERSISTENT_COHOMOLOGY_BASE_H_

#include <stdexcept>  // std::out_of_range
#include <cstddef>
#include <vector>
#include <initializer_list>
#include <map>
#include <limits>  // std::numeric_limits

#if BOOST_VERSION >= 108100
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <boost/unordered_map.hpp>
#endif

#include <gudhi/Debug_utils.h>
#include <gudhi/Persistent_cohomology/Simple_union_find.h>
#include <gudhi/Persistent_cohomology/Compressed_annotation_matrix.h>

namespace Gudhi {

namespace persistent_cohomology {

template <class CoefficientField, class KeyHandleMap>
class Persistent_cohomology_base {
 public:
  using Cell_key = typename KeyHandleMap::Cell_key;       /**< Data stored for each simplex. */
  using Cell_handle = typename KeyHandleMap::Cell_handle; /**< Handle to specify a simplex. */
  using Arith_element = typename CoefficientField::Element;

 private:
  using A_ds_type = std::vector<std::pair<Cell_key, Arith_element>>;

 public:
  Persistent_cohomology_base(KeyHandleMap* key_handle_map, std::size_t num_cells)
      : key_handle_map_(key_handle_map), coeff_field_(), cam_(num_cells), transverse_idx_(num_cells) {
    if (num_cells > std::numeric_limits<Cell_key>::max()) {
      // num_simplices must be strictly lower than the limit, because a value is reserved for null_key.
      throw std::out_of_range("The number of simplices is more than key type numeric limit.");
    }
  };

  // TODO: try to not duplicate coeff_field_? As it can contain tables.
  // Even though probably negligible compared to the matrix?
  void init_coefficients(int characteristic) {
    coeff_field_.init(characteristic);
    cam_.set_characteristic(characteristic);
  }

  void init_coefficients(int min, int max) {
    coeff_field_.init(min, max);
    cam_.set_characteristic(min, max);
  }

  template <class Cell_boundary_range = std::initializer_list<Cell_handle>, class F>
  void update_cohomology_groups(Cell_handle sigma, const Cell_boundary_range& boundary, bool kill_only, F&& add_pair) {
    // Compute the annotation of the boundary of sigma:
    std::map<Cell_key, Arith_element> map_a_ds;
    _annotation_of_the_boundary(map_a_ds, boundary);
    // Update the cohomology groups:
    if (map_a_ds.empty()) {  // sigma is a creator in all fields represented in coeff_field_
      if (!kill_only) {
        _create_cocycle(sigma, coeff_field_.multiplicative_identity(), coeff_field_.characteristic());
      }
    } else {  // sigma is a destructor in at least a field in coeff_field_
      // Convert map_a_ds to a vector
      A_ds_type a_ds;  // admits reverse iterators
      for (const auto& map_a_ds_ref : map_a_ds) {
        a_ds.push_back(std::pair<Cell_key, Arith_element>(map_a_ds_ref.first, map_a_ds_ref.second));
      }

      Arith_element inv_x, characteristic;
      Arith_element prod = coeff_field_.characteristic();  // Product of characteristic of the fields
      for (auto a_ds_rit = a_ds.rbegin(); (a_ds_rit != a_ds.rend()) && (prod != coeff_field_.multiplicative_identity());
           ++a_ds_rit) {
        std::tie(inv_x, characteristic) = coeff_field_.inverse(a_ds_rit->second, prod);

        if (inv_x != coeff_field_.additive_identity()) {
          _destroy_cocycle(sigma, a_ds, a_ds_rit->first, inv_x, characteristic, std::forward<F>(add_pair));
          prod /= characteristic;
        }
      }
      if (prod != coeff_field_.multiplicative_identity() && !kill_only) {
        _create_cocycle(sigma, coeff_field_.multiplicative_identity(prod), prod);
      }
    }
  }

  template <class F>
  void finalize_computation(F&& add_pair) {
    for (const auto& cocycle : transverse_idx_) {
      std::forward<F>(add_pair)(key_handle_map_->cell(cocycle.first), key_handle_map_->null_cell(), cocycle.second);
    }
  }

 protected:
  KeyHandleMap* key_handle_map_;
  CoefficientField coeff_field_;

  void _create_cocycle(Cell_handle sigma, Arith_element x, Arith_element characteristic) {
    Cell_key key = key_handle_map_->key(sigma);
    cam_.insert_column(key, x);
    transverse_idx_[key] = characteristic;  // insert the new row
  }

 private:
  using Cam = Compressed_annotation_matrix<Cell_key, CoefficientField>;
  using Column = typename Cam::Column;

  Cam cam_;
#if BOOST_VERSION >= 108100
  boost::unordered_flat_map<Cell_key, Arith_element> transverse_idx_;
#else
  boost::unordered_map<Cell_key, Arith_element> transverse_idx_;
#endif

  /*
   * Compute the annotation of the boundary of a simplex.
   */
  template <class Cell_boundary_range>
  void _annotation_of_the_boundary(std::map<Cell_key, Arith_element>& map_a_ds, const Cell_boundary_range& boundary) {
    // traverses the boundary of sigma, keeps track of the annotation vectors,
    // with multiplicity. We used to sum the coefficients directly in
    // annotations_in_boundary by using a map, we now do it later.
    using annotation_t = std::pair<Column const*, int>;
    thread_local std::vector<annotation_t> annotations_in_boundary;
    annotations_in_boundary.clear();
    int sign = -1;  // \in {-1,1} provides the sign in the alternate sum in the boundary.
    Cell_key key;

    for (const auto& sh : boundary) {
      key = key_handle_map_->key(sh);
      if (key != key_handle_map_->null_key()) {  // A simplex with null_key is a killer, and have null annotation
        // Find its annotation vector
        const auto& col = cam_.get_column(key);
        if (!col.is_empty()) {
          // and insert it in annotations_in_boundary with multiplicative factor "sign".
          annotations_in_boundary.emplace_back(&col, sign);
        }
      }
      sign = -sign;
    }
    // Place identical annotations consecutively so we can easily sum their multiplicities.
    std::sort(annotations_in_boundary.begin(), annotations_in_boundary.end(),
              [](annotation_t const& a, annotation_t const& b) { return a.first < b.first; });

    // Sum the annotations with multiplicity, using a map<key,coeff>
    // to represent a sparse vector.
    std::pair<typename std::map<Cell_key, Arith_element>::iterator, bool> result_insert_a_ds;

    for (auto ann_it = annotations_in_boundary.begin(); ann_it != annotations_in_boundary.end(); /**/) {
      Column const* col = ann_it->first;
      int mult = ann_it->second;
      ++ann_it;
      while (ann_it != annotations_in_boundary.end() && ann_it->first == col) {
        mult += ann_it->second;
        ++ann_it;
      }
      // The following test is just a heuristic, it is not required, and it is fine that is misses p == 0.
      if (mult != coeff_field_.additive_identity()) {  // For all columns in the boundary,
        for (auto&& entry_ref : *col) {                // insert every cell in map_a_ds with multiplicity
          Arith_element w_y = coeff_field_.times(entry_ref.get_element(), mult);  // coefficient * multiplicity

          if (w_y != coeff_field_.additive_identity()) {  // if != 0
            result_insert_a_ds = map_a_ds.emplace(entry_ref.get_row_index(), w_y);
            if (!(result_insert_a_ds.second)) {  // if cell_ref.key_ already a Key in map_a_ds
              result_insert_a_ds.first->second = coeff_field_.plus_equal(result_insert_a_ds.first->second, w_y);
              if (result_insert_a_ds.first->second == coeff_field_.additive_identity()) {
                map_a_ds.erase(result_insert_a_ds.first);
              }
            }
          }
        }
      }
    }
  }

  /*  \brief Destroy a cocycle class.
   *
   * The cocycle class is destroyed by the insertion of sigma.
   * The methods proceeds to a reduction of the matrix representing
   * the cohomology groups using Gauss pivoting. The reduction zeros-out
   * the row containing the entry with highest key in
   * a_ds, the annotation of the boundary of simplex sigma. This key
   * is "death_key".*/
  template <class F>
  void _destroy_cocycle(Cell_handle sigma, A_ds_type const& a_ds, Cell_key death_key, Arith_element inv_x,
                        Arith_element characteristic, F&& add_pair) {
    // Create a finite persistent interval for which the interval exists
    std::forward<F>(add_pair)(key_handle_map_->cell(death_key), sigma, characteristic);

    const auto& death_key_row = cam_.get_row(death_key);  // Find the beginning of the row.
    auto death_cocycle_coeff = transverse_idx_.find(death_key);

    auto row_entry_it = death_key_row.begin();

    while (row_entry_it != death_key_row.end()) {  // Traverse all entries in the row at index death_key.
      const auto& entry = *row_entry_it;
      Arith_element w = coeff_field_.times_minus(inv_x, entry.get_element());

      if (w != coeff_field_.additive_identity()) {
        ++row_entry_it;  // has to be done before column addition to avoid invalidating the iterator
        cam_.multiply_source_and_add_to(w, a_ds, entry.get_column_index());
      } else {
        ++row_entry_it;
      }  // If w == 0, pass.
    }

    // Because it is a killer simplex, set the data of sigma to null_key().
    if (characteristic == coeff_field_.characteristic()) {
      key_handle_map_->assign_key(sigma, key_handle_map_->null_key());
    }
    if (death_cocycle_coeff->second == characteristic) {
      cam_.erase_empty_row(death_key);
      transverse_idx_.erase(death_cocycle_coeff);
    } else {
      death_cocycle_coeff->second /= characteristic;
    }
  }
};

template <class CoefficientField, class KeyHandleMap>
class Polytopal_persistent_cohomology_base : public Persistent_cohomology_base<CoefficientField, KeyHandleMap> {
 public:
  using Base = Persistent_cohomology_base<CoefficientField, KeyHandleMap>;
  using Cell_key = typename Base::Cell_key;       /**< Data stored for each simplex. */
  using Cell_handle = typename Base::Cell_handle; /**< Handle to specify a simplex. */

  Polytopal_persistent_cohomology_base(KeyHandleMap* key_handle_map, std::size_t num_cells, std::size_t num_vertices)
      : Base(key_handle_map, num_cells),
        current_vertex_index_(0),
        vertex_sets_(num_vertices),
        vertices_(num_vertices),
        zero_cocycles_() {};

  template <class Cell_boundary_range = std::initializer_list<Cell_handle>, class F>
  void update_cohomology_groups(Cell_handle sigma, const Cell_boundary_range& boundary, bool kill_only, F&& add_pair) {
    std::vector<Cell_handle> b(boundary.begin(), boundary.end());
    GUDHI_CHECK(b.size() != 1, std::invalid_argument("Boundary of polytopial cell cannot be of size 1."));
    std::size_t dim = b.empty() ? 0 : b.size() - 1;  // true for at least dim 0 and 1
    switch (dim) {
      case 0:
        _update_vertex_cohomology_groups(sigma);
        break;
      case 1:
        _update_edge_cohomology_groups(sigma, Base::key_handle_map_->key(b[0]), Base::key_handle_map_->key(b[1]),
                                       kill_only, std::forward<F>(add_pair));
        break;
      default:
        Base::update_cohomology_groups(sigma, b, kill_only, std::forward<F>(add_pair));
        break;
    }
  }

  template <class F>
  void finalize_computation(F&& add_pair) {
    for (Cell_key i = 0; i < static_cast<Cell_key>(vertices_.size()); ++i) {
      if (vertex_sets_.get_parent(i) == i  // root of its tree
          && zero_cocycles_.find(i) == zero_cocycles_.end()) {
        add_pair(vertices_[i], Base::key_handle_map_->null_cell(), Base::coeff_field_.characteristic());
      }
    }
    for (const auto& zero_idx : zero_cocycles_) {
      add_pair(vertices_[zero_idx.second], Base::key_handle_map_->null_cell(), Base::coeff_field_.characteristic());
    }
    Base::finalize_computation(std::forward<F>(add_pair));
  }

 private:
  Cell_key current_vertex_index_;
  /*  Disjoint sets data structure to link the model of KeyHandleMap
   * with the compressed annotation matrix.
   * ds_rank_ is a property map Cell_key -> int, ds_parent_ is a property map
   * Cell_key -> simplex_key_t */
  Simple_union_find<Cell_key> vertex_sets_;
  // maps vertex position wrt other vertices to position in the whole filtration.
  std::vector<Cell_handle> vertices_;
#if BOOST_VERSION >= 108100
  /*  Dictionary establishing the correspondence between the Cell_key of
   * the root vertex in the union-find ds and the Cell_key of the vertex which
   * created the connected component as a 0-dimension homology feature.*/
  boost::unordered_flat_map<Cell_key, Cell_key> zero_cocycles_;
#else
  boost::unordered_map<Simplex_key, Simplex_key> zero_cocycles_;
#endif

  void _update_vertex_cohomology_groups(Cell_handle vertex) {
    Base::key_handle_map_->assign_key(vertex, current_vertex_index_);
    vertex_sets_.make_set(current_vertex_index_);
    vertices_[current_vertex_index_] = vertex;
    ++current_vertex_index_;
  }

  template <class F>
  void _update_edge_cohomology_groups(Cell_handle sigma, Cell_key u, Cell_key v, bool kill_only, F&& add_pair) {
    Cell_key ku = vertex_sets_.find_set(u);
    Cell_key kv = vertex_sets_.find_set(v);

    if (ku != kv) {  // Destroy a connected component
      vertex_sets_.merge_sets(ku, kv);
      // Keys of the simplices which created the connected components containing
      // respectively u and v.
      Cell_key idx_co_u, idx_co_v;
      auto map_it_u = zero_cocycles_.find(ku);
      // If the index of the cocycle representing the class is already ku.
      if (map_it_u == zero_cocycles_.end()) {
        idx_co_u = ku;
      } else {
        idx_co_u = map_it_u->second;
      }

      auto map_it_v = zero_cocycles_.find(kv);
      // If the index of the cocycle representing the class is already kv.
      if (map_it_v == zero_cocycles_.end()) {
        idx_co_v = kv;
      } else {
        idx_co_v = map_it_v->second;
      }

      if (idx_co_u < idx_co_v) {  // Kill cocycle [idx_co_v], which is younger.
        std::forward<F>(add_pair)(vertices_[idx_co_v], sigma, Base::coeff_field_.characteristic());
        // Maintain the index of the 0-cocycle alive.
        if (kv != idx_co_v) {
          zero_cocycles_.erase(map_it_v);
        }
        if (kv == vertex_sets_.find_set(kv)) {
          if (ku != idx_co_u) {
            zero_cocycles_.erase(map_it_u);
          }
          zero_cocycles_[kv] = idx_co_u;
        }
      } else {  // Kill cocycle [idx_co_u], which is younger.
        std::forward<F>(add_pair)(vertices_[idx_co_u], sigma, Base::coeff_field_.characteristic());
        // Maintain the index of the 0-cocycle alive.
        if (ku != idx_co_u) {
          zero_cocycles_.erase(map_it_u);
        }
        if (ku == vertex_sets_.find_set(ku)) {
          if (kv != idx_co_v) {
            zero_cocycles_.erase(map_it_v);
          }
          zero_cocycles_[ku] = idx_co_v;
        }
      }
      Base::key_handle_map_->assign_key(sigma, Base::key_handle_map_->null_key());
    } else {
      // If ku == kv, same connected component: create a 1-cocycle class.
      if (!kill_only) {
        Base::_create_cocycle(sigma, Base::coeff_field_.multiplicative_identity(), Base::coeff_field_.characteristic());
      }
    }
  }
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // GUDHI_PERSISTENT_COHOMOLOGY_BASE_H_
