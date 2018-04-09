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

#ifndef PERSISTENT_COHOMOLOGY_H_
#define PERSISTENT_COHOMOLOGY_H_

#include <gudhi/Persistent_cohomology/Persistent_cohomology_column.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>
#include <gudhi/Simple_object_pool.h>

#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/intrusive/list.hpp>

#include <map>
#include <utility>
#include <list>
#include <vector>
#include <set>
#include <fstream>  // std::ofstream
#include <limits>  // for numeric_limits<>
#include <tuple>
#include <algorithm>
#include <string>
#include <stdexcept>  // for std::out_of_range

namespace Gudhi {

namespace persistent_cohomology {

/** \brief Computes the persistent cohomology of a filtered complex.
 *
 * \ingroup persistent_cohomology
 * 
 * The computation is implemented with a Compressed Annotation Matrix
 * (CAM)\cite DBLP:conf/esa/BoissonnatDM13,
 * and is adapted to the computation of Multi-Field Persistent Homology (MF)
 * \cite boissonnat:hal-00922572 .
 *
 * \implements PersistentHomology
 *
 */
// TODO(CM): Memory allocation policy: classic, use a mempool, etc.
template<class FilteredComplex, class CoefficientField>
class Persistent_cohomology {
 public:
  typedef FilteredComplex Complex_ds;
  // Data attached to each simplex to interface with a Property Map.
  typedef typename Complex_ds::Simplex_key Simplex_key;
  typedef typename Complex_ds::Simplex_handle Simplex_handle;
  typedef typename Complex_ds::Filtration_value Filtration_value;
  typedef typename CoefficientField::Element Arith_element;
  // Compressed Annotation Matrix types:
  // Column type
  typedef Persistent_cohomology_column<Simplex_key, Arith_element> Column;  // contains 1 set_hook
  // Cell type
  typedef typename Column::Cell Cell;   // contains 2 list_hooks
  // Remark: constant_time_size must be false because base_hook_cam_h has auto_unlink link_mode
  typedef boost::intrusive::list<Cell,
      boost::intrusive::constant_time_size<false>,
      boost::intrusive::base_hook<base_hook_cam_h> > Hcell;

  typedef boost::intrusive::set<Column,
      boost::intrusive::constant_time_size<false> > Cam;
  // Sparse column type for the annotation of the boundary of an element.
  typedef std::vector<std::pair<Simplex_key, Arith_element> > A_ds_type;
  // Persistent interval type. The Arith_element field is used for the multi-field framework.
  typedef std::tuple<Simplex_handle, Simplex_handle, Arith_element> Persistent_interval;

  /** \brief Initializes the Persistent_cohomology class.
   *
   * @param[in] cpx Complex for which the persistent homology is computed.
   * cpx is a model of FilteredComplex
   * @exception std::out_of_range In case the number of simplices is more than Simplex_key type numeric limit.
   */
  explicit Persistent_cohomology(Complex_ds& cpx)
      : cpx_(&cpx),
        dim_max_(cpx.dimension()),                       // upper bound on the dimension of the simplices
        coeff_field_(),                                  // initialize the field coefficient structure.
        num_simplices_(cpx_->num_simplices()),           // num_simplices save to avoid to call thrice the function
        ds_rank_(num_simplices_),                        // union-find
        ds_parent_(num_simplices_),                      // union-find
        ds_repr_(num_simplices_, NULL),                  // union-find -> annotation vectors
        dsets_(&ds_rank_[0], &ds_parent_[0]),            // union-find
        cam_(),                                          // collection of annotation vectors
        zero_cocycles_(),                                // union-find -> Simplex_key of creator for 0-homology
        transverse_idx_(),                               // key -> row
        persistent_pairs_(),
        interval_length_policy(&cpx, 0),
        column_pool_(),  // memory pools for the CAM
        cell_pool_() {
    if (cpx_->num_simplices() > std::numeric_limits<Simplex_key>::max()) {
      // num_simplices must be strictly lower than the limit, because a value is reserved for null_key.
      throw std::out_of_range("The number of simplices is more than Simplex_key type numeric limit.");
    }
    Simplex_key idx_fil = 0;
    for (auto sh : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(sh, idx_fil);
      ++idx_fil;
      dsets_.make_set(cpx_->key(sh));
    }
  }

  /** \brief Initializes the Persistent_cohomology class.
   *
   * @param[in] cpx Complex for which the persistent homology is compiuted.
   * cpx is a model of FilteredComplex
   *
   * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the
   *                                complex is computed. If false, it is ignored. Default is false.
   */
  Persistent_cohomology(Complex_ds& cpx, bool persistence_dim_max)
      : Persistent_cohomology(cpx) {
    if (persistence_dim_max) {
      ++dim_max_;
    }
  }

  ~Persistent_cohomology() {
    // Clean the transversal lists
    for (auto & transverse_ref : transverse_idx_) {
      // Destruct all the cells
      transverse_ref.second.row_->clear_and_dispose([&](Cell*p){p->~Cell();});
      delete transverse_ref.second.row_;
    }
  }

 private:
  struct length_interval {
    length_interval(Complex_ds * cpx, Filtration_value min_length)
        : cpx_(cpx),
          min_length_(min_length) {
    }

    bool operator()(Simplex_handle sh1, Simplex_handle sh2) {
      return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_;
    }

    void set_length(Filtration_value new_length) {
      min_length_ = new_length;
    }

    Complex_ds * cpx_;
    Filtration_value min_length_;
  };

 public:
  /** \brief Initializes the coefficient field.*/
  void init_coefficients(int charac) {
    coeff_field_.init(charac);
  }
  /** \brief Initializes the coefficient field for multi-field persistent homology.*/
  void init_coefficients(int charac_min, int charac_max) {
    coeff_field_.init(charac_min, charac_max);
  }

  /** SK TEMP increment the dimension by 1
   */
  void increment_dimension() {
    dim_max_++;
  }

  /** \brief Compute the persistent homology of the filtered simplicial
   * complex.
   *
   * @param[in] min_interval_length the computation discards all intervals of length
   *                                less or equal than min_interval_length
   *
   * Assumes that the filtration provided by the simplicial complex is
   * valid. Undefined behavior otherwise. */
  void compute_persistent_cohomology(Filtration_value min_interval_length = 0) {
    interval_length_policy.set_length(min_interval_length);
    // Compute all finite intervals
    for (auto sh : cpx_->filtration_simplex_range()) {
      int dim_simplex = cpx_->dimension(sh);
      switch (dim_simplex) {
        case 0:
          break;
        case 1:
          update_cohomology_groups_edge(sh);
          break;
        default:
          update_cohomology_groups(sh, dim_simplex);
          break;
      }
    }
    // Compute infinite intervals of dimension 0
    Simplex_key key;
    for (auto v_sh : cpx_->skeleton_simplex_range(0)) {  // for all 0-dimensional simplices
      key = cpx_->key(v_sh);

      if (ds_parent_[key] == key  // root of its tree
      && zero_cocycles_.find(key) == zero_cocycles_.end()) {
        persistent_pairs_.emplace_back(
            cpx_->simplex(key), cpx_->null_simplex(), coeff_field_.characteristic());
      }
    }
    for (auto zero_idx : zero_cocycles_) {
      persistent_pairs_.emplace_back(
          cpx_->simplex(zero_idx.second), cpx_->null_simplex(), coeff_field_.characteristic());
    }
    // Compute infinite interval of dimension > 0
    for (auto cocycle : transverse_idx_) {
      persistent_pairs_.emplace_back(
          cpx_->simplex(cocycle.first), cpx_->null_simplex(), cocycle.second.characteristics_);
    }
  }

 private:
  /** \brief Update the cohomology groups under the insertion of an edge.
   *
   * The 0-homology is maintained with a simple Union-Find data structure, which
   * explains the existance of a specific function of edge insertions. */
  void update_cohomology_groups_edge(Simplex_handle sigma) {
    Simplex_handle u, v;
    boost::tie(u, v) = cpx_->endpoints(sigma);

    Simplex_key ku = dsets_.find_set(cpx_->key(u));
    Simplex_key kv = dsets_.find_set(cpx_->key(v));

    if (ku != kv) {        // Destroy a connected component
      dsets_.link(ku, kv);
      // Keys of the simplices which created the connected components containing
      // respectively u and v.
      Simplex_key idx_coc_u, idx_coc_v;
      auto map_it_u = zero_cocycles_.find(ku);
      // If the index of the cocycle representing the class is already ku.
      if (map_it_u == zero_cocycles_.end()) {
        idx_coc_u = ku;
      } else {
        idx_coc_u = map_it_u->second;
      }

      auto map_it_v = zero_cocycles_.find(kv);
      // If the index of the cocycle representing the class is already kv.
      if (map_it_v == zero_cocycles_.end()) {
        idx_coc_v = kv;
      } else {
        idx_coc_v = map_it_v->second;
      }

      if (cpx_->filtration(cpx_->simplex(idx_coc_u))
          < cpx_->filtration(cpx_->simplex(idx_coc_v))) {  // Kill cocycle [idx_coc_v], which is younger.
        if (interval_length_policy(cpx_->simplex(idx_coc_v), sigma)) {
          persistent_pairs_.emplace_back(
              cpx_->simplex(idx_coc_v), sigma, coeff_field_.characteristic());
        }
        // Maintain the index of the 0-cocycle alive.
        if (kv != idx_coc_v) {
          zero_cocycles_.erase(map_it_v);
        }
        if (kv == dsets_.find_set(kv)) {
          if (ku != idx_coc_u) {
            zero_cocycles_.erase(map_it_u);
          }
          zero_cocycles_[kv] = idx_coc_u;
        }
      } else {  // Kill cocycle [idx_coc_u], which is younger.
        if (interval_length_policy(cpx_->simplex(idx_coc_u), sigma)) {
          persistent_pairs_.emplace_back(
              cpx_->simplex(idx_coc_u), sigma, coeff_field_.characteristic());
        }
        // Maintain the index of the 0-cocycle alive.
        if (ku != idx_coc_u) {
          zero_cocycles_.erase(map_it_u);
        }
        if (ku == dsets_.find_set(ku)) {
          if (kv != idx_coc_v) {
            zero_cocycles_.erase(map_it_v);
          }
          zero_cocycles_[ku] = idx_coc_v;
        }
      }
      cpx_->assign_key(sigma, cpx_->null_key());
    } else {  // If ku == kv, same connected component: create a 1-cocycle class.
      create_cocycle(sigma, coeff_field_.multiplicative_identity(), coeff_field_.characteristic());
    }
  }

  /*
   * Compute the annotation of the boundary of a simplex.
   */
  void annotation_of_the_boundary(
      std::map<Simplex_key, Arith_element> & map_a_ds, Simplex_handle sigma,
      int dim_sigma) {
    // traverses the boundary of sigma, keeps track of the annotation vectors,
    // with multiplicity. We used to sum the coefficients directly in
    // annotations_in_boundary by using a map, we now do it later.
    typedef std::pair<Column *, int> annotation_t;
    thread_local std::vector<annotation_t> annotations_in_boundary;
    annotations_in_boundary.clear();
    int sign = 1 - 2 * (dim_sigma % 2);  // \in {-1,1} provides the sign in the
                                         // alternate sum in the boundary.
    Simplex_key key;
    Column * curr_col;

    for (auto sh : cpx_->boundary_simplex_range(sigma)) {
      key = cpx_->key(sh);
      if (key != cpx_->null_key()) {  // A simplex with null_key is a killer, and have null annotation
        // Find its annotation vector
        curr_col = ds_repr_[dsets_.find_set(key)];
        if (curr_col != NULL) {  // and insert it in annotations_in_boundary with multyiplicative factor "sign".
          annotations_in_boundary.emplace_back(curr_col, sign);
        }
      }
      sign = -sign;
    }
    // Place identical annotations consecutively so we can easily sum their multiplicities.
    std::sort(annotations_in_boundary.begin(), annotations_in_boundary.end(),
              [](annotation_t const& a, annotation_t const& b) { return a.first < b.first; });

    // Sum the annotations with multiplicity, using a map<key,coeff>
    // to represent a sparse vector.
    std::pair<typename std::map<Simplex_key, Arith_element>::iterator, bool> result_insert_a_ds;

    for (auto ann_it = annotations_in_boundary.begin(); ann_it != annotations_in_boundary.end(); /**/) {
      Column* col = ann_it->first;
      int mult = ann_it->second;
      while (++ann_it != annotations_in_boundary.end() && ann_it->first == col) {
        mult += ann_it->second;
      }
      // The following test is just a heuristic, it is not required, and it is fine that is misses p == 0.
      if (mult != coeff_field_.additive_identity()) {  // For all columns in the boundary,
        for (auto cell_ref : col->col_) {  // insert every cell in map_a_ds with multiplicity
          Arith_element w_y = coeff_field_.times(cell_ref.coefficient_, mult);  // coefficient * multiplicity

          if (w_y != coeff_field_.additive_identity()) {  // if != 0
            result_insert_a_ds = map_a_ds.insert(std::pair<Simplex_key, Arith_element>(cell_ref.key_, w_y));
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

  /*
   * Update the cohomology groups under the insertion of a simplex.
   */
  void update_cohomology_groups(Simplex_handle sigma, int dim_sigma) {
// Compute the annotation of the boundary of sigma:
    std::map<Simplex_key, Arith_element> map_a_ds;
    annotation_of_the_boundary(map_a_ds, sigma, dim_sigma);
// Update the cohomology groups:
    if (map_a_ds.empty()) {  // sigma is a creator in all fields represented in coeff_field_
      if (dim_sigma < dim_max_) {
        create_cocycle(sigma, coeff_field_.multiplicative_identity(),
                       coeff_field_.characteristic());
      }
    } else {        // sigma is a destructor in at least a field in coeff_field_
      // Convert map_a_ds to a vector
      A_ds_type a_ds;  // admits reverse iterators
      for (auto map_a_ds_ref : map_a_ds) {
        a_ds.push_back(
            std::pair<Simplex_key, Arith_element>(map_a_ds_ref.first,
                                                  map_a_ds_ref.second));
      }

      Arith_element inv_x, charac;
      Arith_element prod = coeff_field_.characteristic();  // Product of characteristic of the fields
      for (auto a_ds_rit = a_ds.rbegin();
          (a_ds_rit != a_ds.rend())
              && (prod != coeff_field_.multiplicative_identity()); ++a_ds_rit) {
        std::tie(inv_x, charac) = coeff_field_.inverse(a_ds_rit->second, prod);

        if (inv_x != coeff_field_.additive_identity()) {
          destroy_cocycle(sigma, a_ds, a_ds_rit->first, inv_x, charac);
          prod /= charac;
        }
      }
      if (prod != coeff_field_.multiplicative_identity()
          && dim_sigma < dim_max_) {
        create_cocycle(sigma, coeff_field_.multiplicative_identity(prod), prod);
      }
    }
  }

  /*  \brief Create a new cocycle class.
   *
   * The class is created by the insertion of the simplex sigma.
   * The methods adds a cocycle, representing the new cocycle class,
   * to the matrix representing the cohomology groups.
   * The new cocycle has value 0 on every simplex except on sigma
   * where it worths 1.*/
  void create_cocycle(Simplex_handle sigma, Arith_element x,
                      Arith_element charac) {
    Simplex_key key = cpx_->key(sigma);
    // Create a column containing only one cell,
    Column * new_col = column_pool_.construct(key);
    Cell * new_cell = cell_pool_.construct(key, x, new_col);
    new_col->col_.push_back(*new_cell);
    // and insert it in the matrix, in constant time thanks to the hint cam_.end().
    // Indeed *new_col has the biggest lexicographic value because key is the
    // biggest key used so far.
    cam_.insert(cam_.end(), *new_col);
    // Update the disjoint sets data structure.
    Hcell * new_hcell = new Hcell;
    new_hcell->push_back(*new_cell);
    transverse_idx_[key] = cocycle(charac, new_hcell);  // insert the new row
    ds_repr_[key] = new_col;
  }

  /*  \brief Destroy a cocycle class.
   *
   * The cocycle class is destroyed by the insertion of sigma.
   * The methods proceeds to a reduction of the matrix representing
   * the cohomology groups using Gauss pivoting. The reduction zeros-out
   * the row containing the cell with highest key in
   * a_ds, the annotation of the boundary of simplex sigma. This key
   * is "death_key".*/
  void destroy_cocycle(Simplex_handle sigma, A_ds_type const& a_ds,
                       Simplex_key death_key, Arith_element inv_x,
                       Arith_element charac) {
    // Create a finite persistent interval for which the interval exists
    if (interval_length_policy(cpx_->simplex(death_key), sigma)) {
      persistent_pairs_.emplace_back(cpx_->simplex(death_key)  // creator
          , sigma                                              // destructor
          , charac);                                           // fields
    }

    auto death_key_row = transverse_idx_.find(death_key);  // Find the beginning of the row.
    std::pair<typename Cam::iterator, bool> result_insert_cam;

    auto row_cell_it = death_key_row->second.row_->begin();

    while (row_cell_it != death_key_row->second.row_->end()) {  // Traverse all cells in
      // the row at index death_key.
      Arith_element w = coeff_field_.times_minus(inv_x, row_cell_it->coefficient_);

      if (w != coeff_field_.additive_identity()) {
        Column * curr_col = row_cell_it->self_col_;
        ++row_cell_it;
        // Disconnect the column from the rows in the CAM.
        for (auto& col_cell : curr_col->col_) {
          col_cell.base_hook_cam_h::unlink();
        }

        // Remove the column from the CAM before modifying its value
        cam_.erase(cam_.iterator_to(*curr_col));
        // Proceed to the reduction of the column
        plus_equal_column(*curr_col, a_ds, w);

        if (curr_col->col_.empty()) {  // If the column is null
          ds_repr_[curr_col->class_key_] = NULL;
          column_pool_.destroy(curr_col);  // delete curr_col;
        } else {
          // Find whether the column obtained is already in the CAM
          result_insert_cam = cam_.insert(*curr_col);
          if (result_insert_cam.second) {  // If it was not in the CAM before: insertion has succeeded
            for (auto& col_cell : curr_col->col_) {
              // re-establish the row links
              transverse_idx_[col_cell.key_].row_->push_front(col_cell);
            }
          } else {  // There is already an identical column in the CAM:
            // merge two disjoint sets.
            dsets_.link(curr_col->class_key_,
                        result_insert_cam.first->class_key_);

            Simplex_key key_tmp = dsets_.find_set(curr_col->class_key_);
            ds_repr_[key_tmp] = &(*(result_insert_cam.first));
            result_insert_cam.first->class_key_ = key_tmp;
            // intrusive containers don't own their elements, we have to release them manually
            curr_col->col_.clear_and_dispose([&](Cell*p){cell_pool_.destroy(p);});
            column_pool_.destroy(curr_col);  // delete curr_col;
          }
        }
      } else {
        ++row_cell_it;
      }  // If w == 0, pass.
    }

    // Because it is a killer simplex, set the data of sigma to null_key().
    if (charac == coeff_field_.characteristic()) {
      cpx_->assign_key(sigma, cpx_->null_key());
    }
    if (death_key_row->second.characteristics_ == charac) {
      delete death_key_row->second.row_;
      transverse_idx_.erase(death_key_row);
    } else {
      death_key_row->second.characteristics_ /= charac;
    }
  }

  /*
   * Assign:    target <- target + w * other.
   */
  void plus_equal_column(Column & target, A_ds_type const& other  // value_type is pair<Simplex_key,Arith_element>
                         , Arith_element w) {
    auto target_it = target.col_.begin();
    auto other_it = other.begin();
    while (target_it != target.col_.end() && other_it != other.end()) {
      if (target_it->key_ < other_it->first) {
        ++target_it;
      } else {
        if (target_it->key_ > other_it->first) {
          Cell * cell_tmp = cell_pool_.construct(Cell(other_it->first   // key
              , coeff_field_.additive_identity(), &target));

          cell_tmp->coefficient_ = coeff_field_.plus_times_equal(cell_tmp->coefficient_, other_it->second, w);

          target.col_.insert(target_it, *cell_tmp);

          ++other_it;
        } else {  // it1->key == it2->key
          // target_it->coefficient_ <- target_it->coefficient_ + other_it->second * w
          target_it->coefficient_ = coeff_field_.plus_times_equal(target_it->coefficient_, other_it->second, w);
          if (target_it->coefficient_ == coeff_field_.additive_identity()) {
            auto tmp_it = target_it;
            ++target_it;
            ++other_it;   // iterators remain valid
            Cell * tmp_cell_ptr = &(*tmp_it);
            target.col_.erase(tmp_it);  // removed from column

            cell_pool_.destroy(tmp_cell_ptr);  // delete from memory
          } else {
            ++target_it;
            ++other_it;
          }
        }
      }
    }
    while (other_it != other.end()) {
      Cell * cell_tmp = cell_pool_.construct(Cell(other_it->first, coeff_field_.additive_identity(), &target));
      cell_tmp->coefficient_ = coeff_field_.plus_times_equal(cell_tmp->coefficient_, other_it->second, w);
      target.col_.insert(target.col_.end(), *cell_tmp);

      ++other_it;
    }
  }

  /*
   * Compare two intervals by length.
   */
  struct cmp_intervals_by_length {
    explicit cmp_intervals_by_length(Complex_ds * sc)
        : sc_(sc) {
    }
    bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
      return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
          > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
    }
    Complex_ds * sc_;
  };

 public:
  /** \brief Output the persistence diagram in ostream.
   *
   * The file format is the following:
   *    p1*...*pr   dim b d
   *
   * where "dim" is the dimension of the homological feature,
   * b and d are respectively the birth and death of the feature and
   * p1*...*pr is the product of prime numbers pi such that the homology
   * feature exists in homology with Z/piZ coefficients.
   */
  void output_diagram(std::ostream& ostream = std::cout) {
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(persistent_pairs_), std::end(persistent_pairs_), cmp);
    bool has_infinity = std::numeric_limits<Filtration_value>::has_infinity;
    for (auto pair : persistent_pairs_) {
      // Special case on windows, inf is "1.#INF" (cf. unitary tests and R package TDA)
      if (has_infinity && cpx_->filtration(get<1>(pair)) == std::numeric_limits<Filtration_value>::infinity()) {
        ostream << get<2>(pair) << "  " << cpx_->dimension(get<0>(pair)) << " "
          << cpx_->filtration(get<0>(pair)) << " inf " << std::endl;
      } else {
        ostream << get<2>(pair) << "  " << cpx_->dimension(get<0>(pair)) << " "
          << cpx_->filtration(get<0>(pair)) << " "
          << cpx_->filtration(get<1>(pair)) << " " << std::endl;
      }
    }
  }

  void write_output_diagram(std::string diagram_name) {
    std::ofstream diagram_out(diagram_name.c_str());
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(persistent_pairs_), std::end(persistent_pairs_), cmp);
    bool has_infinity = std::numeric_limits<Filtration_value>::has_infinity;
    for (auto pair : persistent_pairs_) {
      // Special case on windows, inf is "1.#INF"
      if (has_infinity && cpx_->filtration(get<1>(pair)) == std::numeric_limits<Filtration_value>::infinity()) {
        diagram_out << cpx_->dimension(get<0>(pair)) << " "
              << cpx_->filtration(get<0>(pair)) << " inf" << std::endl;
      } else {
        diagram_out << cpx_->dimension(get<0>(pair)) << " "
              << cpx_->filtration(get<0>(pair)) << " "
              << cpx_->filtration(get<1>(pair)) << std::endl;
      }
    }
  }

  /** @brief Returns Betti numbers.
   * @return A vector of Betti numbers.
   */
  std::vector<int> betti_numbers() const {
    // Init Betti numbers vector with zeros until Simplicial complex dimension
    std::vector<int> betti_numbers(dim_max_, 0);

    for (auto pair : persistent_pairs_) {
      // Count never ended persistence intervals
      if (cpx_->null_simplex() == get<1>(pair)) {
        // Increment corresponding betti number
        betti_numbers[cpx_->dimension(get<0>(pair))] += 1;
      }
    }
    return betti_numbers;
  }

  /** @brief Returns the Betti number of the dimension passed by parameter.
   * @param[in] dimension The Betti number dimension to get.
   * @return Betti number of the given dimension
   *
   */
  int betti_number(int dimension) const {
    int betti_number = 0;

    for (auto pair : persistent_pairs_) {
      // Count never ended persistence intervals
      if (cpx_->null_simplex() == get<1>(pair)) {
        if (cpx_->dimension(get<0>(pair)) == dimension) {
          // Increment betti number found
          ++betti_number;
        }
      }
    }
    return betti_number;
  }

  /** @brief Returns the persistent Betti numbers.
   * @param[in] from The persistence birth limit to be added in the number \f$(persistent birth \leq from)\f$.
   * @param[in] to The persistence death limit to be added in the number  \f$(persistent death > to)\f$.
   * @return A vector of persistent Betti numbers.
   */
  std::vector<int> persistent_betti_numbers(Filtration_value from, Filtration_value to) const {
    // Init Betti numbers vector with zeros until Simplicial complex dimension
    std::vector<int> betti_numbers(dim_max_, 0);
    for (auto pair : persistent_pairs_) {
      // Count persistence intervals that covers the given interval
      // null_simplex test : if the function is called with to=+infinity, we still get something useful. And it will
      // still work if we change the complex filtration function to reject null simplices.
      if (cpx_->filtration(get<0>(pair)) <= from &&
          (get<1>(pair) == cpx_->null_simplex() || cpx_->filtration(get<1>(pair)) > to)) {
        // Increment corresponding betti number
        betti_numbers[cpx_->dimension(get<0>(pair))] += 1;
      }
    }
    return betti_numbers;
  }

  /** @brief Returns the persistent Betti number of the dimension passed by parameter.
   * @param[in] dimension The Betti number dimension to get.
   * @param[in] from The persistence birth limit to be added in the number \f$(persistent birth \leq from)\f$.
   * @param[in] to The persistence death limit to be added in the number  \f$(persistent death > to)\f$.
   * @return Persistent Betti number of the given dimension
   */
  int persistent_betti_number(int dimension, Filtration_value from, Filtration_value to) const {
    int betti_number = 0;

    for (auto pair : persistent_pairs_) {
      // Count persistence intervals that covers the given interval
      // null_simplex test : if the function is called with to=+infinity, we still get something useful. And it will
      // still work if we change the complex filtration function to reject null simplices.
      if (cpx_->filtration(get<0>(pair)) <= from &&
          (get<1>(pair) == cpx_->null_simplex() || cpx_->filtration(get<1>(pair)) > to)) {
        if (cpx_->dimension(get<0>(pair)) == dimension) {
          // Increment betti number found
          ++betti_number;
        }
      }
    }
    return betti_number;
  }

  /** @brief Returns the persistent pairs.
   * @return Persistent pairs
   *
   */
  const std::vector<Persistent_interval>& get_persistent_pairs() const {
    return persistent_pairs_;
  }

  /** @brief Returns persistence intervals for a given dimension.
   * @param[in] dimension Dimension to get the birth and death pairs from.
   * @return A vector of persistence intervals (birth and death) on a fixed dimension.
   */
  std::vector< std::pair< Filtration_value , Filtration_value > >
  intervals_in_dimension(int dimension) {
    std::vector< std::pair< Filtration_value , Filtration_value > > result;
    // auto && pair, to avoid unnecessary copying
    for (auto && pair : persistent_pairs_) {
      if (cpx_->dimension(get<0>(pair)) == dimension) {
        result.emplace_back(cpx_->filtration(get<0>(pair)), cpx_->filtration(get<1>(pair)));
      }
    }
    return result;
  }

 private:
  /*
   * Structure representing a cocycle.
   */
  struct cocycle {
    cocycle()
        : row_(nullptr),
          characteristics_() {
    }
    cocycle(Arith_element characteristics, Hcell * row)
        : row_(row),
          characteristics_(characteristics) {
    }

    Hcell * row_;                    // points to the corresponding row in the CAM
    Arith_element characteristics_;  // product of field characteristics for which the cocycle exist
  };

 public:
  Complex_ds * cpx_;
  int dim_max_;
  CoefficientField coeff_field_;
  size_t num_simplices_;

  /*  Disjoint sets data structure to link the model of FilteredComplex
   * with the compressed annotation matrix.
   * ds_rank_ is a property map Simplex_key -> int, ds_parent_ is a property map
   * Simplex_key -> simplex_key_t */
  std::vector<int> ds_rank_;
  std::vector<Simplex_key> ds_parent_;
  std::vector<Column *> ds_repr_;
  boost::disjoint_sets<int *, Simplex_key *> dsets_;
  /* The compressed annotation matrix fields.*/
  Cam cam_;
  /*  Dictionary establishing the correspondance between the Simplex_key of
   * the root vertex in the union-find ds and the Simplex_key of the vertex which
   * created the connected component as a 0-dimension homology feature.*/
  std::map<Simplex_key, Simplex_key> zero_cocycles_;
  /*  Key -> row. */
  std::map<Simplex_key, cocycle> transverse_idx_;
  /* Persistent intervals. */
  std::vector<Persistent_interval> persistent_pairs_;
  length_interval interval_length_policy;

  Simple_object_pool<Column> column_pool_;
  Simple_object_pool<Cell> cell_pool_;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_H_
