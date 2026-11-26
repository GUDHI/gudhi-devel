/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2026/11 Hannah Schreiber: Re-factorization with base classes for future generalization
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef GUDHI_PERSISTENT_COHOMOLOGY_H_
#define GUDHI_PERSISTENT_COHOMOLOGY_H_

#include <iostream>
#include <utility>
#include <vector>
#include <fstream>  // std::ofstream
#include <tuple>
#include <algorithm>
#include <string>

#include <gudhi/Persistent_cohomology/Field_Zp.h>  //only here for external use....
#include <gudhi/Persistent_cohomology/Persistent_cohomology_base.h>

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
template <class FilteredComplex, class CoefficientField>
class Persistent_cohomology {
 public:
  // Data attached to each simplex to interface with a Property Map.

  /** \brief Data stored for each simplex. */
  using Simplex_key = typename FilteredComplex::Simplex_key;
  /** \brief Handle to specify a simplex. */
  using Simplex_handle = typename FilteredComplex::Simplex_handle;
  /** \brief Type for the value of the filtration function. */
  using Filtration_value = typename FilteredComplex::Filtration_value;
  /** \brief Type of element of the field. */
  using Arith_element = typename CoefficientField::Element;
  /** \brief Type for birth and death FilteredComplex::Simplex_handle.
   * The Arith_element field is used for the multi-field framework. */
  using Persistent_interval = std::tuple<Simplex_handle, Simplex_handle, Arith_element>;

  /** \brief Initializes the Persistent_cohomology class.
   *
   * @param[in] cpx Complex for which the persistent homology is computed.
   * cpx is a model of FilteredComplex
   *
   * @param[in] persistence_dim_max if true, the persistent homology for the maximal dimension in the
   *                                complex is computed. If false, it is ignored. Default is false.
   *
   * @exception std::out_of_range In case the number of simplices is more than Simplex_key type numeric limit.
   */
  explicit Persistent_cohomology(FilteredComplex& cpx, bool persistence_dim_max = false)
      : cpx_(&cpx),
        key_handle_map_(cpx_),
        dim_max_(cpx.dimension()),  // upper bound on the dimension of the simplices
        base_(&key_handle_map_, cpx.num_simplices(), cpx.num_vertices()),
        persistent_pairs_(),
        interval_length_policy(&cpx, 0) {
    if (persistence_dim_max) {
      ++dim_max_;
    }
  }

 private:
  struct length_interval {
    length_interval(FilteredComplex* cpx, Filtration_value min_length) : cpx_(cpx), min_length_(min_length) {}

    bool operator()(Simplex_handle sh1, Simplex_handle sh2) {
      return cpx_->filtration(sh2) - cpx_->filtration(sh1) > min_length_;
    }

    void set_length(Filtration_value new_length) { min_length_ = new_length; }

    FilteredComplex* cpx_;
    Filtration_value min_length_;
  };

  /*
   * Compare two intervals by length.
   */
  struct cmp_intervals_by_length {
    explicit cmp_intervals_by_length(FilteredComplex* sc) : sc_(sc) {}
    bool operator()(const Persistent_interval& p1, const Persistent_interval& p2) {
      return (sc_->filtration(get<1>(p1)) - sc_->filtration(get<0>(p1)) >
              sc_->filtration(get<1>(p2)) - sc_->filtration(get<0>(p2)));
    }
    FilteredComplex* sc_;
  };

  struct Key_handle_map {
    using Cell_key = Simplex_key;
    using Cell_handle = Simplex_handle;

    explicit Key_handle_map(FilteredComplex* sc) : sc_(sc) {}

    Cell_handle cell(Cell_key key) const { return sc_->simplex(key); }

    Cell_handle null_cell() const { return sc_->null_simplex(); }

    void assign_key(Cell_handle sh, Cell_key key) { sc_->assign_key(sh, key); }

    Cell_key key(Cell_handle sh) const { return sc_->key(sh); }

    Cell_key null_key() const { return sc_->null_key(); }

    FilteredComplex* sc_;
  };

 public:
  /** \brief Initializes the coefficient field.*/
  void init_coefficients(int charac) { base_.init_coefficients(charac); }
  /** \brief Initializes the coefficient field for multi-field persistent homology.*/
  void init_coefficients(int charac_min, int charac_max) { base_.init_coefficients(charac_min, charac_max); }

  /** \brief Compute the persistent homology of the filtered simplicial
   * complex.
   *
   * @param[in] min_interval_length the computation discards all intervals of length
   *                                less or equal than min_interval_length
   *
   * Assumes that the filtration provided by the simplicial complex is
   * valid. Undefined behavior otherwise. */
  void compute_persistent_cohomology(Filtration_value min_interval_length = 0) {
    if (dim_max_ <= 0) return;  // --------->>

    interval_length_policy.set_length(min_interval_length);
    Simplex_key idx_fil = -1;
    for (auto sh : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(sh, ++idx_fil);
      int dim_simplex = cpx_->dimension(sh);
      base_.update_cohomology_groups(
          sh, cpx_->boundary_simplex_range(sh), dim_max_ == dim_simplex,
          [this](Simplex_handle birth, Simplex_handle death, const Arith_element& characteristic) {
            if (interval_length_policy(birth, death)) {
              persistent_pairs_.emplace_back(birth, death, characteristic);
            }
          });
    }
    base_.finalize_computation([this](Simplex_handle birth, Simplex_handle death, const Arith_element& characteristic) {
      persistent_pairs_.emplace_back(birth, death, characteristic);
    });
  }

  /**
   * @private
   * Temporary patch to make the computation of more general FilteredComplex possible
   * Will be removed when a better solution was agreed on.
   */
  void compute_persistent_cohomology_without_optimizations(Filtration_value min_interval_length = 0) {
    if (dim_max_ <= 0) return;  // --------->>

    interval_length_policy.set_length(min_interval_length);
    Simplex_key idx_fil = -1;
    for (auto sh : cpx_->filtration_simplex_range()) {
      cpx_->assign_key(sh, ++idx_fil);
      int dim_simplex = cpx_->dimension(sh);
      base_.Base::update_cohomology_groups(
          sh, cpx_->boundary_simplex_range(sh), dim_max_ == dim_simplex,
          [this](Simplex_handle birth, Simplex_handle death, const Arith_element& characteristic) {
            if (interval_length_policy(birth, death)) {
              persistent_pairs_.emplace_back(birth, death, characteristic);
            }
          });
    }
    base_.Base::finalize_computation([this](Simplex_handle birth, Simplex_handle death, const Arith_element& characteristic) {
      persistent_pairs_.emplace_back(birth, death, characteristic);
    });
  }

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
    for (auto pair : persistent_pairs_) {
      ostream << get<2>(pair) << "  " << cpx_->dimension(get<0>(pair)) << " " << cpx_->filtration(get<0>(pair)) << " "
              << cpx_->filtration(get<1>(pair)) << " " << std::endl;
    }
  }

  void write_output_diagram(const std::string& diagram_name) {
    std::ofstream diagram_out(diagram_name.c_str());
    diagram_out.exceptions(std::ofstream::failbit);
    cmp_intervals_by_length cmp(cpx_);
    std::sort(std::begin(persistent_pairs_), std::end(persistent_pairs_), cmp);
    for (auto pair : persistent_pairs_) {
      diagram_out << cpx_->dimension(get<0>(pair)) << " " << cpx_->filtration(get<0>(pair)) << " "
                  << cpx_->filtration(get<1>(pair)) << std::endl;
    }
  }

  /** @brief Returns Betti numbers.
   * @return A vector of Betti numbers.
   */
  [[nodiscard]] std::vector<int> betti_numbers() const {
    // Init Betti numbers vector with zeros until Simplicial complex dimension and don't allocate a vector of negative
    // size for an empty complex
    std::vector<int> betti_numbers(std::max(dim_max_, 0));

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
  [[nodiscard]] int betti_number(int dimension) const {
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
    // Init Betti numbers vector with zeros until Simplicial complex dimension and don't allocate a vector of negative
    // size for an empty complex
    std::vector<int> betti_numbers(std::max(dim_max_, 0));
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

  /** @brief Returns a list of persistence birth and death FilteredComplex::Simplex_handle pairs.
   * @return A list of Persistent_cohomology::Persistent_interval
   */
  const std::vector<Persistent_interval>& get_persistent_pairs() const { return persistent_pairs_; }

  /** @brief Returns persistence intervals for a given dimension.
   * @param[in] dimension Dimension to get the birth and death pairs from.
   * @return A vector of persistence intervals (birth and death) on a fixed dimension.
   */
  std::vector<std::pair<Filtration_value, Filtration_value>> intervals_in_dimension(int dimension) {
    std::vector<std::pair<Filtration_value, Filtration_value>> result;
    // auto && pair, to avoid unnecessary copying
    for (auto&& pair : persistent_pairs_) {
      if (cpx_->dimension(get<0>(pair)) == dimension) {
        result.emplace_back(cpx_->filtration(get<0>(pair)), cpx_->filtration(get<1>(pair)));
      }
    }
    return result;
  }

 private:
  FilteredComplex* cpx_;
  Key_handle_map key_handle_map_;
  int dim_max_;
  Polytopal_persistent_cohomology_base<CoefficientField, Key_handle_map> base_;

  /* Persistent intervals. */
  std::vector<Persistent_interval> persistent_pairs_;
  length_interval interval_length_policy;
};

}  // namespace persistent_cohomology

}  // namespace Gudhi

#endif  // GUDHI_PERSISTENT_COHOMOLOGY_H_
