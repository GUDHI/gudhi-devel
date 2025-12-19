/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Persistence_interface_matrix.h
 * @author Hannah Schreiber
 * @brief Contains the @ref Gudhi::multi_persistence::Persistence_interface_matrix class.
 */

#ifndef MP_PERSISTENCE_INTERFACE_MATRIX_H_INCLUDED
#define MP_PERSISTENCE_INTERFACE_MATRIX_H_INCLUDED

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <boost/range/iterator_range_core.hpp>

#include <gudhi/Matrix.h>
#include <gudhi/Debug_utils.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @class Persistence_interface_matrix Persistence_interface_matrix.h \
 * gudhi/Multi_persistence/Persistence_interface_matrix.h
 * @ingroup multi_persistence
 *
 * @brief Interface respecting the @ref PersistenceAlgorithm concept to use @ref Slicer with the homology, vineyard
 * and representative cycle algorithms implemented in @ref Gudhi::persistence_matrix::Matrix.
 *
 * @tparam PosIdxPersistenceMatrixOptions Options respecting the
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions concept such that, either
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::column_indexation_type "column_indexation_type" is
 * @ref Gudhi::persistence_matrix::Column_indexation_types::POSITION "POSITION", or,
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::is_of_boundary_type "is_of_boundary_type" is true and
 * @ref Gudhi::persistence_matrix::PersistenceMatrixOptions::column_indexation_type "column_indexation_type" is
 * @ref Gudhi::persistence_matrix::Column_indexation_types::CONTAINER "CONTAINER".
 */
template <class PosIdxPersistenceMatrixOptions>
class Persistence_interface_matrix
{
 public:
  using Options = PosIdxPersistenceMatrixOptions;
  using Matrix = Gudhi::persistence_matrix::Matrix<Options>; /**< Complex type */
  using Dimension = typename Options::Dimension;             /**< Dimension type */
  using Index = typename Options::Index;                     /**< Index type */
  using Map = std::vector<Index>;                            /**< Map type */
  using Bar = typename Matrix::Bar;                          /**< Bar type */
  using Cycle = typename Matrix::Cycle;                      /**< Cycle type */
  template<class Complex>
  using As_type = Persistence_interface_matrix<PosIdxPersistenceMatrixOptions>; /**< This type. */

  class Barcode_iterator
      : public boost::iterator_facade<Barcode_iterator, const Bar, boost::random_access_traversal_tag>
  {
   private:
    using Base = boost::iterator_facade<Barcode_iterator, const Bar, boost::random_access_traversal_tag>;

   public:
    using Barcode = typename Matrix::Barcode;

    using difference_type = typename Base::difference_type;
    using size_type = std::size_t;
    using const_reference = const Bar&;

    Barcode_iterator(const Barcode& barcode, const Map& permutation, size_type pos)
        : barcode_(&barcode),
          permutation_(&permutation),
          currPos_(pos > barcode.size() ? barcode.size() : pos),
          currBar_()
    {
      // end otherwise
      if (currPos_ < barcode.size()) _update_bar();
    }

    // necessary for boost::iterator_range to be able to use operator().
    // operator[] not possible in this particular case
    Bar operator[](difference_type n) const
    {
      const auto& b = (*barcode_)[currPos_ + n];  // n is out of range if currPos_ + n is out of range from barcode_
      return Bar(_posToIndex(b.birth), b.death == Bar::inf ? b.death : _posToIndex(b.death), b.dim);
    }

   private:
    using Pos_index = typename std::tuple_element<2, Bar>::type;

    friend class boost::iterator_core_access;

    bool equal(Barcode_iterator const& other) const
    {
      return barcode_ == other.barcode_ && permutation_ == other.permutation_ && currPos_ == other.currPos_;
    }

    const_reference dereference() const { return currBar_; }

    void increment()
    {
      ++currPos_;
      if (currPos_ < barcode_->size()) _update_bar();
    }

    void decrement()
    {
      --currPos_;
      if (currPos_ < barcode_->size()) _update_bar();
    }

    void advance(difference_type n)
    {
      currPos_ += n;
      if (currPos_ < barcode_->size()) _update_bar();
    }

    difference_type distance_to(const Barcode_iterator& other) const { return other.currPos_ - currPos_; }

    Pos_index _posToIndex(Pos_index pos) const { return (*permutation_)[pos]; }

    void _update_bar()
    {
      const auto& b = (*barcode_)[currPos_];
      currBar_.dim = b.dim;
      currBar_.birth = _posToIndex(b.birth);
      currBar_.death = b.death == Bar::inf ? b.death : _posToIndex(b.death);
    }

    Barcode const* barcode_;
    Map const* permutation_;
    size_type currPos_;
    Bar currBar_;
  };

  class Cycles_iterator
      : public boost::iterator_facade<Cycles_iterator, const Cycle, boost::random_access_traversal_tag>
  {
   private:
    using Base = boost::iterator_facade<Barcode_iterator, const Bar, boost::random_access_traversal_tag>;

   public:
    using Cycles = std::vector<Cycle>;

    using difference_type = typename Base::difference_type;
    using size_type = std::size_t;
    using const_reference = const Cycle&;

    Cycles_iterator(const Cycles& cycles, const Map& permutation, Map const* idToPos, size_type pos)
        : cycles_(&cycles),
          permutation_(&permutation),
          idToPos_(idToPos),
          currPos_(pos > cycles.size() ? cycles.size() : pos),
          currCycle_()
    {
      if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
        GUDHI_CHECK(idToPos_ != nullptr, "ID to position map has to be set for chain matrices using vineyard.");
      }
      // end otherwise
      if (currPos_ < cycles.size()) _update_cycle();
    }

    // necessary for boost::iterator_range to be able to use operator().
    // operator[] not possible in this particular case
    Cycle operator[](difference_type n) const
    {
      Cycle res((*cycles_)[currPos_ + n]);
      for (auto& id : res) {
        if constexpr (PosIdxPersistenceMatrixOptions::is_z2) {
          id = _id_to_index(id);
        } else {
          id.first = _id_to_index(id.first);
        }
      }
      return res;
    }

   private:
    friend class boost::iterator_core_access;

    bool equal(Cycles_iterator const& other) const
    {
      return cycles_ == other.cycles_ && permutation_ == other.permutation_ && currPos_ == other.currPos_;
    }

    const_reference dereference() const { return currCycle_; }

    void increment()
    {
      ++currPos_;
      if (currPos_ < cycles_->size()) _update_cycle();
    }

    void decrement()
    {
      --currPos_;
      if (currPos_ < cycles_->size()) _update_cycle();
    }

    void advance(difference_type n)
    {
      currPos_ += n;
      if (currPos_ < cycles_->size()) _update_cycle();
    }

    difference_type distance_to(const Cycles_iterator& other) const { return other.currPos_ - currPos_; }

    Index _id_to_index(Index id) const
    {
      if constexpr (Options::is_of_boundary_type || !Options::has_vine_update) {
        // works for RU because id == pos, but does not work for chain with vine
        // we need a id to pos map in that case
        return (*permutation_)[id];
      } else {
        return (*permutation_)[(*idToPos_)[id]];
      }
    }

    void _update_cycle()
    {
      const auto& c = (*cycles_)[currPos_];
      currCycle_.resize(c.size());
      for (size_type i = 0; i < c.size(); ++i) {
        if constexpr (PosIdxPersistenceMatrixOptions::is_z2) {
          currCycle_[i] = _id_to_index(c[i]);
        } else {
          currCycle_[i].first = _id_to_index(c[i]);
        }
      }
    }

    Cycles const* cycles_;
    Map const* permutation_;
    Map const* idToPos_;
    size_type currPos_;
    Cycle currCycle_;
  };

  using Barcode = boost::iterator_range<Barcode_iterator>; /**< Barcode type */
  using Cycles = boost::iterator_range<Cycles_iterator>;   /**< Cycle container type */

  static constexpr const auto nullDeath = Bar::inf;
  /**
   * @brief True if and only if PosIdxPersistenceMatrixOptions::has_vine_update is true.
   */
  static constexpr const bool is_vine = Options::has_vine_update;
  /**
   * @brief True if and only if PosIdxPersistenceMatrixOptions::can_retrieve_representative_cycles is true.
   */
  static constexpr const bool has_rep_cycles = Options::can_retrieve_representative_cycles;

  Persistence_interface_matrix() : permutation_(nullptr) {}

  // `permutation` is assumed to have stable size, i.e., its address never changes
  template <class Complex>
  Persistence_interface_matrix(const Complex& cpx, const Map& permutation)
      : matrix_(permutation.size()), permutation_(&permutation)
  {
    static_assert(
        Options::column_indexation_type == Gudhi::persistence_matrix::Column_indexation_types::POSITION ||
            (Options::is_of_boundary_type &&
             Options::column_indexation_type == Gudhi::persistence_matrix::Column_indexation_types::CONTAINER),
        "Matrix has a non supported index scheme.");

    _initialize(cpx);
  }

  Persistence_interface_matrix(const Persistence_interface_matrix& other) = delete;

  // permutation is assumed to be the same than from the copied object, just its address can change
  Persistence_interface_matrix(const Persistence_interface_matrix& other, const Map& permutation)
      : matrix_(other.matrix_), permutation_(other.is_initialized() ? &permutation : nullptr), idToPos_(other.idToPos_)
  {
    GUDHI_CHECK(!other.is_initialized() || permutation == *other.permutation_,
                "Only the address of the permutation vector is allowed to change, not its content.");
  }

  Persistence_interface_matrix(Persistence_interface_matrix&& other) = delete;

  // permutation is assumed to be the same than from the moved object, just its address can change
  Persistence_interface_matrix(Persistence_interface_matrix&& other, const Map& permutation)
      : matrix_(std::move(other.matrix_)),
        permutation_(other.is_initialized() ? &permutation : nullptr),
        idToPos_(std::move(other.idToPos_))
  {
    other.permutation_ = nullptr;
  }

  ~Persistence_interface_matrix() = default;

  Persistence_interface_matrix& operator=(const Persistence_interface_matrix& other) = delete;
  Persistence_interface_matrix& operator=(Persistence_interface_matrix&& other) noexcept = delete;

  // TODO: swap?

  template <class Complex>
  void reinitialize(const Complex& cpx, const Map& permutation)
  {
    matrix_ = Matrix(permutation.size());
    permutation_ = &permutation;
    _initialize(cpx);
  }

  void reset()
  {
    matrix_ = Matrix();
    permutation_ = nullptr;
    if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
      if (idToPos_) idToPos_->clear();
    }
  }

  [[nodiscard]] bool is_initialized() const { return permutation_ != nullptr; }

  Dimension get_dimension(Index i) const
  {
    GUDHI_CHECK(is_initialized(), "Dimension can not be computed uninitialized.");
    return matrix_.get_column_dimension(i);
  }

  Barcode get_barcode()
  {
    GUDHI_CHECK(is_initialized(), "Barcode can not be computed uninitialized.");

    const auto& barcode = matrix_.get_current_barcode();
    return Barcode(Barcode_iterator(barcode, *permutation_, 0),
                   Barcode_iterator(barcode, *permutation_, barcode.size()));
  }

  void vine_swap(Index i)
  {
    static_assert(is_vine, "`vine_swap` is not enabled with the given options.");
    GUDHI_CHECK(is_initialized(), "Vineyard can not be computed uninitialized.");

    if constexpr (!Options::is_of_boundary_type) {
      auto id1 = matrix_.get_pivot(i);
      auto id2 = matrix_.get_pivot(i + 1);
      std::swap((*idToPos_)[id1], (*idToPos_)[id2]);
    }
    matrix_.vine_swap(i);
  }

  Cycles get_representative_cycles(bool update)
  {
    static_assert(has_rep_cycles, "`get_representative_cycles` is not enabled with the given options.");
    GUDHI_CHECK(is_initialized(), "Representative cycles can not be computed uninitialized.");

    if (update) matrix_.update_representative_cycles();

    const auto& cycles = matrix_.get_representative_cycles();
    Map* idToPosPtr = nullptr;
    if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
      idToPosPtr = &(*idToPos_);
    }
    return Cycles(Cycles_iterator(cycles, *permutation_, idToPosPtr, 0),
                  Cycles_iterator(cycles, *permutation_, idToPosPtr, cycles.size()));
  }

  Cycle get_representative_cycle(Index barcodeIndex, bool update)
  {
    static_assert(has_rep_cycles, "`get_representative_cycle` is not enabled with the given options.");
    GUDHI_CHECK(is_initialized(), "Representative cycles can not be computed uninitialized.");

    auto id_to_index = [&](Index id) -> Index {
      if constexpr (Options::is_of_boundary_type || !Options::has_vine_update) {
        // works for RU because id == pos, but does not work for chain with vine
        // we need a id to pos map in that case
        return (*permutation_)[id];
      } else {
        return (*permutation_)[(*idToPos_)[id]];
      }
    };

    if (update) matrix_.update_representative_cycles();

    const auto& c = matrix_.get_representative_cycle(matrix_.get_current_barcode()[barcodeIndex]);

    Cycle cycle(c.size());
    for (Index i = 0; i < c.size(); ++i) {
      if constexpr (PosIdxPersistenceMatrixOptions::is_z2) {
        cycle[i] = id_to_index(c[i]);
      } else {
        cycle[i].first = id_to_index(c[i]);
      }
    }
    return cycle;
  }

  /**
   * @brief Outstream operator.
   */
  friend std::ostream& operator<<(std::ostream& stream, Persistence_interface_matrix& pers)
  {
    stream << "Matrix:\n";
    stream << "[\n";
    for (auto i = 0U; i < pers.matrix_.get_number_of_columns(); i++) {
      stream << "[";
      for (const auto& v : pers.matrix_.get_column(i)) stream << v << ", ";
      stream << "]\n";
    }
    stream << "]\n";
    stream << "Permutation:\n";
    if (pers.permutation_ != nullptr) {
      for (auto v : *pers.permutation_) {
        stream << v << " ";
      }
    }
    stream << "\n";
    if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
      stream << "ID to position map:\n";
      for (auto v : *pers.idToPos_) {
        stream << v << " ";
      }
      stream << "\n";
    }

    return stream;
  }

 private:
  Matrix matrix_;
  Map const* permutation_;
  std::optional<Map> idToPos_;

  template <class Complex>
  void _initialize(const Complex& cpx)
  {
    if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
      idToPos_.emplace();
      idToPos_->reserve(permutation_->size());
    }
    const auto& boundaries = cpx.get_boundaries();
    const auto& dimensions = cpx.get_dimensions();
    // simplex IDs need to be increasing in order, so the original ones cannot be used
    Map permutationInv(cpx.get_number_of_cycle_generators());
    Map translated_boundary;
    std::size_t id = 0;
    for (auto i : *permutation_) {
      permutationInv[i] = id;  // permutation is assumed to be a valid filtration
      translated_boundary.resize(boundaries[i].size());
      for (std::size_t j = 0; j < boundaries[i].size(); ++j) {
        translated_boundary[j] = permutationInv[boundaries[i][j]];
      }
      std::sort(translated_boundary.begin(), translated_boundary.end());
      matrix_.insert_boundary(id, translated_boundary, dimensions[i]);
      if constexpr (Options::has_vine_update && !Options::is_of_boundary_type) {
        idToPos_->push_back(id);
      }
      ++id;  // IDs corresponds to the indices in permutation_
    }
  }
};

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_PERSISTENCE_INTERFACE_MATRIX_H_INCLUDED
