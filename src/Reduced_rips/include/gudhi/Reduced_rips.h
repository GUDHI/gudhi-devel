/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Thomas Burnett, Musashi Koyama
 *
 *    Copyright (C) 2026 Thomas Burnett, Musashi Koyama
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file Reduced_rips.h
 * @author Thomas Burnett, Musashi Koyama
 * @brief The public Reduced_rips class: degree-1 Vietoris-Rips persistence from a point cloud or a distance
 * matrix, via the Reduced Vietoris-Rips filtration.
 */

#ifndef REDUCED_RIPS_H_
#define REDUCED_RIPS_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Reduced_rips/Euclidean_geometry.h>
#include <gudhi/Reduced_rips/Euclidean_kd_tree.h>
#include <gudhi/Reduced_rips/Helpers.h>
#include <gudhi/Reduced_rips/Matrix_geometry.h>
#include <gudhi/Reduced_rips/Persistence_engine.h>

namespace Gudhi {

namespace reduced_rips {

/**
 * @class Reduced_rips Reduced_rips.h gudhi/Reduced_rips.h
 * @brief Degree-1 Vietoris-Rips persistent homology, computed via the Reduced Vietoris-Rips filtration, from
 * either a Euclidean point cloud or a bare symmetric distance matrix.
 *
 * @ingroup reduced_rips
 *
 * @details
 * Rather than building the full Vietoris-Rips complex, this class computes the relative neighborhood graph
 * and reduces only the relevant 2-simplices (those certified by the lune construction of the reference
 * paper). This scales to far larger inputs than building the full complex would, though only degree 1 is
 * computed.
 *
 * Two inputs are accepted, via two named factories (@ref from_points and @ref from_distance_matrix) that share one
 * metric reduction core:
 * - a **Euclidean point cloud** (the primary constructor): the relative neighborhood graph is computed with
 *   a Delaunay triangulation in ambient dimension 2 and 3 (selected at run time from the point dimension)
 *   and a direct O(n^2) construction in higher dimension; neighbor queries use a kd-tree or a brute-force
 *   scan (see #Search). The lens-ball / wide-angle certificates of the reference paper accelerate the lune
 *   computation;
 * - an **arbitrary symmetric distance matrix** (@ref from_distance_matrix): the same reduction, driven purely by
 *   the supplied distances. No coordinates are available, so the relative neighborhood graph uses the
 *   dimension-free O(n^2) construction, neighbor queries scan matrix rows, and the lune connected components
 *   are found by the exact union-find.
 *
 * The reduced filtration is exact for any symmetric matrix of non-negative dissimilarities.
 *
 * @tparam Filtration_value_ Arithmetic type used for the distance/filtration arithmetic and the birth/death
 * values of the output barcode (`double` by default). Point coordinates are still stored as `double`, and
 * the Euclidean path's spatial acceleration (the CGAL kd-tree and 2D/3D Delaunay) searches in `double`; but
 * every candidate it returns is re-tested with an exact distance in `Filtration_value_`, so the filtration
 * values themselves carry this type's precision.
 *
 * @tparam Index_ Unsigned integer type used internally to store point indices and 1-simplex ids (`std::uint32_t`
 * by default). It must be wide enough to number the 1-simplices, which can far exceed the point count (the worst case
 * is n(n-1)/2). The default allows ~4.29 billion of them while halving the reduction's neighbor lists and columns
 * versus a 64-bit index; the factories throw `std::invalid_argument` if the point count itself does not fit, and the
 * computation throws `std::overflow_error` in the event that the processed edges outgrow it.
 */
template <typename Filtration_value_ = double, typename Index_ = std::uint32_t>
class Reduced_rips {
 public:
  /** @brief Type used to store filtration / persistence values. */
  using Filtration_value = Filtration_value_;
  /** @brief Unsigned integer type used to store point indices and 1-simplex ids internally. */
  using Index = Index_;
  /** @brief A persistence bar as a `{birth, death}` array of distances (`bar[0]` birth, `bar[1]` death) */
  using Persistence_interval = std::array<Filtration_value, 2>;

  /** @brief Strategy for the spatial neighbor queries (k-nearest and radius search). Applies to the
   * Euclidean point-cloud constructor only; the distance-matrix path always scans matrix rows.
   *
   * - `kd_tree`: a CGAL Epick_d kd-tree. Its reach depends on the *ambient* dimension; the right choice in
   *   low dimension, where the tree can prune well.
   * - `brute_force`: a flat O(n) scan per query. No traversal overhead, and faster than the kd-tree once its
   *   pruning collapses in high ambient dimension, at the cost of O(n^2) total search work.
   * - `automatic`: pick `kd_tree` for ambient dimension <= 3, `brute_force` otherwise.
   */
  enum class Search : std::uint8_t { automatic, kd_tree, brute_force };

  /** @brief Constructs an empty diagram. Assign the result of a @ref from_points or @ref from_distance_matrix factory
   * to populate it. */
  Reduced_rips() = default;

  /** @brief Builds the degree-1 persistence from a range of points. The barcode is computed eagerly and then
   * returned by #persistence().
   *
   * @tparam PointRange A forward range whose elements are themselves forward ranges of coordinates convertible
   * to `double` (e.g. `std::vector<std::vector<double>>`). All points must share the same dimension. The range
   * is read once; nothing is retained by reference.
   *
   * @param[in] points Range of points, as above.
   * @param[in] num_neighbors Initial neighbor budget per point for the heap seeding. Pass 0 (the default) to
   * use `sqrt(n)`.
   * @param[in] search Spatial-search strategy (see #Search); `automatic` by default.
   *
   * @exception std::invalid_argument In debug mode, if points have differing dimension.
   */
  template <typename PointRange>
  static Reduced_rips from_points(const PointRange& points, unsigned int num_neighbors = 0,
                                  Search search = Search::automatic) {
    static_assert(std::is_floating_point_v<Filtration_value>,
                  "Reduced_rips::from_points requires a floating-point Filtration_value. Use from_distance_matrix "
                  "for exact scalar types.");
    Reduced_rips rr;
    rr.num_neighbors_ = num_neighbors;
    rr.search_ = search;
    rr.ingest_points(points);
    if (rr.n_ >= 2) rr.compute_from_points();
    return rr;
  }

  /** @brief Builds the degree-1 persistence from an arbitrary symmetric distance matrix. The barcode is
   * computed eagerly and then returned by #persistence().
   *
   * Accepts either a full symmetric matrix (`matrix[i][j]` for all `i, j`) or a lower-triangular one
   * (`matrix[i]` of length `i`, holding `matrix[i][j]` for `j < i`); only the entries below the diagonal of
   * each row are read, so the two layouts are handled uniformly. The diagonal is taken to be zero and entries
   * must be non-negative dissimilarities.
   *
   * @tparam DistanceMatrix A forward range of rows, each row a forward range of distances convertible to
   * `double`, laid out as above. The range is read once; nothing is retained by reference.
   *
   * @param[in] matrix Range of ranges of distances, as above.
   * @param[in] num_neighbors Initial neighbor budget per point for the heap seeding. Pass 0 (the default) to
   * use `sqrt(n)`.
   * @param[in] search Ignored (matrix neighbor queries always scan rows).
   *
   * @exception std::invalid_argument In debug mode, if a row is too short to supply its lower-triangle
   * distances, or if a distance is negative (or NaN).
   */
  template <typename DistanceMatrix>
  static Reduced_rips from_distance_matrix(const DistanceMatrix& matrix, unsigned int num_neighbors = 0,
                                           Search search = Search::automatic) {
    Reduced_rips rr;
    rr.num_neighbors_ = num_neighbors;
    rr.search_ = search;
    rr.ingest_distance_matrix(matrix);
    if (rr.n_ >= 2) rr.compute_from_matrix();
    return rr;
  }

  /** @brief Returns the degree-1 persistence barcode, as (birth, death) pairs of distances in ascending order
   * of death (the order in which their cycles are filled). Computed eagerly at construction.
   */
  [[nodiscard]] const std::vector<Persistence_interval>& persistence() const { return barcode_; }

  /** @brief Ambient dimension of the input points, or 0 for the distance-matrix input (no coordinates). */
  [[nodiscard]] std::size_t dimension() const { return dim_; }

  /** @brief Number of distinct 1-simplices whose lune was evaluated and applied (diagnostic). */
  [[nodiscard]] std::size_t num_one_simplices() const { return num_one_simplices_; }
  /** @brief Number of 2-simplex columns formed during the reduction (diagnostic). */
  [[nodiscard]] std::size_t num_two_simplices() const { return num_two_simplices_; }
  /** @brief Number of recorded (non-apparent) persistent pairs (diagnostic). */
  [[nodiscard]] std::size_t num_persistence_pairs() const { return num_persistence_pairs_; }

 private:
  [[nodiscard]] detail::Cloud cloud() const { return {coords_.data(), dim_, n_}; }

  // Reads a range of points into coords_ as a flat row-major buffer, setting dim_ and n_.
  template <typename PointRange>
  void ingest_points(const PointRange& points) {
    auto it = std::begin(points), end = std::end(points);
    if (it == end) return;  // no points: empty barcode
    dim_ = std::distance(std::begin(*it), std::end(*it));
    if (dim_ == 0) return;                                                     // zero-dimensional points: empty barcode
    coords_.reserve(static_cast<std::size_t>(std::distance(it, end)) * dim_);  // forward range: multipass is fine
    for (; it != end; ++it, ++n_) {
      GUDHI_CHECK_code(std::size_t before = coords_.size());
      coords_.insert(coords_.end(), std::begin(*it), std::end(*it));
      GUDHI_CHECK(coords_.size() - before == dim_,
                  std::invalid_argument("Reduced_rips: all points must share one dimension"));
    }
  }

  // Resolves Search::automatic to a concrete strategy: kd-tree in low ambient dimension, brute force above.
  [[nodiscard]] bool use_brute_force() const {
    if (search_ == Search::kd_tree) return false;
    if (search_ == Search::brute_force) return true;
    return dim_ >= 4;
  }

  // Reads a full or lower-triangular symmetric distance matrix into matrix_ as a flat n by n row-major buffer
  // of the supplied distances, kept as-is. Only the entries below the diagonal of each row are read, matrix[i][j]
  // for j < i, so the full and lower-triangular layouts are handled uniformly in a single forward pass per row;
  // the upper triangle is mirrored and the diagonal is zero.
  template <typename DistanceMatrix>
  void ingest_distance_matrix(const DistanceMatrix& matrix) {
    n_ = std::distance(std::begin(matrix), std::end(matrix));
    if (n_ < 2) return;  // fewer than two points: empty barcode
    matrix_.assign(n_ * n_, Filtration_value(0));
    std::size_t i = 0;
    for (const auto& row : matrix) {
      std::size_t j = 0;
      auto rit = std::begin(row);
      const auto rend = std::end(row);
      for (; j < i && rit != rend; ++j, ++rit) {
        Filtration_value d = *rit;  // lower triangle: distance between i and j
        GUDHI_CHECK(d >= Filtration_value(0), std::invalid_argument("Reduced_rips: distances must be non-negative"));
        matrix_[(i * n_) + j] = d;
        matrix_[(j * n_) + i] = d;
      }
      GUDHI_CHECK(j >= i, std::invalid_argument("Reduced_rips: distance matrix row is too short"));
      ++i;
    }
  }

  // Point-cloud path: build the kd-tree (or brute-force) geometry over the ingested coordinates and compute.
  void compute_from_points() {
    detail::Cloud pm = cloud();
    Euclidean_kd_tree<Index> kd_tree(pm, use_brute_force());
    Euclidean_geometry<Filtration_value, Index> geom(pm, kd_tree);
    compute_impl(geom);
  }

  // Distance-matrix path: build the matrix geometry over the ingested distances and compute.
  void compute_from_matrix() {
    Matrix_geometry<Filtration_value, Index> geom(std::move(matrix_), n_);
    compute_impl(geom);
  }

  template <class Geom>
  void compute_impl(Geom& geom) {
    Persistence_engine<Geom> engine(geom, num_neighbors_);
    engine.run();

    const auto& bars = engine.barcode();
    barcode_.reserve(bars.size());
    // The engine emits birth/death on the geometry's edge-length scale; to_distance maps each back to a real
    // distance (sqrt for the Euclidean policy, the identity for the matrix policy).
    for (const auto& bar : bars) barcode_.push_back({Geom::to_distance(bar.first), Geom::to_distance(bar.second)});
    num_one_simplices_ = engine.committed_edges();
    num_two_simplices_ = engine.columns_formed();
    num_persistence_pairs_ = engine.deaths();
  }

  std::vector<double> coords_;            // flat row-major point storage, n_ points x dim_ coordinates (cloud input)
  std::vector<Filtration_value> matrix_;  // flat n_ x n_ row-major distances (distance-matrix input)
  std::size_t dim_ = 0;
  std::size_t n_ = 0;
  unsigned int num_neighbors_ = 0;
  Search search_ = Search::automatic;
  std::vector<Persistence_interval> barcode_;
  std::size_t num_one_simplices_ = 0;
  std::size_t num_two_simplices_ = 0;
  std::size_t num_persistence_pairs_ = 0;
};

}  // namespace reduced_rips

}  // namespace Gudhi

#endif  // REDUCED_RIPS_H_
