/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FREUDENTHAL_TRIANGULATION_H_
#define FREUDENTHAL_TRIANGULATION_H_

#include <vector>
#include <algorithm>  // for std::sort
#include <cmath>      // for std::floor
#include <numeric>    // for std::iota
#include <cstdlib>    // for std::size_t

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include <gudhi/Permutahedral_representation.h>
#include <gudhi/Debug_utils.h>  // for GUDHI_CHECK

namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Freudenthal_triangulation
 * \brief A class that stores any affine transformation of the Freudenthal-Kuhn
 * triangulation.
 *
 * \ingroup coxeter_triangulation
 *
 * \details The data structure is a record that consists of a matrix
 * that represents the linear transformation of the Freudenthal-Kuhn triangulation
 * and a vector that represents the offset.
 *
 * \tparam Permutahedral_representation_ Type of a simplex given by a permutahedral representation.
 * Needs to be a model of SimplexInCoxeterTriangulation.
 */
template <class Permutahedral_representation_ =
              Permutahedral_representation<std::vector<int>, std::vector<std::vector<std::size_t> > > >
class Freudenthal_triangulation {
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

 public:
  /** \brief Type of the simplices in the triangulation. */
  using Simplex_handle = Permutahedral_representation_;

  /** \brief Type of the vertices in the triangulation. */
  using Vertex_handle = typename Permutahedral_representation_::Vertex;

  /** \brief Constructor of the Freudenthal-Kuhn triangulation of a given dimension.
   * @param[in] dimension The dimension of the triangulation.
   */
  Freudenthal_triangulation(std::size_t dimension)
      : Freudenthal_triangulation(dimension, Matrix::Identity(dimension, dimension), Vector::Zero(dimension)) {
    is_freudenthal_ = true;
  }

  /** \brief Constructor of the Freudenthal-Kuhn triangulation of a given dimension under
   * a linear transformation by a given matrix.
   * @param[in] dimension The dimension of the triangulation.
   * @param[in] matrix The matrix that defines the linear transformation.
   * Needs to be invertible.
   */
  Freudenthal_triangulation(std::size_t dimension, const Matrix& matrix)
      : Freudenthal_triangulation(dimension, matrix, Vector::Zero(dimension)) {}

  /** \brief Constructor of the Freudenthal-Kuhn triangulation of a given dimension under
   * an affine transformation by a given matrix and a translation vector.
   * @param[in] dimension The dimension of the triangulation.
   * @param[in] matrix The matrix that defines the linear transformation.
   * Needs to be invertible.
   * @param[in] offset The offset vector.
   * 
   * @exception std::invalid_argument In debug mode, if offset size is different from dimension.
   */
  Freudenthal_triangulation(unsigned dimension, const Matrix& matrix, const Vector& offset)
      : dimension_(dimension),
        matrix_(matrix),
        offset_(offset),
        colpivhouseholderqr_(matrix_.colPivHouseholderQr()),
        is_freudenthal_(false) {
    GUDHI_CHECK(dimension == offset_.size(), std::invalid_argument("Offset must be of size 'dimension'"));
  }

  /** \brief Dimension of the triangulation. */
  unsigned dimension() const { return dimension_; }

  /** \brief Matrix that defines the linear transformation of the triangulation. */
  const Matrix& matrix() const { return matrix_; }

  /** \brief Vector that defines the offset of the triangulation. */
  const Vector& offset() const { return offset_; }

  /** \brief Change the linear transformation matrix to a given value.
   *  @param[in] matrix New value of the linear transformation matrix.
   */
  void change_matrix(const Eigen::MatrixXd& matrix) {
    matrix_ = matrix;
    colpivhouseholderqr_ = matrix.colPivHouseholderQr();
    is_freudenthal_ = false;
  }

  /** \brief Change the offset vector to a given value.
   *  @param[in] offset New value of the offset vector.
   */
  void change_offset(const Eigen::VectorXd& offset) {
    offset_ = offset;
    is_freudenthal_ = false;
  }

  /** \brief Returns the permutahedral representation of the simplex in the
   *  triangulation that contains a given query point.
   * \details Using the additional parameter scale, the search can be done in a
   * triangulation that shares the origin, but is scaled by a given factor.
   * This parameter can be useful to simulate the point location in a subdivided
   * triangulation.
   * The returned simplex is always minimal by inclusion.
   *
   * \tparam Point_d A class that represents a point in d-dimensional Euclidean space.
   * The coordinates should be random-accessible. Needs to provide the method size().
   *
   * @param[in] point The query point.
   * @param[in] scale The scale of the triangulation.
   * 
   * @exception std::invalid_argument In debug mode, if point dimension is different from triangulation one.
   */
  template <class Point_d>
  Simplex_handle locate_point(const Point_d& point, double scale = 1) const {
    using Ordered_set_partition = typename Simplex_handle::OrderedSetPartition;
    using Part = typename Ordered_set_partition::value_type;
    unsigned d = point.size();
    GUDHI_CHECK(d == dimension_,
                std::invalid_argument("The point must be of the same dimension as the triangulation"));
    double error = 1e-9;
    Simplex_handle output;
    std::vector<double> z;
    if (is_freudenthal_) {
      for (std::size_t i = 0; i < d; i++) {
        double x_i = scale * point[i];
        int y_i = std::floor(x_i);
        output.vertex().push_back(y_i);
        z.push_back(x_i - y_i);
      }
    } else {
      Eigen::VectorXd p_vect(d);
      for (std::size_t i = 0; i < d; i++) p_vect(i) = point[i];
      Eigen::VectorXd x_vect = colpivhouseholderqr_.solve(p_vect - offset_);
      for (std::size_t i = 0; i < d; i++) {
        double x_i = scale * x_vect(i);
        int y_i = std::floor(x_i);
        output.vertex().push_back(y_i);
        z.push_back(x_i - y_i);
      }
    }
    z.push_back(0);
    Part indices(d + 1);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&z](std::size_t i1, std::size_t i2) { return z[i1] > z[i2]; });

    output.partition().push_back(Part(1, indices[0]));
    for (std::size_t i = 1; i <= d; ++i)
      if (z[indices[i - 1]] > z[indices[i]] + error)
        output.partition().push_back(Part(1, indices[i]));
      else
        output.partition().back().push_back(indices[i]);
    return output;
  }

  /** \brief Returns the Cartesian coordinates of the given vertex.
   * \details Using the additional parameter scale, the search can be done in a
   * triangulation that shares the origin, but is scaled by a given factor.
   * This parameter can be useful to simulate the computation of Cartesian coordinates
   * of a vertex in a subdivided triangulation.
   * @param[in] vertex The query vertex.
   * @param[in] scale The scale of the triangulation.
   */
  Eigen::VectorXd cartesian_coordinates(const Vertex_handle& vertex, double scale = 1) const {
    Eigen::VectorXd v_vect(dimension_);
    for (std::size_t j = 0; j < dimension_; j++) v_vect(j) = vertex[j] / scale;
    return matrix_ * v_vect + offset_;
  }

  /** \brief Returns the Cartesian coordinates of the barycenter of a given simplex.
   * \details Using the additional parameter scale, the search can be done in a
   * triangulation that shares the origin, but is scaled by a given factor.
   * This parameter can be useful to simulate the computation of Cartesian coordinates
   * of the barycenter of a simplex in a subdivided triangulation.
   * @param[in] simplex The query simplex.
   * @param[in] scale The scale of the triangulation.
   */
  Eigen::VectorXd barycenter(const Simplex_handle& simplex, double scale = 1) const {
    Eigen::VectorXd res_vector(dimension_);
    res_vector.setZero(dimension_, 1);
    for (auto v : simplex.vertex_range()) {
      res_vector += cartesian_coordinates(v, scale);
    }
    return (1. / (simplex.dimension() + 1)) * res_vector;
  }

 protected:
  unsigned dimension_;
  Matrix matrix_;
  Vector offset_;
  Eigen::ColPivHouseholderQR<Matrix> colpivhouseholderqr_;
  bool is_freudenthal_;
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif
