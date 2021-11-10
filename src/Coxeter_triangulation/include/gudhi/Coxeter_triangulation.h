/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef COXETER_TRIANGULATION_H_
#define COXETER_TRIANGULATION_H_

#include <vector>
#include <cmath>  // for std::sqrt

#include <boost/range/iterator_range.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include <gudhi/Freudenthal_triangulation.h>
#include <gudhi/Permutahedral_representation.h>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 * \class Coxeter_triangulation
 * \brief A class that stores Coxeter triangulation of type \f$\tilde{A}_d\f$.
 * This triangulation has the greatest simplex quality out of all linear transformations
 * of the Freudenthal-Kuhn triangulation.
 *
 * \ingroup coxeter_triangulation
 *
 * \tparam Permutahedral_representation_ Type of a simplex given by a permutahedral representation.
 *  Needs to be a model of SimplexInCoxeterTriangulation.
 */
template <class Permutahedral_representation_ =
              Permutahedral_representation<std::vector<int>, std::vector<std::vector<std::size_t> > > >
class Coxeter_triangulation : public Freudenthal_triangulation<Permutahedral_representation_> {
  using Matrix = Eigen::MatrixXd;

  Matrix root_matrix(unsigned d) {
    Matrix cartan(Matrix::Identity(d, d));
    for (unsigned i = 1; i < d; i++) {
      cartan(i - 1, i) = -0.5;
      cartan(i, i - 1) = -0.5;
    }
    Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
    Eigen::VectorXd sqrt_diag(d);
    for (unsigned i = 0; i < d; ++i) sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);

    Matrix lower(Matrix::Ones(d, d));
    lower = lower.triangularView<Eigen::Lower>();

    Matrix result = (lower * saes.eigenvectors() * sqrt_diag.asDiagonal()).inverse();
    return result;
  }

 public:
  /** \brief Constructor of Coxeter triangulation of a given dimension.
   * @param[in] dimension The dimension of the triangulation.
   */
  Coxeter_triangulation(std::size_t dimension)
      : Freudenthal_triangulation<Permutahedral_representation_>(dimension, root_matrix(dimension)) {}
};

}  // namespace coxeter_triangulation

}  // namespace Gudhi

#endif
