/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef FUNCTIONS_FUNCTION_AFFINE_PLANE_IN_RD_H_
#define FUNCTIONS_FUNCTION_AFFINE_PLANE_IN_RD_H_

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** 
 * \class Function_affine_plane_in_Rd 
 * \brief A class for the function that defines an m-dimensional implicit affine plane 
 * embedded in d-dimensional Euclidean space.
 *
 * \ingroup coxeter_triangulation
 */
struct Function_affine_plane_in_Rd {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd x = p - off_;
    return normal_matrix_.transpose() * x;
  }
  
  /** \brief Returns the domain dimension. Same as the ambient dimension of the sphere. */
  std::size_t amb_d() const {return d_;};

  /** \brief Returns the codomain dimension. Same as the codimension of the sphere. */
  std::size_t cod_d() const {return k_;};

  /** \brief Returns a point on the affine plane. */
  Eigen::VectorXd seed() const {
    return off_;
  }
  
  /** 
   * \brief Constructor of the function that defines an m-dimensional implicit affine
   * plane in the d-dimensional Euclidean space.
   *
   * @param[in] normal_matrix A normal matrix of the affine plane. The number of rows should
   * correspond to the ambient dimension, the number of columns should corespond to 
   * the size of the normal basis (codimension) 
   * @param[in] offset The offset vector of the affine plane.
   */
  Function_affine_plane_in_Rd(const Eigen::MatrixXd& normal_matrix,
			      const Eigen::VectorXd& offset)
    : normal_matrix_(normal_matrix),
      d_(normal_matrix.rows()),
      k_(normal_matrix.cols()),
      m_(d_ - k_),
      off_(offset) {
    normal_matrix_.colwise().normalize();
  }

  /** 
   * \brief Constructor of the function that defines an m-dimensional implicit affine
   * plane in the d-dimensional Euclidean space that passes through origin.
   *
   * @param[in] normal A normal matrix of the affine plane. The number of rows should
   * correspond to the ambient dimension, the number of columns should corespond to 
   * the size of the normal basis (codimension) 
   */
  Function_affine_plane_in_Rd(const Eigen::MatrixXd& normal_matrix)
    : normal_matrix_(normal_matrix),
      d_(normal_matrix.rows()),
      k_(normal_matrix.cols()),
      m_(d_ - k_),
      off_(Eigen::VectorXd::Zero(d_)) {
    normal_matrix_.colwise().normalize();
  }

private:
  Eigen::MatrixXd normal_matrix_;
  std::size_t d_, k_, m_;
  Eigen::VectorXd off_;  
};

} // namespace coxeter_triangulation

} // namespace Gudhi


#endif
