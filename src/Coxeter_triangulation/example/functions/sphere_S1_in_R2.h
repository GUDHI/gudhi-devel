#ifndef SPHERE_S1_IN_R2_H_
#define SPHERE_S1_IN_R2_H_

#include <gudhi_patches/CGAL/Epick_d.h>
#include <gudhi/random_point_generators.h> // construct_point

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

struct Function_S1_in_R2 {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef typename Kernel::Point_d Point_d;

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0)-off_[0], y = p(1)-off_[1];
    Eigen::VectorXd coords(cod_d());
    coords(0) = x*x + y*y - r_*r_;
    return coords;
  }

  std::size_t amb_d() const {return 2;};
  std::size_t cod_d() const {return 1;};
  
  Function_S1_in_R2(double r) : r_(r), off_(*CGAL::Random_points_on_sphere_d<Point_d>(3, 0.001)) {}
  double r_;
  Point_d off_; // random perturbation of the center
};

#endif
