#ifndef CHAIR_IN_R3_H_
#define CHAIR_IN_R3_H_

#include <cstdlib>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>



struct Function_chair_in_R3 {
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef typename Kernel::Point_d Point_d;

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0)-off_[0], y = p(1)-off_[1], z = p(2)-off_[2];
    Eigen::VectorXd coords(cod_d());
    coords(0) = std::pow(x*x + y*y + z*z - a_*k_*k_, 2) - b_*((z-k_)*(z-k_) - 2*x*x)*((z+k_)*(z+k_) - 2*y*y);
    return coords;
  }

  std::size_t amb_d() const {return 3;}
  std::size_t cod_d() const {return 1;}
  Eigen::VectorXd seed() const {
    // double z0 = std::sqrt(k_/(1.0-b_)) * std::sqrt(a_-b_ + std::sqrt((a_-b_)*(a_-b_) - (1.0-b_)*(a_*a_ - b_)*k_*k_));
    double t1 = a_-b_;
    double det = t1*t1 - (1.0 - b_)*(a_*a_ - b_);
    double z0 = k_*std::sqrt((t1+std::sqrt(det))/(1-b_));
    return Eigen::Vector3d(off_[0], off_[1], z0+off_[2]);
  }
  
  Function_chair_in_R3(double a = 0.8, double b = 0.4, double k = 1.0) :
    a_(a), b_(b), k_(k), 
    off_(*CGAL::Random_points_on_sphere_d<Point_d>(3, 0.001)) {}
  double a_, b_, k_;
  Point_d off_; // random perturbation of the center
};

#endif

// (x^2 + y^2 + z^2 - a*k^2)^2 - b*((z-k)^2 - 2*x^2)*((z+k)^2 - 2*y^2)
// sqrt(k/(1-b))*sqrt(a-b + sqrt((a-b)^2 - (1-b)*(a^2 - b)*k^2))
