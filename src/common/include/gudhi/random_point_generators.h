/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clement Jamin
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef RANDOM_POINT_GENERATORS_H_
#define RANDOM_POINT_GENERATORS_H_

#include <CGAL/number_utils.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <vector>  // for vector<>

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Alpha_complex_3d is only available for CGAL >= 4.11
#endif

namespace Gudhi {

///////////////////////////////////////////////////////////////////////////////
// Note: All these functions have been tested with the CGAL::Epick_d kernel
///////////////////////////////////////////////////////////////////////////////

// construct_point: dim 2

template <typename Kernel>
typename Kernel::Point_d construct_point(const Kernel &k,
                                         typename Kernel::FT x1, typename Kernel::FT x2) {
  typename Kernel::FT tab[2];
  tab[0] = x1;
  tab[1] = x2;
  return k.construct_point_d_object()(2, &tab[0], &tab[2]);
}

// construct_point: dim 3

template <typename Kernel>
typename Kernel::Point_d construct_point(const Kernel &k,
                                         typename Kernel::FT x1, typename Kernel::FT x2, typename Kernel::FT x3) {
  typename Kernel::FT tab[3];
  tab[0] = x1;
  tab[1] = x2;
  tab[2] = x3;
  return k.construct_point_d_object()(3, &tab[0], &tab[3]);
}

// construct_point: dim 4

template <typename Kernel>
typename Kernel::Point_d construct_point(const Kernel &k,
                                         typename Kernel::FT x1, typename Kernel::FT x2, typename Kernel::FT x3,
                                         typename Kernel::FT x4) {
  typename Kernel::FT tab[4];
  tab[0] = x1;
  tab[1] = x2;
  tab[2] = x3;
  tab[3] = x4;
  return k.construct_point_d_object()(4, &tab[0], &tab[4]);
}

// construct_point: dim 5

template <typename Kernel>
typename Kernel::Point_d construct_point(const Kernel &k,
                                         typename Kernel::FT x1, typename Kernel::FT x2, typename Kernel::FT x3,
                                         typename Kernel::FT x4, typename Kernel::FT x5) {
  typename Kernel::FT tab[5];
  tab[0] = x1;
  tab[1] = x2;
  tab[2] = x3;
  tab[3] = x4;
  tab[4] = x5;
  return k.construct_point_d_object()(5, &tab[0], &tab[5]);
}

// construct_point: dim 6

template <typename Kernel>
typename Kernel::Point_d construct_point(const Kernel &k,
                                         typename Kernel::FT x1, typename Kernel::FT x2, typename Kernel::FT x3,
                                         typename Kernel::FT x4, typename Kernel::FT x5, typename Kernel::FT x6) {
  typename Kernel::FT tab[6];
  tab[0] = x1;
  tab[1] = x2;
  tab[2] = x3;
  tab[3] = x4;
  tab[4] = x5;
  tab[5] = x6;
  return k.construct_point_d_object()(6, &tab[0], &tab[6]);
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_plane(std::size_t num_points, int intrinsic_dim,
                                                               int ambient_dim,
                                                               double coord_min = -5., double coord_max = 5.) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    std::vector<FT> pt(ambient_dim, FT(0));
    for (int j = 0; j < intrinsic_dim; ++j)
      pt[j] = rng.get_double(coord_min, coord_max);

    Point p = k.construct_point_d_object()(ambient_dim, pt.begin(), pt.end());
    points.push_back(p);
    ++i;
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_moment_curve(std::size_t num_points, int dim,
                                                                      typename Kernel::FT min_x,
                                                                      typename Kernel::FT max_x) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    FT x = rng.get_double(min_x, max_x);
    std::vector<FT> coords;
    coords.reserve(dim);
    for (int p = 1; p <= dim; ++p)
      coords.push_back(std::pow(CGAL::to_double(x), p));
    Point p = k.construct_point_d_object()(
                                           dim, coords.begin(), coords.end());
    points.push_back(p);
    ++i;
  }
  return points;
}


// R = big radius, r = small radius
template <typename Kernel/*, typename TC_basis*/>
std::vector<typename Kernel::Point_d> generate_points_on_torus_3D(std::size_t num_points, double R, double r,
                                                                  bool uniform = false) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    FT u, v;
    if (uniform) {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    } else {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    Point p = construct_point(k,
                              (R + r * std::cos(u)) * std::cos(v),
                              (R + r * std::cos(u)) * std::sin(v),
                              r * std::sin(u));
    points.push_back(p);
    ++i;
  }
  return points;
}

// "Private" function used by generate_points_on_torus_d
template <typename Kernel, typename OutputIterator>
static void generate_uniform_points_on_torus_d(const Kernel &k, int dim, std::size_t num_slices,
                                               OutputIterator out,
                                               double radius_noise_percentage = 0.,
                                               std::vector<typename Kernel::FT> current_point =
                                                       std::vector<typename Kernel::FT>()) {
  CGAL::Random rng;
  int point_size = static_cast<int>(current_point.size());
  if (point_size == 2 * dim) {
    *out++ = k.construct_point_d_object()(point_size, current_point.begin(), current_point.end());
  } else {
    for (std::size_t slice_idx = 0; slice_idx < num_slices; ++slice_idx) {
      double radius_noise_ratio = 1.;
      if (radius_noise_percentage > 0.) {
        radius_noise_ratio = rng.get_double(
                                            (100. - radius_noise_percentage) / 100.,
                                            (100. + radius_noise_percentage) / 100.);
      }
      std::vector<typename Kernel::FT> cp2 = current_point;
      double alpha = 6.2832 * slice_idx / num_slices;
      cp2.push_back(radius_noise_ratio * std::cos(alpha));
      cp2.push_back(radius_noise_ratio * std::sin(alpha));
      generate_uniform_points_on_torus_d(
                                         k, dim, num_slices, out, radius_noise_percentage, cp2);
    }
  }
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_torus_d(std::size_t num_points, int dim, bool uniform = false,
                                                                 double radius_noise_percentage = 0.) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;

  std::vector<Point> points;
  points.reserve(num_points);
  if (uniform) {
    std::size_t num_slices = (std::size_t)std::pow(num_points, 1. / dim);
    generate_uniform_points_on_torus_d(
                                       k, dim, num_slices, std::back_inserter(points), radius_noise_percentage);
  } else {
    for (std::size_t i = 0; i < num_points;) {
      double radius_noise_ratio = 1.;
      if (radius_noise_percentage > 0.) {
        radius_noise_ratio = rng.get_double(
                                            (100. - radius_noise_percentage) / 100.,
                                            (100. + radius_noise_percentage) / 100.);
      }
      std::vector<typename Kernel::FT> pt;
      pt.reserve(dim * 2);
      for (int curdim = 0; curdim < dim; ++curdim) {
        FT alpha = rng.get_double(0, 6.2832);
        pt.push_back(radius_noise_ratio * std::cos(alpha));
        pt.push_back(radius_noise_ratio * std::sin(alpha));
      }

      Point p = k.construct_point_d_object()(pt.begin(), pt.end());
      points.push_back(p);
      ++i;
    }
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_sphere_d(std::size_t num_points, int dim, double radius,
                                                                  double radius_noise_percentage = 0.) {
  typedef typename Kernel::Point_d Point;
  Kernel k;
  CGAL::Random rng;
  CGAL::Random_points_on_sphere_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    Point p = *generator++;
    if (radius_noise_percentage > 0.) {
      double radius_noise_ratio = rng.get_double(
                                                 (100. - radius_noise_percentage) / 100.,
                                                 (100. + radius_noise_percentage) / 100.);

      typename Kernel::Point_to_vector_d k_pt_to_vec =
          k.point_to_vector_d_object();
      typename Kernel::Vector_to_point_d k_vec_to_pt =
          k.vector_to_point_d_object();
      typename Kernel::Scaled_vector_d k_scaled_vec =
          k.scaled_vector_d_object();
      p = k_vec_to_pt(k_scaled_vec(k_pt_to_vec(p), radius_noise_ratio));
    }
    points.push_back(p);
    ++i;
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_in_ball_d(std::size_t num_points, int dim, double radius) {
  typedef typename Kernel::Point_d Point;
  Kernel k;
  CGAL::Random rng;
  CGAL::Random_points_in_ball_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    Point p = *generator++;
    points.push_back(p);
    ++i;
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_in_cube_d(std::size_t num_points, int dim, double radius) {
  typedef typename Kernel::Point_d Point;
  Kernel k;
  CGAL::Random rng;
  CGAL::Random_points_in_cube_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    Point p = *generator++;
    points.push_back(p);
    ++i;
  }
  return points;
}

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_two_spheres_d(std::size_t num_points, int dim, double radius,
                                                                       double distance_between_centers,
                                                                       double radius_noise_percentage = 0.) {
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::Vector_d Vector;
  Kernel k;
  CGAL::Random rng;
  CGAL::Random_points_on_sphere_d<Point> generator(dim, radius);
  std::vector<Point> points;
  points.reserve(num_points);

  std::vector<FT> t(dim, FT(0));
  t[0] = distance_between_centers;
  Vector c1_to_c2(t.begin(), t.end());

  for (std::size_t i = 0; i < num_points;) {
    Point p = *generator++;
    if (radius_noise_percentage > 0.) {
      double radius_noise_ratio = rng.get_double(
                                                 (100. - radius_noise_percentage) / 100.,
                                                 (100. + radius_noise_percentage) / 100.);

      typename Kernel::Point_to_vector_d k_pt_to_vec =
          k.point_to_vector_d_object();
      typename Kernel::Vector_to_point_d k_vec_to_pt =
          k.vector_to_point_d_object();
      typename Kernel::Scaled_vector_d k_scaled_vec =
          k.scaled_vector_d_object();
      p = k_vec_to_pt(k_scaled_vec(k_pt_to_vec(p), radius_noise_ratio));
    }

    typename Kernel::Translated_point_d k_transl =
        k.translated_point_d_object();
    Point p2 = k_transl(p, c1_to_c2);
    points.push_back(p);
    points.push_back(p2);
    i += 2;
  }
  return points;
}

// Product of a 3-sphere and a circle => d = 3 / D = 5

template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_3sphere_and_circle(std::size_t num_points,
                                                                            double sphere_radius) {
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_d Point;
  Kernel k;
  CGAL::Random rng;
  CGAL::Random_points_on_sphere_d<Point> generator(3, sphere_radius);
  std::vector<Point> points;
  points.reserve(num_points);

  typename Kernel::Compute_coordinate_d k_coord =
      k.compute_coordinate_d_object();
  for (std::size_t i = 0; i < num_points;) {
    Point p_sphere = *generator++;  // First 3 coords

    FT alpha = rng.get_double(0, 6.2832);
    std::vector<FT> pt(5);
    pt[0] = k_coord(p_sphere, 0);
    pt[1] = k_coord(p_sphere, 1);
    pt[2] = k_coord(p_sphere, 2);
    pt[3] = std::cos(alpha);
    pt[4] = std::sin(alpha);
    Point p(pt.begin(), pt.end());
    points.push_back(p);
    ++i;
  }
  return points;
}

// a = big radius, b = small radius
template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_klein_bottle_3D(std::size_t num_points, double a, double b,
                                                                         bool uniform = false) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    FT u, v;
    if (uniform) {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    } else {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    double tmp = cos(u / 2) * sin(v) - sin(u / 2) * sin(2. * v);
    Point p = construct_point(k,
                              (a + b * tmp) * cos(u),
                              (a + b * tmp) * sin(u),
                              b * (sin(u / 2) * sin(v) + cos(u / 2) * sin(2. * v)));
    points.push_back(p);
    ++i;
  }
  return points;
}

// a = big radius, b = small radius
template <typename Kernel>
std::vector<typename Kernel::Point_d> generate_points_on_klein_bottle_4D(std::size_t num_points, double a, double b,
                                                                         double noise = 0., bool uniform = false) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    FT u, v;
    if (uniform) {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    } else {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    Point p = construct_point(k,
                              (a + b * cos(v)) * cos(u) + (noise == 0. ? 0. : rng.get_double(0, noise)),
                              (a + b * cos(v)) * sin(u) + (noise == 0. ? 0. : rng.get_double(0, noise)),
                              b * sin(v) * cos(u / 2) + (noise == 0. ? 0. : rng.get_double(0, noise)),
                              b * sin(v) * sin(u / 2) + (noise == 0. ? 0. : rng.get_double(0, noise)));
    points.push_back(p);
    ++i;
  }
  return points;
}


// a = big radius, b = small radius

template <typename Kernel>
std::vector<typename Kernel::Point_d>
generate_points_on_klein_bottle_variant_5D(
                                           std::size_t num_points, double a, double b, bool uniform = false) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::FT FT;
  Kernel k;
  CGAL::Random rng;

  // if uniform
  std::size_t num_lines = (std::size_t)sqrt(num_points);

  std::vector<Point> points;
  points.reserve(num_points);
  for (std::size_t i = 0; i < num_points;) {
    FT u, v;
    if (uniform) {
      std::size_t k1 = i / num_lines;
      std::size_t k2 = i % num_lines;
      u = 6.2832 * k1 / num_lines;
      v = 6.2832 * k2 / num_lines;
    } else {
      u = rng.get_double(0, 6.2832);
      v = rng.get_double(0, 6.2832);
    }
    FT x1 = (a + b * cos(v)) * cos(u);
    FT x2 = (a + b * cos(v)) * sin(u);
    FT x3 = b * sin(v) * cos(u / 2);
    FT x4 = b * sin(v) * sin(u / 2);
    FT x5 = x1 + x2 + x3 + x4;

    Point p = construct_point(k, x1, x2, x3, x4, x5);
    points.push_back(p);
    ++i;
  }
  return points;
}

}  // namespace Gudhi

#endif  // RANDOM_POINT_GENERATORS_H_
