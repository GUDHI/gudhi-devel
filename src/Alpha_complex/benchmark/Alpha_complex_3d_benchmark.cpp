#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Alpha_complex.h>
// to construct a simplex_tree from alpha complex
#include <gudhi/Simplex_tree.h>
#include <gudhi/random_point_generators.h>
#include <gudhi/Clock.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits
#include <fstream>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/Random.h>

std::ofstream results_csv("results.csv");

template <typename Kernel>
void benchmark_points_on_torus_dD(const std::string& msg) {
  std::cout << "+ " << msg << std::endl;

  results_csv << "\"" << msg << "\";" << std::endl;
  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"alpha_creation_time(sec.)\";"
              << "\"simplex_creation_time(sec.)\";" << std::endl;

  using K = CGAL::Epick_d<CGAL::Dimension_tag<3>>;
  for (int nb_points = 1000; nb_points <= 125000; nb_points *= 5) {
    std::cout << "  Alpha complex dD on torus with " << nb_points << " points." << std::endl;
    std::vector<K::Point_d> points_on_torus = Gudhi::generate_points_on_torus_3D<K>(nb_points, 1.0, 0.5);
    std::vector<typename Kernel::Point_d> points;

    for (auto p : points_on_torus) {
      points.push_back(typename Kernel::Point_d(p.begin(), p.end()));
    }

    Gudhi::Clock ac_create_clock("    benchmark_points_on_torus_dD - Alpha complex 3d creation");
    ac_create_clock.begin();
    Gudhi::alpha_complex::Alpha_complex<Kernel> alpha_complex_from_points(points);
    ac_create_clock.end();
    std::cout << ac_create_clock;

    Gudhi::Simplex_tree<> simplex;
    Gudhi::Clock st_create_clock("    benchmark_points_on_torus_dD - simplex tree creation");
    st_create_clock.begin();
    alpha_complex_from_points.create_complex(simplex);
    st_create_clock.end();
    std::cout << st_create_clock;

    results_csv << nb_points << ";" << simplex.num_simplices() << ";" << ac_create_clock.num_seconds() << ";"
                << st_create_clock.num_seconds() << ";" << std::endl;

    std::cout << "    benchmark_points_on_torus_dD - nb simplices = " << simplex.num_simplices() << std::endl;
  }
}

template <typename Alpha_complex_3d>
void benchmark_points_on_torus_3D(const std::string& msg) {
  using K = CGAL::Epick_d<CGAL::Dimension_tag<3>>;
  std::cout << "+ " << msg << std::endl;

  results_csv << "\"" << msg << "\";" << std::endl;
  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"alpha_creation_time(sec.)\";"
              << "\"simplex_creation_time(sec.)\";" << std::endl;

  for (int nb_points = 1000; nb_points <= 125000; nb_points *= 5) {
    std::cout << "  Alpha complex 3d on torus with " << nb_points << " points." << std::endl;
    std::vector<K::Point_d> points_on_torus = Gudhi::generate_points_on_torus_3D<K>(nb_points, 1.0, 0.5);
    std::vector<typename Alpha_complex_3d::Point_3> points;

    for (auto p : points_on_torus) {
      points.push_back(typename Alpha_complex_3d::Point_3(p[0], p[1], p[2]));
    }

    Gudhi::Clock ac_create_clock("    benchmark_points_on_torus_3D - Alpha complex 3d creation");
    ac_create_clock.begin();
    Alpha_complex_3d alpha_complex_from_points(points);
    ac_create_clock.end();
    std::cout << ac_create_clock;

    Gudhi::Simplex_tree<> simplex;
    Gudhi::Clock st_create_clock("    benchmark_points_on_torus_3D - simplex tree creation");
    st_create_clock.begin();
    alpha_complex_from_points.create_complex(simplex);
    st_create_clock.end();
    std::cout << st_create_clock;

    results_csv << nb_points << ";" << simplex.num_simplices() << ";" << ac_create_clock.num_seconds() << ";"
                << st_create_clock.num_seconds() << ";" << std::endl;

    std::cout << "    benchmark_points_on_torus_3D - nb simplices = " << simplex.num_simplices() << std::endl;
  }
}

template <typename Weighted_alpha_complex_3d>
void benchmark_weighted_points_on_torus_3D(const std::string& msg) {
  using K = CGAL::Epick_d<CGAL::Dimension_tag<3>>;

  std::cout << "+ " << msg << std::endl;

  results_csv << "\"" << msg << "\";" << std::endl;
  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"alpha_creation_time(sec.)\";"
              << "\"simplex_creation_time(sec.)\";" << std::endl;

  CGAL::Random random(8);

  for (int nb_points = 1000; nb_points <= 125000; nb_points *= 5) {
    std::cout << "  Alpha complex 3d on torus with " << nb_points << " points." << std::endl;
    std::vector<K::Point_d> points_on_torus = Gudhi::generate_points_on_torus_3D<K>(nb_points, 1.0, 0.5);

    using Point = typename Weighted_alpha_complex_3d::Point_3;
    using Weighted_point = typename Weighted_alpha_complex_3d::Triangulation_3::Weighted_point;

    std::vector<Weighted_point> points;

    for (auto p : points_on_torus) {
      points.push_back(Weighted_point(Point(p[0], p[1], p[2]), 0.9 + random.get_double(0., 0.01)));
    }

    Gudhi::Clock ac_create_clock("    benchmark_weighted_points_on_torus_3D - Alpha complex 3d creation");
    ac_create_clock.begin();
    Weighted_alpha_complex_3d alpha_complex_from_points(points);
    ac_create_clock.end();
    std::cout << ac_create_clock;

    Gudhi::Simplex_tree<> simplex;
    Gudhi::Clock st_create_clock("    benchmark_weighted_points_on_torus_3D - simplex tree creation");
    st_create_clock.begin();
    alpha_complex_from_points.create_complex(simplex);
    st_create_clock.end();
    std::cout << st_create_clock;

    results_csv << nb_points << ";" << simplex.num_simplices() << ";" << ac_create_clock.num_seconds() << ";"
                << st_create_clock.num_seconds() << ";" << std::endl;

    std::cout << "    benchmark_weighted_points_on_torus_3D - nb simplices = " << simplex.num_simplices() << std::endl;
  }
}

template <typename Periodic_alpha_complex_3d>
void benchmark_periodic_points(const std::string& msg) {
  std::cout << "+ " << msg << std::endl;

  results_csv << "\"" << msg << "\";" << std::endl;
  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"alpha_creation_time(sec.)\";"
              << "\"simplex_creation_time(sec.)\";" << std::endl;

  CGAL::Random random(8);

  for (double nb_points = 10.; nb_points <= 40.; nb_points += 10.) {
    std::cout << "  Periodic alpha complex 3d with " << nb_points * nb_points * nb_points << " points." << std::endl;
    using Point = typename Periodic_alpha_complex_3d::Point_3;
    std::vector<Point> points;

    for (double i = 0; i < nb_points; i++) {
      for (double j = 0; j < nb_points; j++) {
        for (double k = 0; k < nb_points; k++) {
          points.push_back(
              Point(i + random.get_double(0., 0.1), j + random.get_double(0., 0.1), k + random.get_double(0., 0.1)));
        }
      }
    }

    Gudhi::Clock ac_create_clock("    benchmark_periodic_points - Alpha complex 3d creation");
    ac_create_clock.begin();
    Periodic_alpha_complex_3d alpha_complex_from_points(points, 0., 0., 0., nb_points, nb_points, nb_points);
    ac_create_clock.end();
    std::cout << ac_create_clock;

    Gudhi::Simplex_tree<> simplex;
    Gudhi::Clock st_create_clock("    benchmark_periodic_points - simplex tree creation");
    st_create_clock.begin();
    alpha_complex_from_points.create_complex(simplex);
    st_create_clock.end();
    std::cout << st_create_clock;

    results_csv << nb_points * nb_points * nb_points << ";" << simplex.num_simplices() << ";"
                << ac_create_clock.num_seconds() << ";" << st_create_clock.num_seconds() << ";" << std::endl;

    std::cout << "    benchmark_periodic_points - nb simplices = " << simplex.num_simplices() << std::endl;
  }
}

template <typename Weighted_periodic_alpha_complex_3d>
void benchmark_weighted_periodic_points(const std::string& msg) {
  std::cout << "+ " << msg << std::endl;

  results_csv << "\"" << msg << "\";" << std::endl;
  results_csv << "\"nb_points\";"
              << "\"nb_simplices\";"
              << "\"alpha_creation_time(sec.)\";"
              << "\"simplex_creation_time(sec.)\";" << std::endl;

  CGAL::Random random(8);

  for (double nb_points = 10.; nb_points <= 40.; nb_points += 10.) {
    std::cout << "  Weighted periodic alpha complex 3d with " << nb_points * nb_points * nb_points << " points."
              << std::endl;

    using Point = typename Weighted_periodic_alpha_complex_3d::Point_3;
    using Weighted_point = typename Weighted_periodic_alpha_complex_3d::Triangulation_3::Weighted_point;
    std::vector<Weighted_point> points;

    for (double i = 0; i < nb_points; i++) {
      for (double j = 0; j < nb_points; j++) {
        for (double k = 0; k < nb_points; k++) {
          points.push_back(Weighted_point(
              Point(i + random.get_double(0., 0.1), j + random.get_double(0., 0.1), k + random.get_double(0., 0.1)),
              random.get_double(0., (nb_points * nb_points) / 64.)));
        }
      }
    }

    Gudhi::Clock ac_create_clock("    benchmark_weighted_periodic_points - Alpha complex 3d creation");
    ac_create_clock.begin();
    Weighted_periodic_alpha_complex_3d alpha_complex_from_points(points, 0., 0., 0., nb_points, nb_points, nb_points);
    ac_create_clock.end();
    std::cout << ac_create_clock;

    Gudhi::Simplex_tree<> simplex;
    Gudhi::Clock st_create_clock("    benchmark_weighted_periodic_points - simplex tree creation");
    st_create_clock.begin();
    alpha_complex_from_points.create_complex(simplex);
    st_create_clock.end();
    std::cout << st_create_clock;

    results_csv << nb_points * nb_points * nb_points << ";" << simplex.num_simplices() << ";"
                << ac_create_clock.num_seconds() << ";" << st_create_clock.num_seconds() << ";" << std::endl;

    std::cout << "    benchmark_weighted_periodic_points - nb simplices = " << simplex.num_simplices() << std::endl;
  }
}

int main(int argc, char** argv) {
  benchmark_points_on_torus_dD<CGAL::Epick_d<CGAL::Dimension_tag<3>>>("Fast static dimension version");
  benchmark_points_on_torus_dD<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>("Fast dynamic dimension version");
  benchmark_points_on_torus_dD<CGAL::Epeck_d<CGAL::Dimension_tag<3>>>("Exact static dimension version");
  benchmark_points_on_torus_dD<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>("Exact dynamic dimension version");

  benchmark_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, false>>("Fast version");
  benchmark_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>>("Safe version");
  benchmark_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>>("Exact version");

  benchmark_weighted_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, false>>("Fast version");
  benchmark_weighted_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>>("Safe version");
  benchmark_weighted_points_on_torus_3D<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>>("Exact version");

  benchmark_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, true>>("Fast version");
  benchmark_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>>("Safe version");
  benchmark_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, true>>("Exact version");

  benchmark_weighted_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, true>>("Fast version");
  benchmark_weighted_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, true>>("Safe version");
  benchmark_weighted_periodic_points<
      Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, true>>("Exact version");

  return 0;
}
