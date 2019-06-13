/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef MODEL_MODEL_H_
#define MODEL_MODEL_H_

#include <gudhi/Clock.h>
#include <gudhi/Skeleton_blocker/Skeleton_blocker_simple_geometric_traits.h>
#include <gudhi/Skeleton_blocker_geometric_complex.h>
#include <gudhi/Off_reader.h>

#include <CGAL/Euclidean_distance.h>

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "utils/UI_utils.h"
#include "utils/Lloyd_builder.h"
#include "utils/Rips_builder.h"
#include "utils/K_nearest_builder.h"
#include "utils/Vertex_collapsor.h"
#include "utils/Edge_collapsor.h"
#include "utils/Edge_contractor.h"
#include "utils/Persistence_compute.h"
#include "utils/Critical_points.h"
#include "utils/Is_manifold.h"

#include "Complex_typedefs.h"

template<typename Complex>
class CGAL_geometric_flag_complex_wrapper {
  Complex& complex_;
  typedef typename Complex::Vertex_handle Vertex_handle;
  typedef typename Complex::Point Point;

  const bool load_only_points_;

 public:
  CGAL_geometric_flag_complex_wrapper(Complex& complex, bool load_only_points = false) :
      complex_(complex),
      load_only_points_(load_only_points) { }

  void init(int dim, int num_vertices, int num_max_faces, int num_edges) const { }

  void point(const std::vector<double>& coords) {
    Point p(coords.size(), coords.begin(), coords.end());
    complex_.add_vertex(p);
  }

  void maximal_face(std::vector<int> vertices) {
    if (!load_only_points_) {
      // std::cout << "size:" << vertices.size() << std::endl;
      for (std::size_t i = 0; i < vertices.size(); ++i)
        for (std::size_t j = i + 1; j < vertices.size(); ++j)
          complex_.add_edge_without_blockers(Vertex_handle(vertices[i]), Vertex_handle(vertices[j]));
    }
  }

  void done() const { }
};

class Model {
 public:
  Complex complex_;
  typedef Complex::Vertex_handle Vertex_handle;

  Model() : complex_() { }

 public:
  void off_file_open(const std::string& name_file) {
    UIDBGMSG("load off file", name_file);
    complex_.clear();
    CGAL_geometric_flag_complex_wrapper<Complex> read_wraper(complex_, false);
    Gudhi::read_off(name_file, read_wraper);
  }

  void off_points_open(const std::string& name_file) {
    UIDBGMSG("load off points", name_file);
    complex_.clear();
    CGAL_geometric_flag_complex_wrapper<Complex> read_wraper(complex_, true);
    Gudhi::read_off(name_file, read_wraper);
  }

  void off_file_save(const std::string& name_file) {
    UIDBG("save off file");
    UIDBG("save off off_points_save");
    std::ofstream file(name_file);
    if (file.is_open()) {
      file << "OFF\n";
      file << complex_.num_vertices() << " " << complex_.num_edges() << " 0\n";
      for (auto v : complex_.vertex_range()) {
        const auto& pt(complex_.point(v));
        for (auto it = pt.cartesian_begin(); it != pt.cartesian_end(); ++it)
          file << *it << " ";
        file << std::endl;
      }
      for (auto e : complex_.edge_range())
        file << "2 " << complex_.first_vertex(e) << " " << complex_.second_vertex(e) << "\n";
      file.close();
    } else {
      std::cerr << "Could not open file " << name_file << std::endl;
    }
  }

  void off_points_save(const std::string& name_file) {
    UIDBG("save off off_points_save");
    std::ofstream file(name_file);
    if (file.is_open()) {
      file << "OFF\n";
      file << complex_.num_vertices() << " 0 0\n";
      for (auto v : complex_.vertex_range()) {
        const auto& pt(complex_.point(v));
        for (auto it = pt.cartesian_begin(); it != pt.cartesian_end(); ++it)
          file << *it << " ";
        file << std::endl;
      }
      file.close();
    } else {
      std::cerr << "Could not open file " << name_file << std::endl;
    }
  }

  // point sets operations
  void uniform_noise(double amplitude) {
    UIDBG("unif noise");
    for (auto v : complex_.vertex_range())
      complex_.point(v) = add_uniform_noise(complex_.point(v), amplitude);
  }

 private:
  Point add_uniform_noise(const Point& point, double amplitude) {
    std::vector<double> new_point(point.dimension());
    for (int i = 0; i < point.dimension(); ++i) {
      new_point[i] = point[i] + (rand() % 2 - .5) * amplitude;
    }
    return Point(point.dimension(), new_point.begin(), new_point.end());
  }

 public:
  void lloyd(int num_iterations, int num_closest_neighbors) {
    UIDBG("lloyd");
    Lloyd_builder<Complex> lloyd_builder(complex_, 1);
  }

  double squared_eucl_distance(const Point& p1, const Point& p2) const {
    return Geometry_trait::Squared_distance_d()(p1, p2);
  }

  // complex operations from points

  void build_rips(double alpha) {
    UIDBG("build_rips");
    Rips_builder<Complex> rips_builder(complex_, alpha);
  }

  void build_k_nearest_neighbors(unsigned k) {
    UIDBG("build_k_nearest");
    complex_.keep_only_vertices();
    K_nearest_builder<Complex> k_nearest_builder(complex_, k);
  }

  void build_delaunay() {
    UIDBG("build_delaunay");
    complex_.keep_only_vertices();
  }

  void contract_edges(unsigned num_contractions) {
    Gudhi::Clock c;
    Edge_contractor<Complex> contractor(complex_, num_contractions);
    std::cout << "Time to simplify: " << c.num_seconds() << "s" << std::endl;
  }

  void collapse_vertices(unsigned num_collapses) {
    auto old_num_vertices = complex_.num_vertices();
    Vertex_collapsor<Complex> collapsor(complex_, complex_.num_vertices());
    UIDBGMSG("num vertices collapsed:", old_num_vertices - complex_.num_vertices());
  }

  void collapse_edges(unsigned num_collapses) {
    Edge_collapsor<Complex> collapsor(complex_, num_collapses);
  }

  void show_graph_stats() {
    std::cout << "++++++ Graph stats +++++++" << std::endl;
    std::cout << "Num vertices : " << complex_.num_vertices() << std::endl;
    std::cout << "Num edges : " << complex_.num_edges() << std::endl;
    std::cout << "Num connected components : " << complex_.num_connected_components() << std::endl;
    std::cout << "Min/avg/max degree : " << min_degree() << "/" << avg_degree() << "/" << max_degree() << std::endl;
    std::cout << "Num connected components : " << complex_.num_connected_components() << std::endl;
    std::cout << "Num connected components : " << complex_.num_connected_components() << std::endl;
    std::cout << "+++++++++++++++++++++++++" << std::endl;
  }

 private:
  int min_degree() const {
    int res = (std::numeric_limits<int>::max)();
    for (auto v : complex_.vertex_range())
      res = (std::min)(res, complex_.degree(v));
    return res;
  }

  int max_degree() const {
    int res = 0;
    for (auto v : complex_.vertex_range())
      res = (std::max)(res, complex_.degree(v));
    return res;
  }

  int avg_degree() const {
    int res = 0;
    for (auto v : complex_.vertex_range())
      res += complex_.degree(v);
    return res / complex_.num_vertices();
  }

 public:
  void show_complex_stats() {
    std::cout << "++++++ Mesh stats +++++++" << std::endl;
    std::cout << "Num vertices : " << complex_.num_vertices() << std::endl;
    std::cout << "Num edges : " << complex_.num_edges() << std::endl;
    std::cout << "Num connected components : " << complex_.num_connected_components() << std::endl;
    std::cout << "+++++++++++++++++++++++++" << std::endl;
  }

  void show_complex_dimension() {
    unsigned num_simplices = 0;
    int euler = 0;
    int dimension = 0;
    Gudhi::Clock clock;
    for (const auto &s : complex_.complex_simplex_range()) {
      num_simplices++;
      dimension = (std::max)(s.dimension(), dimension);
      if (s.dimension() % 2 == 0)
        euler += 1;
      else
        euler -= 1;
    }
    clock.end();
    std::cout << "++++++ Mesh dimension +++++++" << std::endl;
    std::cout << "Dimension : " << dimension << std::endl;
    std::cout << "Euler characteristic : " << euler << std::endl;
    std::cout << "Num simplices : " << num_simplices << std::endl;
    std::cout << "Total time: " << clock << std::endl;
    std::cout << "Time per simplex: " << clock.num_seconds() / num_simplices << " s" << std::endl;
    std::cout << "+++++++++++++++++++++++++" << std::endl;
  }

  void show_homology_group() {
#ifdef _WIN32
    std::cout << "Works only on linux x64 for the moment\n";
#else
    Gudhi::Clock clock;
    run_chomp();
    clock.end();
#endif
  }

  void show_euler_characteristic() {
    unsigned num_simplices = 0;
    int euler = 0;
    int dimension = 0;
    for (const auto &s : complex_.complex_simplex_range()) {
      num_simplices++;
      dimension = (std::max)(s.dimension(), dimension);
      if (s.dimension() % 2 == 0)
        euler += 1;
      else
        euler -= 1;
    }
    std::cout << "Saw " << num_simplices << " simplices with maximum dimension " << dimension << std::endl;
    std::cout << "The euler characteristic is : " << euler << std::endl;
  }

  void show_persistence(int p, double threshold, int max_dim, double min_pers) {
    Persistence_compute<Complex> persistence(complex_, std::cout, Persistence_params(p, threshold, max_dim, min_pers));
  }

  void show_critical_points(double max_distance) {
    Critical_points<Complex> critical_points(complex_, std::cout, max_distance);
  }

  void show_is_manifold() {
    unsigned dim;
    bool is_manifold;
    Is_manifold<Complex> test_manifold(complex_, dim, is_manifold);

    if (is_manifold) {
      std::cout << "The complex is a " << dim << "-manifold\n";
    } else {
      if (dim < 4) {
        std::cout << "The complex has dimension greater than " << dim << " and is not a manifold\n";
      } else {
        std::cout << "The complex has dimension>=4 and may or may not be a manifold\n";
      }
    }
  }

 private:
  void run_chomp() {
    save_complex_in_file_for_chomp();
    std::cout << "Call CHOMP library\n";
    int returnValue = system("utils/homsimpl chomp.sim");
    std::cout << "CHOMP returns" << returnValue << std::endl;
  }

  void save_complex_in_file_for_chomp() {
    std::ofstream file;
    file.open("chomp.sim");
    for (const auto &s : complex_.complex_simplex_range()) {
      bool first = true;
      file << "(";
      for (auto x : s) {
        if (first)
          first = false;
        else
          file << ",";
        file << x;
      }
      file << ")\n";
    }
  }

 public:
  unsigned num_vertices() const {
    return complex_.num_vertices();
  }

  unsigned num_edges() const {
    return complex_.num_edges();
  }
};

#endif  // MODEL_MODEL_H_
