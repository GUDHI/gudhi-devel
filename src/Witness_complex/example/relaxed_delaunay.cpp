/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2016  INRIA Sophia Antipolis-Méditerranée (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gudhi/Simplex_tree.h>
#include <gudhi/Relaxed_witness_complex.h>
#include <gudhi/Dim_lists.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Persistent_cohomology.h>
#include "Landmark_choice_random_knn.h"
#include "Landmark_choice_sparsification.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <algorithm>
#include <set>
#include <queue>
#include <iterator>
#include <string>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
//#include <CGAL/Sphere_d.h>

#include "generators.h"
#include "output.h"

using namespace Gudhi;
using namespace Gudhi::witness_complex;
using namespace Gudhi::persistent_cohomology;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef K::Sphere_d Sphere_d;
typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;

typedef std::vector<Point_d> Point_Vector;
typedef Relaxed_witness_complex< Simplex_tree<> > RelaxedWitnessComplex;
typedef Simplex_tree<>::Simplex_handle Simplex_handle;



/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust(std::string file_name, std::vector< std::vector< double > > & points) {
  std::ifstream in_file(file_name.c_str(), std::ios::in);
  if (!in_file.is_open()) {
    std::cerr << "Unable to open file " << file_name << std::endl;
    return;
  }
  std::string line;
  double x;
  while (getline(in_file, line)) {
    std::vector< double > point;
    std::istringstream iss(line);
    while (iss >> x) {
      point.push_back(x);
    }
    if (point.size() != 1)
      points.push_back(point);
  }
  in_file.close();
}

int main (int argc, char * const argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
              << " 1 file_name alpha limD\n";
    return 0;
  }
  std::string file_name = argv[1];
  double alpha2  = atof(argv[2]);
  int limD = atoi(argv[3]);
  
  // Read points
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  generate_points_random_box(point_vector, 200, 2);
  write_points(file_name, point_vector);
  
  std::cout << "The file contains " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  // 1. Compute Delaunay centers
  Delaunay_triangulation delaunay(point_vector[0].size());
  delaunay.insert(point_vector.begin(), point_vector.end());
  Point_Vector del_centers;
  for (auto f_it = delaunay.full_cells_begin(); f_it != delaunay.full_cells_end(); ++f_it) {
    if (delaunay.is_infinite(f_it))
      continue;
    Point_Vector vertices;
    for (auto v_it = f_it->vertices_begin(); v_it != f_it->vertices_end(); ++v_it)
      vertices.push_back((*v_it)->point());
    Sphere_d sphere(vertices.begin(), vertices.end());
    del_centers.push_back(sphere.center());
  }
  std::cout << "Delaunay center count: " << del_centers.size() << ".\n";

  // 2. Build Relaxed Witness Complex
  std::vector<std::vector<int>> knn;
  std::vector<std::vector<double>> distances;
  Gudhi::witness_complex::build_distance_matrix(del_centers,   // aka witnesses
                                                point_vector,  // aka landmarks
                                                alpha2,
                                                limD,
                                                knn,
                                                distances);
  
  write_wl("wl_distances.txt", distances);
  Simplex_tree<> simplex_tree;
  Gudhi::witness_complex::Relaxed_witness_complex<Simplex_tree<>> rwc(distances,
                                                                      knn,
                                                                      simplex_tree,
                                                                      point_vector.size(),
                                                                      alpha2,
                                                                      limD);
  std::vector<int> dim_simplices(limD+1);
  for (auto sh: simplex_tree.complex_simplex_range()) {
    dim_simplices[simplex_tree.dimension(sh)]++;
  }
  for (unsigned i =0; i != dim_simplices.size(); ++i)
    std::cout << "dim[" << i << "]: " << dim_simplices[i] << " simplices.\n";

  std::vector<int> landmarks_ind;
  for (unsigned i = 0; i < point_vector.size(); ++i)
    landmarks_ind.push_back(i);
  write_witness_mesh(point_vector, landmarks_ind, simplex_tree, simplex_tree.complex_simplex_range(), true, true, "relaxed_delaunay.mesh");
  
  // 3. Check if the thing is Relaxed Delaunay
  for (auto sh: simplex_tree.complex_simplex_range()) {
    Point_Vector vertices;
    for (auto v: simplex_tree.simplex_vertex_range(sh))
      vertices.push_back(point_vector[v]);
    Sphere_d sphere(vertices.begin(), vertices.end());
    Point_d center = sphere.center();
    double r2 = sphere.squared_radius();
    typename K::Squared_distance_d dist2;
    std::vector<int> v_inds;
    for (auto v: simplex_tree.simplex_vertex_range(sh))
      v_inds.push_back(v);
    auto range_begin = std::begin(v_inds);
    auto range_end = std::end(v_inds);
    if (simplex_tree.dimension(sh) == (int)point_vector[0].size())
    for (auto v: simplex_tree.complex_vertex_range())
      if (std::find(range_begin, range_end, v) == range_end) {
        if (dist2(point_vector[v], center) < r2 - alpha2)
          std::cout << "WARNING! The vertex " << point_vector[v] << " is inside the (r2-alpha2)-ball (" << center << ", " << r2 << ") distance is " << dist2(point_vector[v], center) << "\n";
      }
  }
}
