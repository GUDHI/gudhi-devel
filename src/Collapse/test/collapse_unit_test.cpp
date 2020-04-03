/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>  // std::pair, std::make_pair
#include <cmath>  // float comparison
#include <limits>
#include <functional>  // greater
#include <tuple>  // std::tie

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "collapse"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

//  ^
// /!\ Nothing else from Simplex_tree shall be included to test includes are well defined.
#include "gudhi/FlagComplexSpMatrix.h"
#include "gudhi/Rips_edge_list.h"

using namespace Gudhi;

// Types definition
using Vector_of_points = std::vector<std::vector<double>>;

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = double;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Rips_edge_list = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;


BOOST_AUTO_TEST_CASE(collapse) {
  typedef size_t Vertex_handle;
  typedef std::vector<std::tuple<Filtration_value, Vertex_handle, Vertex_handle>> Filtered_sorted_edge_list;

  std::size_t number_of_points;
  std::string off_file_points;
  std::string filediag;
  int dim_max;
  int p;
  double min_persistence;

  Map map_empty;

  Distance_matrix sparse_distances;


  Vector_of_points point_vector {{0., 0.},{0., 1.},{1., 0.},{1., 1.}};

  int dimension = point_vector[0].dimension();
  number_of_points = point_vector.size();
  std::cout << "Successfully read " << number_of_points << " point_vector.\n";
  std::cout << "Ambient dimension is " << dimension << ".\n";

  std::cout << "Point Set Generated." << std::endl;

  double threshold = 1.;
  Filtered_sorted_edge_list edge_t;
  std::cout << "Computing the one-skeleton for threshold: " << threshold << std::endl;

  Rips_edge_list Rips_edge_list_from_file(point_vector, threshold, Gudhi::Euclidean_distance());
  Rips_edge_list_from_file.create_edges(edge_t);

  std::cout << "Sorted edge list computed" << std::endl;
  std::cout << "Total number of edges before collapse are: " << edge_t.size() << std::endl;

  if (edge_t.size() <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  // Now we will perform filtered edge collapse to sparsify the edge list edge_t.
  std::cout << "Filtered edge collapse begins" << std::endl;
  FlagComplexSpMatrix mat_filt_edge_coll(number_of_points, edge_t);
  std::cout << "Matrix instansiated" << std::endl;
  Filtered_sorted_edge_list collapse_edges;
  collapse_edges = mat_filt_edge_coll.filtered_edge_collapse();
 
}


