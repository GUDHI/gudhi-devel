/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */
#include <gudhi/Edge_contraction.h>
#include <gudhi/Skeleton_blocker.h>
#include <gudhi/Off_reader.h>
#include <gudhi/Point.h>
#include <gudhi/Clock.h>

#include <iostream>

struct Geometry_trait {
  typedef Point_d Point;
};

using Complex_geometric_traits = Gudhi::skeleton_blocker::Skeleton_blocker_simple_geometric_traits<Geometry_trait>;
using Complex = Gudhi::skeleton_blocker::Skeleton_blocker_geometric_complex< Complex_geometric_traits >;
using Profile = Gudhi::contraction::Edge_profile<Complex>;
using Complex_contractor = Gudhi::contraction::Skeleton_blocker_contractor<Complex>;


template<typename ComplexType>
void build_rips(ComplexType& complex, double offset) {
  if (offset <= 0) return;
  auto vertices = complex.vertex_range();
  for (auto p = vertices.begin(); p != vertices.end(); ++p)
    for (auto q = p; ++q != vertices.end(); /**/) {
      if (squared_dist(complex.point(*p), complex.point(*q)) < 4 * offset * offset)
        complex.add_edge_without_blockers(*p, *q);
    }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage " << argv[0] << " ../../../data/meshes/SO3_10000.off 0.3 to load the file " <<
        "../../data/SO3_10000.off and contract the Rips complex built with paremeter 0.3.\n";
    return -1;
  }

  Complex complex;

  // load only the points
  Gudhi::skeleton_blocker::Skeleton_blocker_off_reader<Complex> off_reader(argv[1], complex, true);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file:" << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Build the Rips complex with " << complex.num_vertices() << " vertices" << std::endl;

  build_rips(complex, atof(argv[2]));

  Gudhi::Clock contraction_chrono("Time to simplify and enumerate simplices");

  std::cout << "Initial complex has " <<
      complex.num_vertices() << " vertices and " <<
      complex.num_edges() << " edges" << std::endl;

  Complex_contractor contractor(complex,
                                new Gudhi::contraction::Edge_length_cost<Profile>,
                                Gudhi::contraction::make_first_vertex_placement<Profile>(),
                                Gudhi::contraction::make_link_valid_contraction<Profile>(),
                                Gudhi::contraction::make_remove_popable_blockers_visitor<Profile>());
  contractor.contract_edges();

  std::cout << "Counting final number of simplices \n";
  unsigned num_simplices = std::distance(complex.complex_simplex_range().begin(), complex.complex_simplex_range().end());

  std::cout << "Final complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_edges() << " edges, " <<
      complex.num_blockers() << " blockers and " <<
      num_simplices << " simplices" << std::endl;

  std::cout << contraction_chrono;

  return EXIT_SUCCESS;
}


