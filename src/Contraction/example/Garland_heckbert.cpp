/* This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014  INRIA Sophia Antipolis-M�diterran�e (France)
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
 * 
 */


#ifndef GARLAND_HECKBERT_H_
#define GARLAND_HECKBERT_H_

#include <gudhi/Point.h>
#include <gudhi/Edge_contraction.h>
#include <gudhi/Skeleton_blocker.h>
#include <gudhi/Off_reader.h>

#include <boost/timer/timer.hpp>
#include <iostream>

#include "Garland_heckbert/Error_quadric.h"

using namespace std;
using namespace Gudhi;
using namespace skbl;
using namespace contraction;

struct Geometry_trait {
  typedef Point_d Point;
};

/**
 * The vertex stored in the complex contains a quadric.
 */
struct Garland_heckbert_traits : public Skeleton_blocker_simple_geometric_traits<Geometry_trait> {
 public:
  struct Garland_heckbert_vertex : public Simple_geometric_vertex {
    Error_quadric<Geometry_trait::Point> quadric;
  };
  typedef Garland_heckbert_vertex Graph_vertex;
};

typedef Skeleton_blocker_geometric_complex< Garland_heckbert_traits > Complex;
typedef Edge_profile<Complex> EdgeProfile;
typedef Skeleton_blocker_contractor<Complex> Complex_contractor;

/**
 * How the new vertex is placed after an edge collapse : here it is placed at
 * the point minimizing the cost of the quadric.
 */
class GH_placement : public Gudhi::contraction::Placement_policy<EdgeProfile> {
  Complex& complex_;

 public:
  typedef Gudhi::contraction::Placement_policy<EdgeProfile>::Placement_type Placement_type;

  GH_placement(Complex& complex) : complex_(complex) { }

  Placement_type operator()(const EdgeProfile& profile) const override {
    auto sum_quad(profile.v0().quadric);
    sum_quad += profile.v1().quadric;

    boost::optional<Point> min_quadric_pt(sum_quad.min_cost());
    if (min_quadric_pt)
      return Placement_type(*min_quadric_pt);
    else
      return profile.p0();
  }
};

/**
 * How much cost an edge collapse : here the costs is given by a quadric
 * which expresses a squared distances with triangles planes.
 */
class GH_cost : public Gudhi::contraction::Cost_policy<EdgeProfile> {
  Complex& complex_;

 public:
  typedef Gudhi::contraction::Cost_policy<EdgeProfile>::Cost_type Cost_type;

  GH_cost(Complex& complex) : complex_(complex) { }

  Cost_type operator()(EdgeProfile const& profile, boost::optional<Point> const& new_point) const override {
    Cost_type res;
    if (new_point) {
      auto sum_quad(profile.v0().quadric);
      sum_quad += profile.v1().quadric;
      res = sum_quad.cost(*new_point);
    }
    return res;
  }
};

/**
 * Visitor that is called at several moment.
 * Here we initializes the quadrics of every vertex at the on_started call back
 * and we update them when contracting an edge (the quadric become the sum of both quadrics).
 */
class GH_visitor : public Gudhi::contraction::Contraction_visitor<EdgeProfile> {
  Complex& complex_;

 public:
  GH_visitor(Complex& complex) : complex_(complex) { }

  // Compute quadrics for every vertex v
  // The quadric of v consists in the sum of quadric
  // of every triangles passing through v weighted by its area

  void on_started(Complex & complex) override {
    for (auto v : complex.vertex_range()) {
      auto & quadric_v(complex[v].quadric);
      for (auto t : complex.triangle_range(v)) {
        auto t_it = t.begin();
        const auto& p0(complex.point(*t_it++));
        const auto& p1(complex.point(*t_it++));
        const auto& p2(complex.point(*t_it++));
        quadric_v += Error_quadric<Point>(p0, p1, p2);
      }
    }
  }

  /**
   * @brief Called when an edge is about to be contracted and replaced by a vertex whose position is *placement.
   */
  void on_contracting(EdgeProfile const &profile, boost::optional< Point > placement)
  override {
    profile.v0().quadric += profile.v1().quadric;
  }
};

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Usage " << argv[0] << " input.off output.off N to load the file input.off, contract N edges and save "
        << "the result to output.off.\n";
    return EXIT_FAILURE;
  }

  Complex complex;

  // load the points
  Skeleton_blocker_off_reader<Complex> off_reader(argv[1], complex);
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read file:" << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Load complex with " << complex.num_vertices() << " vertices" << std::endl;


  int num_contractions = atoi(argv[3]);

  boost::timer::auto_cpu_timer t;

  // constructs the contractor object with Garland Heckbert policies.
  Complex_contractor contractor(complex,
                                new GH_cost(complex),
                                new GH_placement(complex),
                                contraction::make_link_valid_contraction<EdgeProfile>(),
                                new GH_visitor(complex));

  std::cout << "Contract " << num_contractions << " edges" << std::endl;
  contractor.contract_edges(num_contractions);

  std::cout << "Final complex has " <<
      complex.num_vertices() << " vertices, " <<
      complex.num_edges() << " edges and" <<
      complex.num_triangles() << " triangles." << std::endl;

  // write simplified complex
  Skeleton_blocker_off_writer<Complex> off_writer(argv[2], complex);

  return EXIT_SUCCESS;
}



#endif  // GARLAND_HECKBERT_H_
