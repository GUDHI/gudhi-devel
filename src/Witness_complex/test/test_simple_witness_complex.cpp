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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simple_witness_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Witness_complex.h>
#include <gudhi/Strong_witness_complex.h>

#include <iostream>
#include <ctime>
#include <vector>

typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_d Point_d;
typedef Gudhi::witness_complex::Witness_complex<Kernel> WitnessComplex;
typedef Gudhi::witness_complex::Strong_witness_complex<Kernel> StrongWitnessComplex;

/* All landmarks and witnesses are taken on the grid in the following manner.
   LWLWL  
   WW.WW  
   L...L  
   WW.WW  
   LWLWL  

   Witness complex consists of 8 vertices, 12 edges and 4 triangles
 */

BOOST_AUTO_TEST_CASE(simple_witness_complex) {
  Simplex_tree complex, relaxed_complex, strong_relaxed_complex;

  std::vector<Point_d> witnesses, landmarks;

  landmarks.push_back(Point_d(std::vector<FT>{-2,-2}));
  landmarks.push_back(Point_d(std::vector<FT>{-2, 0}));
  landmarks.push_back(Point_d(std::vector<FT>{-2, 2}));
  landmarks.push_back(Point_d(std::vector<FT>{ 0,-2}));
  landmarks.push_back(Point_d(std::vector<FT>{ 0, 2}));
  landmarks.push_back(Point_d(std::vector<FT>{ 2,-2}));
  landmarks.push_back(Point_d(std::vector<FT>{ 2, 0}));
  landmarks.push_back(Point_d(std::vector<FT>{ 2, 2}));
  witnesses.push_back(Point_d(std::vector<FT>{-2,-1}));
  witnesses.push_back(Point_d(std::vector<FT>{-2, 1}));
  witnesses.push_back(Point_d(std::vector<FT>{-1,-2}));
  witnesses.push_back(Point_d(std::vector<FT>{-1,-1}));
  witnesses.push_back(Point_d(std::vector<FT>{-1, 1}));
  witnesses.push_back(Point_d(std::vector<FT>{-1, 2}));
  witnesses.push_back(Point_d(std::vector<FT>{ 1,-2}));
  witnesses.push_back(Point_d(std::vector<FT>{ 1,-1}));
  witnesses.push_back(Point_d(std::vector<FT>{ 1, 1}));
  witnesses.push_back(Point_d(std::vector<FT>{ 1, 2}));
  witnesses.push_back(Point_d(std::vector<FT>{ 2,-1}));
  witnesses.push_back(Point_d(std::vector<FT>{ 2, 1}));
  
  WitnessComplex witness_complex(landmarks.begin(),
                                 landmarks.end(),
                                 witnesses.begin(),
                                 witnesses.end());
  witness_complex.create_complex(complex, 0);

  std::cout << "complex.num_simplices() = " << complex.num_simplices() << std::endl; 
  BOOST_CHECK(complex.num_simplices() == 24);

  witness_complex.create_complex(relaxed_complex, 8.01);

  std::cout << "relaxed_complex.num_simplices() = " << relaxed_complex.num_simplices() << std::endl; 
  BOOST_CHECK(relaxed_complex.num_simplices() == 239);
  
  StrongWitnessComplex strong_witness_complex(landmarks.begin(),
                                              landmarks.end(),
                                              witnesses.begin(),
                                              witnesses.end());

  strong_witness_complex.create_complex(strong_relaxed_complex, 9.1);
    
  std::cout << "strong_relaxed_complex.num_simplices() = " << strong_relaxed_complex.num_simplices() << std::endl; 
  BOOST_CHECK(strong_relaxed_complex.num_simplices() == 239);
  
}
