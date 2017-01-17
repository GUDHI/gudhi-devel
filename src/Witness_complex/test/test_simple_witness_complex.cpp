#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "simple_witness_complex"
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Epick_d.h>

#include <gudhi/Simplex_tree.h>

#include <gudhi/Witness_complex.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/Strong_witness_complex.h>
#include <gudhi/Euclidean_strong_witness_complex.h>

#include <gudhi/Kd_tree_search.h>

#include <iostream>
#include <ctime>
#include <vector>

typedef Gudhi::Simplex_tree<> Simplex_tree;
typedef std::vector< Vertex_handle > typeVectorVertex;
typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_d Point_d;
typedef Gudhi::witness_complex::Euclidean_witness_complex<Kernel> EuclideanWitnessComplex;
typedef Gudhi::witness_complex::Euclidean_strong_witness_complex<Kernel> EuclideanStrongWitnessComplex;

typedef std::vector<Point_d> Point_range;
typedef Gudhi::spatial_searching::Kd_tree_search<Kernel, Point_range> Kd_tree;
typedef Kd_tree::INS_range Nearest_landmark_range; 
typedef std::vector<Nearest_landmark_range> Nearest_landmark_table;
typedef Gudhi::witness_complex::Witness_complex<Nearest_landmark_table> WitnessComplex;


/* All landmarks and witnesses are taken on the grid in the following manner.
   LWLWL  
   WW.WW  
   L...L  
   WW.WW  
   LWLWL  

   Witness complex consists of 8 vertices, 12 edges and 4 triangles
 */

BOOST_AUTO_TEST_CASE(simple_witness_complex) {
  Simplex_tree complex, relaxed_complex, strong_relaxed_complex, strong_relaxed_complex2;
  Simplex_tree complex_ne, relaxed_complex_ne, strong_relaxed_complex_ne, strong_relaxed_complex2_ne;
    
  Point_range witnesses, landmarks;

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

  Kd_tree landmark_tree(landmarks);
  Nearest_landmark_table nearest_landmark_table;
  for (auto w: witnesses)
    nearest_landmark_table.push_back(landmark_tree.query_incremental_nearest_neighbors(w));

  // Weak witness complex: Euclidean version
  EuclideanWitnessComplex eucl_witness_complex(landmarks,
                                               witnesses);
  eucl_witness_complex.create_complex(complex, 0);

  std::cout << "complex.num_simplices() = " << complex.num_simplices() << std::endl; 
  BOOST_CHECK(complex.num_simplices() == 24);

  eucl_witness_complex.create_complex(relaxed_complex, 8.01);

  std::cout << "relaxed_complex.num_simplices() = " << relaxed_complex.num_simplices() << std::endl; 
  BOOST_CHECK(relaxed_complex.num_simplices() == 239);
  // The corner simplex {0,2,5,7} and its cofaces are missing.

  // Weak witness complex: non-Euclidean version
  WitnessComplex witness_complex(nearest_landmark_table);
  witness_complex.create_complex(complex_ne, 0);

  std::cout << "complex.num_simplices() = " << complex_ne.num_simplices() << std::endl; 
  BOOST_CHECK(complex_ne.num_simplices() == 24);

  witness_complex.create_complex(relaxed_complex_ne, 8.01);

  std::cout << "relaxed_complex.num_simplices() = " << relaxed_complex_ne.num_simplices() << std::endl; 
  BOOST_CHECK(relaxed_complex_ne.num_simplices() == 239);
    
  
  
  EuclideanStrongWitnessComplex strong_witness_complex(landmarks,
                                              witnesses);

  strong_witness_complex.create_complex(strong_relaxed_complex, 9.1);
  strong_witness_complex.create_complex(strong_relaxed_complex2, 9.1, 2);
  
  std::cout << "strong_relaxed_complex.num_simplices() = " << strong_relaxed_complex.num_simplices() << std::endl; 
  BOOST_CHECK(strong_relaxed_complex.num_simplices() == 239);

  std::cout << "strong_relaxed_complex2.num_simplices() = " << strong_relaxed_complex2.num_simplices() << std::endl;
  BOOST_CHECK(strong_relaxed_complex2.num_simplices() == 92);
  // 8 vertices, 28 edges, 56 triangles
}
