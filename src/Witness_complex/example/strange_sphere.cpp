/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2015  INRIA Sophia Antipolis-Méditerranée (France)
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
#define BOOST_PARAMETER_MAX_ARITY 12


#include <sys/types.h>
#include <sys/stat.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Strange_witness_complex.h>
#include <gudhi/Landmark_choice_by_random_point.h>
#include <gudhi/reader_utils.h>
#include "output.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <string>
#include <vector>

#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>

#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Random.h>


#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>


#include "generators.h"

using namespace Gudhi;
using namespace Gudhi::witness_complex;

typedef std::vector< Vertex_handle > typeVectorVertex;
typedef Strange_witness_complex< Simplex_tree<> > WitnessComplex;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits<
  FT, Point_d,
  typename K::Cartesian_const_iterator_d,
  typename K::Construct_cartesian_const_iterator_d> Traits_base;
typedef CGAL::Euclidean_distance<Traits_base> Euclidean_distance;

typedef CGAL::Search_traits_adapter<
  std::ptrdiff_t, Point_d*, Traits_base> STraits;
typedef CGAL::Orthogonal_incremental_neighbor_search<STraits, CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Neighbor_search::Distance Distance;
typedef Neighbor_search::iterator KNS_iterator;
typedef Neighbor_search::iterator KNS_range;
typedef boost::container::flat_map<int, int> Point_etiquette_map;
typedef CGAL::Kd_tree<STraits> Tree2;
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;
typedef std::vector<Point_d> Point_Vector;

/** Function that chooses landmarks from W and place it in the kd-tree L.
 *  Note: nbL hould be removed if the code moves to Witness_complex
 */
void landmark_choice(Point_Vector &W, int nbP, int nbL, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  std::cout << "Enter landmark choice to kd tree\n";
  //std::vector<Point_d> landmarks;
  int chosen_landmark;
  //std::pair<Point_etiquette_map::iterator,bool> res = std::make_pair(L_i.begin(),false);
  Point_d* p;
  CGAL::Random rand;
  for (int i = 0; i < nbL; i++)
    {
      //      while (!res.second)
      //  {
      do chosen_landmark = rand.get_int(0,nbP);
      while (std::find(landmarks_ind.begin(), landmarks_ind.end(), chosen_landmark) != landmarks_ind.end());
      //rand++;
      //std::cout << "Chose " << chosen_landmark << std::endl;
      p = &W[chosen_landmark];
      //L_i.emplace(chosen_landmark,i);
      //  }
      landmarks.push_back(*p);
      landmarks_ind.push_back(chosen_landmark);
      //std::cout << "Added landmark " << chosen_landmark << std::endl;
    }
 }

void landmarks_to_witness_complex(Point_Vector &W, int nbL)
{
  //********************Preface: origin point
  unsigned D = W[0].size();
  std::vector<FT> orig_vector;
  for (unsigned i = 0; i < D; i++)
    orig_vector.push_back(0);
  Point_d origin(orig_vector);
  //Distance dist;
  //dist.transformed_distance(0,1);
  //******************** Constructing a WL matrix
  int nbP = W.size();
  Point_Vector landmarks;
  std::vector<int> landmarks_ind;
  landmark_choice(W, nbP, nbL, landmarks, landmarks_ind);
  
  
  STraits traits(&(landmarks[0]));
  Euclidean_distance ed;
  std::vector< std::vector <int> > WL(nbP);
  std::vector< std::vector< typename std::vector<int>::iterator > > ope_limits(nbP);
  Tree L(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nbL),
         typename Tree::Splitter(),
         traits);

  std::cout << "Enter (D+1) nearest landmarks\n";
  //std::cout << "Size of the tree is " << L.size() << std::endl;
  for (int i = 0; i < nbP; i++) {
    //std::cout << "Entered witness number " << i << std::endl;
    Point_d& w = W[i];
    std::queue< typename std::vector<int>::iterator > ope_queue; // queue of points at (1+epsilon) distance to current landmark
    Neighbor_search search(L, w, FT(0), true, CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks[0])));
    Neighbor_search::iterator search_it = search.begin();
      
    //Incremental search and filling WL
    while (WL[i].size() < D)
      WL[i].push_back((search_it++)->first);
  }
  write_wl("WL.txt", WL); //!
  //******************** Constructng a witness complex
  std::cout << "Entered witness complex construction\n";

  Simplex_tree<> simplex_tree;
  WitnessComplex(WL, simplex_tree, nbL, W[0].size()-1);
  //std::cout << simplex_tree << "\n"; //!
  write_witness_mesh(W, landmarks_ind, simplex_tree, false, true);
}



/** Write a gnuplot readable file.
 *  Data range is a random access range of pairs (arg, value)
 */
template < typename Data_range >
void write_data(Data_range & data, std::string filename) {
  std::ofstream ofs(filename, std::ofstream::out);
  for (auto entry : data)
    ofs << entry.first << ", " << entry.second << "\n";
  ofs.close();
}



int main(int argc, char * const argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
        << " nbP nbL dim\n";
    return 0;
  }

  int nbP = atoi(argv[1]);
  int nbL = atoi(argv[2]);
  int dim = atoi(argv[3]);
  clock_t start, end;

  // Construct the Simplex Tree
  Simplex_tree<> simplex_tree;

  std::vector< std::pair<int, double> > l_time;

  // Read the point file
  Point_Vector point_vector;
  generate_points_sphere(point_vector, nbP, dim);
  std::cout << "Successfully generated " << point_vector.size() << " points.\n";
  std::cout << "Ambient dimension is " << point_vector[0].size() << ".\n";

  // Choose landmarks
  start = clock();
  std::vector<std::vector< int > > knn;
  
  landmarks_to_witness_complex(point_vector, nbL);

  end = clock();
  double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
  std::cout << "Witness complex for " << nbL << " landmarks took "
            << time << " s. \n";
  l_time.push_back(std::make_pair(nbP, time));
  //write_data(l_time, "w_time.dat");
}
