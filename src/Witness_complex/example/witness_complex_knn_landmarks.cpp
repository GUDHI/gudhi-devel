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

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>

#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>

//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Witness_complex.h"
#include "gudhi/reader_utils.h"
#include "generators.h"
#include "output.h"
//#include <boost/filesystem.hpp>

//#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

using namespace Gudhi;
//using namespace boost::filesystem;

typedef std::vector< Vertex_handle > typeVectorVertex;

//typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
//typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits<
  FT, Point_d,
  typename K::Cartesian_const_iterator_d,
  typename K::Construct_cartesian_const_iterator_d> Traits_base;
typedef CGAL::Search_traits_adapter<
  std::ptrdiff_t, Point_d*, Traits_base> STraits;
//typedef K TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<STraits> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef K_neighbor_search::Distance Distance;
typedef K_neighbor_search::iterator KNS_iterator;
typedef K_neighbor_search::iterator KNS_range;
typedef boost::container::flat_map<int, int> Point_etiquette_map;

typedef std::vector<Point_d> Point_Vector;

/** Function that chooses landmarks from W and place it in the kd-tree L.
 *  Note: nbL hould be removed if the code moves to Witness_complex
 */
void landmark_choice_to_tree(Point_Vector &W, int nbP, Point_etiquette_map &L_i, int nbL, std::vector< std::vector <int> > &WL)
{
  std::cout << "Enter landmark choice to kd tree\n";
  std::vector<Point_d> landmarks;
  int chosen_landmark;
  //std::pair<Point_etiquette_map::iterator,bool> res = std::make_pair(L_i.begin(),false);
  Point_d* p;
  srand(24660);
  for (int i = 0; i < nbL; i++)
    {
      //      while (!res.second)
      //  {
          chosen_landmark = rand()%nbP;
          p = &W[chosen_landmark];
          //L_i.emplace(chosen_landmark,i);
          //  }
      landmarks.push_back(*p);
      //std::cout << "Added landmark " << chosen_landmark << std::endl;
    }
  Tree L(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nbL),
         typename Tree::Splitter(),
         STraits((Point_d*)&(landmarks[0])));
  /*}


void d_nearest_landmarks(Point_Vector &W, Tree &L, Point_etiquette_map &L_i, std::vector< std::vector <int> > &WL)
{*/
  std::cout << "Enter (D+1) nearest landmarks\n";
  std::cout << "Size of the tree is " << L.size() << std::endl;
//int nbP = W.size();
  int D = W[0].size();
  for (int i = 0; i < nbP; i++)
    {
      //std::cout << "Entered witness number " << i << std::endl;
      Point_d& w = W[i];
      //std::cout << "Safely constructed a point\n";
      //Search D+1 nearest neighbours from the tree of landmarks L 
      K_neighbor_search search(L, w, D+1, FT(0), true,
                               CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,CGAL::Euclidean_distance<Traits_base>>((Point_d*)&(landmarks[0])) );
      //std::cout << "Safely found nearest landmarks\n";
      for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
        {
          //std::cout << "Entered KNN_it with point at distance " << it->second << "\n";
          //Point_etiquette_map::iterator itm = L_i.find(it->first);
          //assert(itm != L_i.end());
          //std::cout << "Entered KNN_it with point at distance " << it->second << "\n";
          WL[i].push_back(it->first);
          //std::cout << i << " " << it->first << ": " << it->second << std::endl; 
        }
    }
}

int main (int argc, char * const argv[])
{
  if (argc != 3)
    {
      std::cerr << "Usage: " << argv[0]
                << " path_to_point_file nbL \n";
      return 0;
    }
  /*
  boost::filesystem::path p;

  for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  std::string file_name   = argv[1];
  int nbL       = atoi(argv[2]);
  
  clock_t start, end;
  //Construct the Simplex Tree
  Witness_complex<> witnessComplex;
 
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  //std::cout << "Successfully read the points\n";
  witnessComplex.setNbL(nbL);
  //  witnessComplex.witness_complex_from_points(point_vector);
  int nbP = point_vector.size();
  std::vector<std::vector< int > > WL(nbP);
  //std::set<int> L;
  Tree L;
  Point_etiquette_map L_i;
  start = clock();
  //witnessComplex.landmark_choice_by_furthest_points(point_vector, point_vector.size(), WL);
  landmark_choice_to_tree(point_vector, nbP, L_i, nbL, WL);
  //d_nearest_landmarks(point_vector, L, L_i, WL);
  end = clock();
  std::cout << "Landmark choice took "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  // Write the WL matrix in a file
  mkdir("output", S_IRWXU);
  const size_t last_slash_idx = file_name.find_last_of("/");
  if (std::string::npos != last_slash_idx)
    {
      file_name.erase(0, last_slash_idx + 1);
    }
  std::string out_file = "output/"+file_name+"_"+argv[2]+".wl";
  write_wl(out_file,WL);
  start = clock();
  witnessComplex.witness_complex(WL);
  //
  end = clock();
  std::cout << "Howdy world! The process took "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  /*
  char buffer[100];
  int i = sprintf(buffer,"%s_%s_result.txt",argv[1],argv[2]);
  if (i >= 0)
    {
      std::string out_file = (std::string)buffer;
      std::ofstream ofs (out_file, std::ofstream::out);
      witnessComplex.st_to_file(ofs);
      ofs.close();
    }
  */

  out_file = "output/"+file_name+"_"+argv[2]+".stree";
  std::ofstream ofs (out_file, std::ofstream::out);
  witnessComplex.st_to_file(ofs);
  ofs.close();

  out_file = "output/"+file_name+"_"+argv[2]+".badlinks";
  std::ofstream ofs2(out_file, std::ofstream::out);
  //witnessComplex.write_bad_links(ofs2);
  ofs2.close();
}
