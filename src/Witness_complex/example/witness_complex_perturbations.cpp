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
#include <set>
#include <iterator>

#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>

//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Witness_complex.h"
#include "gudhi/reader_utils.h"
//#include <boost/filesystem.hpp> 

//#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>

#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Origin.h>

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
typedef CGAL::Kd_tree<STraits> Tree2;

typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;

typedef std::vector<Point_d> Point_Vector;

typedef CGAL::Euclidean_distance<Traits_base> Euclidean_distance;
//typedef K::Equal_d Equal_d;
typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;
/**
 * \brief Customized version of read_points
 * which takes into account a possible nbP first line
 *
 */
inline void
read_points_cust ( std::string file_name , Point_Vector & points)
{  
  std::ifstream in_file (file_name.c_str(),std::ios::in);
  if(!in_file.is_open())
    {
      std::cerr << "Unable to open file " << file_name << std::endl;
      return;
    }
  std::string line;
  double x;
  while( getline ( in_file , line ) )
    {
      std::vector< double > point;
      std::istringstream iss( line );
      while(iss >> x) { point.push_back(x); }
      Point_d p(point.begin(), point.end());
      if (point.size() != 1)
        points.push_back(p);
    }
  in_file.close();
}

/*
void read_points_to_tree (std::string file_name, Tree& tree)
{
  //I assume here that tree is empty
  std::ifstream in_file (file_name.c_str(),std::ios::in);
  if(!in_file.is_open())
    {
      std::cerr << "Unable to open file " << file_name << std::endl;
      return;
    }
  std::string line;
  double x;
   while( getline ( in_file , line ) )
    {
      std::vector<double> coords;
      std::istringstream iss( line );
      while(iss >> x) { coords.push_back(x); }
      if (coords.size() != 1)
        {
          Point_d point(coords.begin(), coords.end());
          tree.insert(point);
        }
    }
  in_file.close();
}
*/

void write_wl( std::string file_name, std::vector< std::vector <int> > & WL)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : WL)
    {
      for (auto l: w)
        ofs << l << " ";
      ofs << "\n";
    }
  ofs.close();
}

void write_points( std::string file_name, std::vector< Point_d > & WL)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : WL)
    {
      for (auto it = w.cartesian_begin(); it != w.cartesian_end(); ++it)
        ofs << *it << " ";
      ofs << "\n";
    }
  ofs.close();
}

void write_edges_gnuplot(std::string file_name, Witness_complex<>& witness_complex, Point_Vector& landmarks)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto u: witness_complex.complex_vertex_range())
    for (auto v: witness_complex.complex_vertex_range())
      {
        typeVectorVertex edge = {u,v};
        if (u < v && witness_complex.find(edge) != witness_complex.null_simplex())   
          {
            for (auto it = landmarks[u].cartesian_begin(); it != landmarks[u].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n";
            for (auto it = landmarks[v].cartesian_begin(); it != landmarks[v].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n\n\n";
          }
    }
  ofs.close();
}
    


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
      landmarks_ind.push_back(chosen_landmark);
      //std::cout << "Added landmark " << chosen_landmark << std::endl;
    }
 }
/*
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
      while (std::count(landmarks_ind.begin(),landmarks_ind.end(),chosen_landmark)!=0)
        chosen_landmark = rand.get_int(0,nbP);
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
*/

int landmark_perturbation(Point_Vector &W, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  //******************** Constructing a WL matrix
  int nbP = W.size();
  int nbL = landmarks.size();
  //Point_Vector landmarks_ = landmarks;
  Euclidean_distance ed;
  //Equal_d ed;
  FT lambda = ed.transformed_distance(landmarks[0],landmarks[1]);
    //FT lambda = 0.1;//Euclidean_distance();
  std::vector< std::vector <int> > WL(nbP);
  Tree L(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nbL),
         typename Tree::Splitter(),
         STraits(&(landmarks[0])));
  /*Tree2 L2(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(nbL),
           typename Tree::Splitter(),
           STraits(&(landmarks[0])));
  */
  std::cout << "Enter (D+1) nearest landmarks\n";
  //std::cout << "Size of the tree is " << L.size() << std::endl;
  int D = W[0].size();
  for (int i = 0; i < nbP; i++)
    {
      //std::cout << "Entered witness number " << i << std::endl;
      Point_d& w = W[i];
      //std::cout << "Safely constructed a point\n";
      ////Search D+1 nearest neighbours from the tree of landmarks L 
      K_neighbor_search search(L, w, D+1, FT(0), true,
                               CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,CGAL::Euclidean_distance<Traits_base>>(&(landmarks[0])) );
      //std::cout << "Safely found nearest landmarks\n";
      for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
        {
          //std::cout << "Entered KNN_it with point at distance " << it->second << "\n";
          //Point_etiquette_map::iterator itm = L_i.find(it->first);
          //assert(itm != L_i.end());
          //std::cout << "Entered KNN_it with point at distance " << it->second << "\n";
          WL[i].push_back(it->first);
          //std::cout << "ITFIRST " << it->first << std::endl;
          //std::cout << i << " " << it->first << ": " << it->second << std::endl; 
        }
      if (i == landmarks_ind[WL[i][0]])
        {
          //std::cout << "'";
          FT dist = ed.transformed_distance(W[i], landmarks[WL[i][1]]);
          if (dist < lambda)
            lambda = dist;
        }
    }
  //std::cout << "\n";
  
  std::string out_file = "wl_result";
  write_wl(out_file,WL);
  
  //******************** Constructng a witness complex
  std::cout << "Entered witness complex construction\n";
  Witness_complex<> witnessComplex;
  witnessComplex.setNbL(nbL);
  witnessComplex.witness_complex(WL);
  //******************** Making a set of bad link landmarks
  std::cout << "Entered bad links\n";
  std::set< int > perturbL;
  int count_badlinks = 0;
  std::cout << "Bad links around ";
  for (auto u: witnessComplex.complex_vertex_range())
    if (!witnessComplex.has_good_link(u))
      {
        //std::cout << "Landmark " << u << " start!" << std::endl;
        //perturbL.insert(u);
        count_badlinks++;
        std::cout << u << " ";
        Point_d& l = landmarks[u];
        Fuzzy_sphere fs(l, sqrt(lambda)*2, 0, STraits(&(landmarks[0])));
        L.search(std::insert_iterator<std::set<int>>(perturbL,perturbL.begin()),fs);
        //L.search(std::inserter(perturbL,perturbL.begin()),fs);
        //L.search(std::ostream_iterator<int>(std::cout,"\n"),fs);
        //std::cout << "PerturbL size is " << perturbL.size() << std::endl;
      }
  std::cout << "\nBad links total: " << count_badlinks << " Points to perturb: " << perturbL.size() << std::endl;
  //std::cout << "landmark[0][0] before" << landmarks[0][0] << std::endl;
  //*********************** Perturb bad link landmarks
  
  for (auto u: perturbL)
    {
      Random_point_iterator rp(D,sqrt(lambda)/2);
      //std::cout << landmarks[u] << std::endl;
      
      std::vector<FT> point;
      for (int i = 0; i < D; i++)
        {
          point.push_back(W[landmarks_ind[u]][i] + (*rp)[i]);
        }
      landmarks[u] = Point_d(point);
      //std::cout << landmarks[u] << std::endl;
    }
  
  //std::cout << "landmark[0][0] after" << landmarks[0][0] << std::endl;
  std::cout << "lambda=" << lambda << std::endl;
  // Write the WL matrix in a file
 
  
  char buffer[100];
  int i = sprintf(buffer,"stree_result.txt");
  
  if (i >= 0)
    {
      std::string out_file = (std::string)buffer;
      std::ofstream ofs (out_file, std::ofstream::out);
      witnessComplex.st_to_file(ofs);
      ofs.close();
    }
  //witnessComplex.write_badlinks("badlinks");
  write_edges_gnuplot("landmarks/edges", witnessComplex, landmarks);
  return count_badlinks;
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
  
  //clock_t start, end;
  //Construct the Simplex Tree
  //Witness_complex<> witnessComplex;
 
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  read_points_cust(file_name, point_vector);
  //std::cout << "Successfully read the points\n";
  //witnessComplex.setNbL(nbL);
  //  witnessComplex.witness_complex_from_points(point_vector);
  int nbP = point_vector.size();
  //std::vector<std::vector< int > > WL(nbP);
  //std::set<int> L;
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  //Point_etiquette_map L_i;
  //start = clock();
  //witnessComplex.landmark_choice_by_furthest_points(point_vector, point_vector.size(), WL);
  landmark_choice(point_vector, nbP, nbL, L, chosen_landmarks);
  int bl = 1;

  mkdir("landmarks", S_IRWXU);
  const size_t last_slash_idx = file_name.find_last_of("/");
  if (std::string::npos != last_slash_idx)
    {
      file_name.erase(0, last_slash_idx + 1);
    }
  write_points("landmarks/initial_pointset",point_vector);
  write_points("landmarks/initial_landmarks",L);
  for (int i = 0; bl != 0; i++)
    {
      std::cout << "========== Start iteration " << i << " ========\n";
      bl = landmark_perturbation(point_vector, L, chosen_landmarks);
      std::ostringstream os(std::ostringstream::ate);;
      os << "landmarks/landmarks0";
      write_points(os.str(),L);
    }
  //end = clock();
  
  /*
  std::cout << "Landmark choice took "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  start = clock();
  witnessComplex.witness_complex(WL);
  //
  end = clock();
  std::cout << "Howdy world! The process took "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  */
  
  /*
  out_file = "output/"+file_name+"_"+argv[2]+".stree";
  std::ofstream ofs (out_file, std::ofstream::out);
  witnessComplex.st_to_file(ofs);
  ofs.close();

  out_file = "output/"+file_name+"_"+argv[2]+".badlinks";
  std::ofstream ofs2(out_file, std::ofstream::out);
  witnessComplex.write_bad_links(ofs2);
  ofs2.close();
  */
}
