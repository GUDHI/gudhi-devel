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
#include <algorithm>
#include <set>
#include <queue>
#include <iterator>

#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>

//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Relaxed_witness_complex.h"
#include "gudhi/reader_utils.h"
#include "gudhi/Collapse/Collapse.h"
//#include <boost/filesystem.hpp> 

//#include <CGAL/Delaunay_triangulation.h>
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

using namespace Gudhi;
//using namespace boost::filesystem;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::FT FT;
typedef K::Point_d Point_d;
typedef CGAL::Search_traits<
  FT, Point_d,
  typename K::Cartesian_const_iterator_d,
  typename K::Construct_cartesian_const_iterator_d> Traits_base;
typedef CGAL::Euclidean_distance<Traits_base> Euclidean_distance;

typedef std::vector< Vertex_handle > typeVectorVertex;

//typedef std::pair<typeVectorVertex, Filtration_value> typeSimplex;
//typedef std::pair< Simplex_tree<>::Simplex_handle, bool > typePairSimplexBool;

typedef CGAL::Search_traits_adapter<
  std::ptrdiff_t, Point_d*, Traits_base> STraits;
//typedef K TreeTraits;
//typedef CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance > Euclidean_adapter;
//typedef CGAL::Kd_tree<STraits> Kd_tree;
typedef CGAL::Orthogonal_incremental_neighbor_search<STraits, CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>> Neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef Neighbor_search::Distance Distance;
typedef Neighbor_search::iterator KNS_iterator;
typedef Neighbor_search::iterator KNS_range;
typedef boost::container::flat_map<int, int> Point_etiquette_map;
typedef CGAL::Kd_tree<STraits> Tree2;

typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;

typedef std::vector<Point_d> Point_Vector;

//typedef K::Equal_d Equal_d;
typedef CGAL::Random_points_in_cube_d<Point_d> Random_cube_iterator;
typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;

bool toric=false;

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


void generate_points_sphere(Point_Vector& W, int nbP, int dim)
{
  CGAL::Random_points_on_sphere_d<Point_d> rp(dim,1);
  for (int i = 0; i < nbP; i++)
    W.push_back(*rp++);
}


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

void write_rl( std::string file_name, std::vector< std::vector <std::vector<int>::iterator> > & rl)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  for (auto w : rl)
    {
      for (auto l: w)
        ofs << *l << " ";
      ofs << "\n";
    }
  ofs.close();
}

std::vector<Point_d> convert_to_torus(std::vector< Point_d>& points)
{
  std::vector< Point_d > points_torus;
  for (auto p: points)
    {
      FT theta = M_PI*p[0];
      FT phi = M_PI*p[1];
      std::vector<FT> p_torus;
      p_torus.push_back((1+0.2*cos(theta))*cos(phi));
      p_torus.push_back((1+0.2*cos(theta))*sin(phi));
      p_torus.push_back(0.2*sin(theta));
      points_torus.push_back(Point_d(p_torus));
    }
  return points_torus;
}


void write_points_torus( std::string file_name, std::vector< Point_d > & points)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  std::vector<Point_d> points_torus = convert_to_torus(points);
  for (auto w : points_torus)
    {
      for (auto it = w.cartesian_begin(); it != w.cartesian_end(); ++it)
        ofs << *it << " ";
      ofs << "\n";
    }
  ofs.close();
}


void write_points( std::string file_name, std::vector< Point_d > & points)
{
  if (toric) write_points_torus(file_name, points);
  else
    {
      std::ofstream ofs (file_name, std::ofstream::out);
      for (auto w : points)
        {
          for (auto it = w.cartesian_begin(); it != w.cartesian_end(); ++it)
            ofs << *it << " ";
          ofs << "\n";
        }
      ofs.close();
    }
}


void write_edges_torus(std::string file_name, Witness_complex<>& witness_complex, Point_Vector& landmarks)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  Point_Vector l_torus = convert_to_torus(landmarks);
  for (auto u: witness_complex.complex_vertex_range())
    for (auto v: witness_complex.complex_vertex_range())
      {
        typeVectorVertex edge = {u,v};
        if (u < v && witness_complex.find(edge) != witness_complex.null_simplex())   
          {
            for (auto it = l_torus[u].cartesian_begin(); it != l_torus[u].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n";
            for (auto it = l_torus[v].cartesian_begin(); it != l_torus[v].cartesian_end(); ++it)
              ofs << *it << " ";
            ofs << "\n\n\n";
          }
    }
  ofs.close();
}

void write_edges(std::string file_name, Witness_complex<>& witness_complex, Point_Vector& landmarks)
{
  std::ofstream ofs (file_name, std::ofstream::out);
  if (toric) write_edges_torus(file_name, witness_complex, landmarks);
  else
    {
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


void landmarks_to_witness_complex(Point_Vector &W, Point_Vector& landmarks, std::vector<int>& landmarks_ind, FT alpha)
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
  int nbL = landmarks.size();
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
  for (int i = 0; i < nbP; i++)
    {
      //std::cout << "Entered witness number " << i << std::endl;
      Point_d& w = W[i];
      std::queue< typename std::vector<int>::iterator > ope_queue; // queue of points at (1+epsilon) distance to current landmark
      Neighbor_search search(L, w, FT(0), true, CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks[0])));
      Neighbor_search::iterator search_it = search.begin();
      
      //Incremental search and filling WL
      while (WL[i].size() < D)
        WL[i].push_back((search_it++)->first);
      FT dtow = ed.transformed_distance(w, landmarks[WL[i][D-1]]);
      while (search_it->second < dtow + alpha)
        WL[i].push_back((search_it++)->first);

      //Filling the (1+epsilon)-limits table
      for (std::vector<int>::iterator wl_it = WL[i].begin(); wl_it != WL[i].end(); ++wl_it)
        {
          ope_queue.push(wl_it);
          FT d_to_curr_l = ed.transformed_distance(w, landmarks[*wl_it]);
          //std::cout << "d_to_curr_l=" << d_to_curr_l << std::endl;
          //std::cout << "d_to_front+alpha=" << d_to_curr_l << std::endl;
          while (d_to_curr_l > alpha + ed.transformed_distance(w, landmarks[*(ope_queue.front())]))
            {
              ope_limits[i].push_back(wl_it);
              ope_queue.pop();
            }
        }
      while (ope_queue.size() > 0)
        {
          ope_limits[i].push_back(WL[i].end());
          ope_queue.pop();
        }
      //std::cout << "Safely constructed a point\n";
      ////Search D+1 nearest neighbours from the tree of landmarks L
      /*
      if (w[0]>0.95)
        std::cout << i << std::endl;
      */
      //K_neighbor_search search(L, w, D, FT(0), true,
      //                         CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks[0])) );
      //std::cout << "Safely found nearest landmarks\n";
      /* 
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
      */
    }
  //std::cout << "\n";
  
  //std::string out_file = "wl_result";
  write_wl("wl_result",WL);
  write_rl("rl_result",ope_limits);
  
  //******************** Constructng a witness complex
  std::cout << "Entered witness complex construction\n";
  Witness_complex<> witnessComplex;
  witnessComplex.setNbL(nbL);
  witnessComplex.relaxed_witness_complex(WL, ope_limits);
  char buffer[100];
  int i = sprintf(buffer,"stree_result.txt");
  
  if (i >= 0)
    {
      std::string out_file = (std::string)buffer;
      std::ofstream ofs (out_file, std::ofstream::out);
      witnessComplex.st_to_file(ofs);
      ofs.close();
    }
  write_edges("landmarks/edges", witnessComplex, landmarks);
  std::cout << Distance().transformed_distance(Point_d(std::vector<double>({0.1,0.1})), Point_d(std::vector<double>({1.9,1.9}))) << std::endl;
}


int main (int argc, char * const argv[])
{
  
  if (argc != 5)
    {
      std::cerr << "Usage: " << argv[0]
                << " nbP nbL dim alpha\n";
      return 0;
    }
  /*
  boost::filesystem::path p;
  for (; argc > 2; --argc, ++argv)
    p /= argv[1];
  */
  
  int nbP       = atoi(argv[1]);
  int nbL       = atoi(argv[2]);
  int dim       = atoi(argv[3]);
  double alpha  = atof(argv[4]);
  //clock_t start, end;
  //Construct the Simplex Tree
  Witness_complex<> witnessComplex;
 
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  //read_points_cust(file_name, point_vector);
  generate_points_sphere(point_vector, nbP, dim);
  /*
  for (auto &p: point_vector)
    {
      assert(std::count(point_vector.begin(),point_vector.end(),p) == 1);
    }
  */
  //std::cout << "Successfully read the points\n";
  //witnessComplex.setNbL(nbL);
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  landmark_choice(point_vector, nbP, nbL, L, chosen_landmarks);
  //start = clock();
  
  write_points("landmarks/initial_pointset",point_vector);
  write_points("landmarks/initial_landmarks",L);
  
  landmarks_to_witness_complex(point_vector, L, chosen_landmarks, alpha);
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
