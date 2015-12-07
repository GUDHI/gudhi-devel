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
#include <iterator>

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
typedef CGAL::Orthogonal_k_neighbor_search<STraits, CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef K_neighbor_search::Distance Distance;
typedef K_neighbor_search::iterator KNS_iterator;
typedef K_neighbor_search::iterator KNS_range;
typedef boost::container::flat_map<int, int> Point_etiquette_map;
typedef CGAL::Kd_tree<STraits> Tree2;

typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;

typedef std::vector<Point_d> Point_Vector;

//typedef K::Equal_d Equal_d;

bool toric=false;

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

/** \brief A test with 600cell, the generalisation of icosaedre in 4d
 */
void landmark_choice_600cell(Point_Vector&W, int nbP, int nbL, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  assert(W[0].size() == 4); //4-dimensionality required
  FT phi = (1+sqrt(5))/2;
  FT phi_1 = FT(1)/phi;
  std::vector<FT> p;
  // 16 vertices
  for (FT a = -0.5; a < 1; a += 1)
    for (FT b = -0.5; b < 1; b += 1)
      for (FT c = -0.5; c < 1; c += 1)
        for (FT d = -0.5; d < 1; d += 1)
          landmarks.push_back(Point_d(std::vector<FT>({a,b,c,d})));
  // 8 vertices
  for (FT a = -0.5; a < 1; a += 1)
    {
      landmarks.push_back(Point_d(std::vector<FT>({a,0,0,0})));
      landmarks.push_back(Point_d(std::vector<FT>({0,a,0,0})));
      landmarks.push_back(Point_d(std::vector<FT>({0,0,a,0})));
      landmarks.push_back(Point_d(std::vector<FT>({0,0,0,a})));
    }
  // 96 vertices
  for (FT a = -phi/2; a < phi; a += phi)
    for (FT b = -0.5; b < 1; b += 1)
      for (FT c = -phi_1/2; c < phi_1; c += phi_1)
        {
          landmarks.push_back(Point_d(std::vector<FT>({a,b,c,0})));
          landmarks.push_back(Point_d(std::vector<FT>({b,a,0,c})));
          landmarks.push_back(Point_d(std::vector<FT>({c,0,a,b})));
          landmarks.push_back(Point_d(std::vector<FT>({0,c,b,a})));
          landmarks.push_back(Point_d(std::vector<FT>({a,c,0,b})));
          landmarks.push_back(Point_d(std::vector<FT>({a,0,b,c})));
          landmarks.push_back(Point_d(std::vector<FT>({c,b,0,a})));
          landmarks.push_back(Point_d(std::vector<FT>({0,b,a,c})));
          landmarks.push_back(Point_d(std::vector<FT>({b,0,c,a})));
          landmarks.push_back(Point_d(std::vector<FT>({0,a,c,b})));
          landmarks.push_back(Point_d(std::vector<FT>({b,c,a,0})));
          landmarks.push_back(Point_d(std::vector<FT>({c,a,b,0})));
        }
  for (int i = 0; i < 120; ++i)
    landmarks_ind.push_back(i);
}

int landmark_perturbation(Point_Vector &W, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  //********************Preface: origin point
  clock_t start, end;
  int D = W[0].size();
  std::vector<FT> orig_vector;
  for (int i=0; i<D; i++)
    orig_vector.push_back(0);
  Point_d origin(orig_vector);
  //Distance dist;
  //dist.transformed_distance(0,1);
  //******************** Constructing a WL matrix
  int nbP = W.size();
  int nbL = landmarks.size();
  Euclidean_distance ed;
  FT lambda = ed.transformed_distance(landmarks[0],landmarks[1]);
  //std::cout << "Lambda=" << lambda << std::endl;
    //FT lambda = 0.1;//Euclidean_distance();
  STraits traits(&(landmarks[0]));
  std::vector< std::vector <int> > WL(nbP);
  Tree L(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nbL),
         typename Tree::Splitter(),
         traits);
  /*Tree2 L2(boost::counting_iterator<std::ptrdiff_t>(0),
           boost::counting_iterator<std::ptrdiff_t>(nbL),
           typename Tree::Splitter(),
           STraits(&(landmarks[0])));
  */
  std::cout << "Enter (D+1) nearest landmarks\n";
  //std::cout << "Size of the tree is " << L.size() << std::endl;
  start = clock();
  for (int i = 0; i < nbP; i++)
    {
      //std::cout << "Entered witness number " << i << std::endl;
      Point_d& w = W[i];
      //std::cout << "Safely constructed a point\n";
      ////Search D+1 nearest neighbours from the tree of landmarks L
      /*
      if (w[0]>0.95)
        std::cout << i << std::endl;
      */
      K_neighbor_search search(L, w, D, FT(0), true,
                               //CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks[0])) );
                               CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks[0])) );
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
  end = clock();
  std::cout << "Landmark choice for " << nbL << " landmarks took "
            << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  std::string out_file = "wl_result";
  write_wl(out_file,WL);
  
  //******************** Constructng a witness complex
  std::cout << "Entered witness complex construction\n";
  Witness_complex<> witnessComplex;
  witnessComplex.setNbL(nbL);
  start = clock();
  witnessComplex.witness_complex(WL);
  //
  end = clock();
  std::cout << "Howdy world! The process took "
       << (double)(end-start)/CLOCKS_PER_SEC << " s. \n";
  //witnessComplex.witness_complex(WL);
  /*
  if (witnessComplex.is_witness_complex(WL))
    std::cout << "!!YES. IT IS A WITNESS COMPLEX!!\n";
  else
    std::cout << "??NO. IT IS NOT A WITNESS COMPLEX??\n";   
  */ 
  //******************** Making a set of bad link landmarks
  std::cout << "Entered bad links\n";
  std::set< int > perturbL;
  int count_badlinks = 0;
  //std::cout << "Bad links around ";
  std::vector< int > count_bad(D);
  std::vector< int > count_good(D);
  for (auto u: witnessComplex.complex_vertex_range())
    if (!witnessComplex.has_good_link(u, count_bad, count_good))
      {
        //std::cout << "Landmark " << u << " start!" << std::endl;
        //perturbL.insert(u);
        count_badlinks++;
        //std::cout << u << " ";
        Point_d& l = landmarks[u];
        Fuzzy_sphere fs(l, sqrt(lambda), 0, traits);
        std::vector<int> curr_perturb;
        L.search(std::insert_iterator<std::vector<int>>(curr_perturb,curr_perturb.begin()),fs);
        for (int i: curr_perturb)
          perturbL.insert(i%nbL);
        //L.search(std::inserter(perturbL,perturbL.begin()),fs);
        //L.search(std::ostream_iterator<int>(std::cout,"\n"),fs);
        //std::cout << "PerturbL size is " << perturbL.size() << std::endl;
      }
  for (unsigned int i = 0; i != count_good.size(); i++)
    if (count_good[i] != 0)
      std::cout << "count_good[" << i << "] = " << count_good[i] << std::endl;
  for (unsigned int i = 0; i != count_bad.size(); i++)
    if (count_bad[i] != 0)
      std::cout << "count_bad[" << i << "] = " << count_bad[i] << std::endl;
  std::cout << "\nBad links total: " << count_badlinks << " Points to perturb: " << perturbL.size() << std::endl;
  //std::cout << "landmark[0][0] before" << landmarks[0][0] << std::endl;
  //*********************** Perturb bad link landmarks
  
  for (auto u: perturbL)
    {
      Random_point_iterator rp(D,sqrt(lambda)/8*nbL/count_badlinks);
      //std::cout << landmarks[u] << std::endl;
      
      std::vector<FT> point;
      for (int i = 0; i < D; i++)
        {
          while (K().squared_distance_d_object()(*rp,origin) < lambda/256)
            rp++;
          FT coord = W[landmarks_ind[u]][i] + (*rp)[i];
          //FT coord = landmarks[u][i] + (*rp)[i];
          if (coord > 1)
            point.push_back(coord-1);
          else if (coord < -1)
            point.push_back(coord+1);
          else
            point.push_back(coord);
        }
      landmarks[u] = Point_d(point);
      //std::cout << landmarks[u] << std::endl;
    }
  
  //std::cout << "landmark[0][0] after" << landmarks[0][0] << std::endl;
  std::cout << "lambda=" << lambda << std::endl;

  //std::cout << "WL size" << WL.size() << std::endl;
  /*
  std::cout << "L:" << std::endl;
  for (int i = 0; i < landmarks.size(); i++)
      std::cout << landmarks[i] << std::endl;
  */
  
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
  return count_badlinks;
}


int main (int argc, char * const argv[])
{
  
  if (argc != 4)
    {
      std::cerr << "Usage: " << argv[0]
                << " nbP nbL dim\n";
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
  //clock_t start, end;
  //Construct the Simplex Tree
  //Witness_complex<> witnessComplex;
 
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
  //  witnessComplex.witness_complex_from_points(point_vector);
  //int nbP = point_vector.size();
  //std::vector<std::vector< int > > WL(nbP);
  //std::set<int> L;
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  //Point_etiquette_map L_i;
  //start = clock();
  //witnessComplex.landmark_choice_by_furthest_points(point_vector, point_vector.size(), WL);
  bool ok=false;
  while (!ok)
    {
      ok = true;
      L = {};
      chosen_landmarks = {};
      landmark_choice(point_vector, nbP, nbL, L, chosen_landmarks);
      //landmark_choice_600cell(point_vector, nbP, nbL, L, chosen_landmarks);
      /*
      for (auto i: chosen_landmarks)
        {
          ok = ok && (std::count(chosen_landmarks.begin(),chosen_landmarks.end(),i) == 1);
          if (!ok) break;
        }
      */
    }
  int bl = nbL, curr_min = bl;
  //write_points("landmarks/initial_pointset",point_vector);
  //write_points("landmarks/initial_landmarks",L);
  
  for (int i = 0; bl > 0; i++)
    //for (int i = 0; i < 1; i++)
    {
      std::cout << "========== Start iteration " << i << "== curr_min(" << curr_min << ")========\n";
      bl=landmark_perturbation(point_vector, L, chosen_landmarks);
      if (bl < curr_min)
        curr_min=bl;
      //write_points("landmarks/landmarks0",L);
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
