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
#include <unistd.h>

//#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/Witness_complex.h"
#include "gudhi/reader_utils.h"
#include "Torus_distance.h" 

#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Kernel_d/Sphere_d.h>

#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Random.h>
#include <CGAL/Delaunay_triangulation.h>


#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

using namespace Gudhi;
//using namespace boost::filesystem;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
//typedef CGAL::Cartesian_d<double> K;
//typedef CGAL::Point_d<K> Point_d;
typedef K::FT FT;
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
//typedef CGAL::Random_points_in_cube_d<CGAL::Point_d<CGAL::Cartesian_d<FT> > > Random_cube_iterator;
typedef CGAL::Random_points_in_cube_d<Point_d> Random_cube_iterator;
typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;

typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;
typedef Delaunay_triangulation::Facet Facet;
typedef CGAL::Sphere_d<K> Sphere_d;

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

void generate_points_random_box(Point_Vector& W, int nbP, int dim)
{
  /*
  Random_cube_iterator rp(dim, 1);
  for (int i = 0; i < nbP; i++)
    {
      std::vector<double> point;
      for (auto it = rp->cartesian_begin(); it != rp->cartesian_end(); ++it)
        point.push_back(*it);
      W.push_back(Point_d(point));
      rp++;
    }
  */
  Random_cube_iterator rp(dim, 1.0);
  for (int i = 0; i < nbP; i++)
    {
      W.push_back(*rp++);
    }
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


void write_points( std::string file_name, std::vector< Point_d > & points)
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

void write_edges(std::string file_name, Witness_complex<>& witness_complex, Point_Vector& landmarks)
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
  int chosen_landmark;
  Point_d* p;
  CGAL::Random rand;
  for (int i = 0; i < nbL; i++)
    {
      //      while (!res.second)
      //  {
      do chosen_landmark = rand.get_int(0,nbP);
      while (std::count(landmarks_ind.begin(),landmarks_ind.end(),chosen_landmark)!=0);
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

void aux_fill_grid(Point_Vector& W, int& width, Point_Vector& landmarks, std::vector<int>& landmarks_ind, std::vector<bool> & curr_pattern)
{
  int D = W[0].size();
  int nb_points = 1;
  for (int i = 0; i < D; ++i)
    nb_points *= width;
  for (int i = 0; i < nb_points; ++i)
    {
      std::vector<double> point;
      int cell_i = i;
      for (int l = 0; l < D; ++l)
        {
          if (curr_pattern[l])
            point.push_back(-1.0+(2.0/width)*(cell_i%width)+(1.0/width));
          else
            point.push_back(-1.0+(2.0/width)*(cell_i%width));
          cell_i /= width;
        }
      landmarks.push_back(Point_d(point));
      landmarks_ind.push_back(0);//landmarks_ind.push_back(W.size());
      //std::cout << "Added point " << W.size() << std::endl;;
      //W.push_back(Point_d(point));
    }
}
  
void aux_put_halves(Point_Vector& W, int& width, Point_Vector& landmarks, std::vector<int>& landmarks_ind, std::vector<bool>& curr_pattern, std::vector<bool>::iterator curr_pattern_it, std::vector<bool>::iterator bool_it, std::vector<bool>::iterator bool_end)
{
  if (curr_pattern_it != curr_pattern.end())
    {
      if (bool_it != bool_end)
        {
          *curr_pattern_it = false;
          aux_put_halves(W, width, landmarks, landmarks_ind, curr_pattern, curr_pattern_it+1, bool_it, bool_end);
          *curr_pattern_it = true;
          aux_put_halves(W, width, landmarks, landmarks_ind, curr_pattern, curr_pattern_it+1, bool_it+1, bool_end);
        }
    }
  else
    if (*bool_it)
      {
        std::cout << "Filling the pattern ";
        for (bool b: curr_pattern)
          if (b) std::cout << '1';
          else   std::cout << '0';
        std::cout << "\n";
        aux_fill_grid(W, width, landmarks, landmarks_ind, curr_pattern);
      }  
}

void landmark_choice_cs(Point_Vector& W, int width, Point_Vector& landmarks, std::vector<int>& landmarks_ind, std::vector<bool>& face_centers)
{
  std::cout << "Enter landmark choice to kd tree\n";
  //int chosen_landmark;
  CGAL::Random rand;
  //To speed things up check the last true in the code and put it as the finishing condition
  unsigned last_true = face_centers.size()-1;
  while (!face_centers[last_true] && last_true != 0)
    last_true--;
  //Recursive procedure to understand where we put +1/2 in centers' coordinates
  std::vector<bool> curr_pattern(W[0].size(), false);
  aux_put_halves(W, width, landmarks, landmarks_ind, curr_pattern, curr_pattern.begin(), face_centers.begin(), face_centers.begin()+(last_true+1));
  std::cout << "The number of landmarks is: " << landmarks.size() << std::endl;

 }

int landmark_perturbation(Point_Vector &W, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  //******************** Preface: origin point
  int D = W[0].size();
  std::vector<FT> orig_vector;
  for (int i=0; i<D; i++)
    orig_vector.push_back(0);
  Point_d origin(orig_vector);

  //******************** Constructing a WL matrix
  int nbP = W.size();
  int nbL = landmarks.size();
  Euclidean_distance ed;
  FT lambda = ed.transformed_distance(landmarks[0],landmarks[1]);
  std::vector<Point_d> landmarks_ext;
  int nb_cells = 1;
  for (int i = 0; i < D; ++i)
    nb_cells *= 3;
  for (int i = 0; i < nb_cells; ++i)
      for (int k = 0; k < nbL; ++k)
        {
          std::vector<double> point;
          int cell_i = i;
          for (int l = 0; l < D; ++l)
            {
              point.push_back(landmarks[k][l] + 2.0*((cell_i%3)-1.0));
              cell_i /= 3;
            }
          landmarks_ext.push_back(point);
        }
  write_points("landmarks/initial_landmarks",landmarks_ext);
  STraits traits(&(landmarks_ext[0]));
  std::vector< std::vector <int> > WL(nbP);

  //********************** Neighbor search in a Kd tree
  Tree L(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nb_cells*nbL),
         typename Tree::Splitter(),
         traits);
  std::cout << "Enter (D+1) nearest landmarks\n";
  for (int i = 0; i < nbP; i++)
    {
      Point_d& w = W[i];
      ////Search D+1 nearest neighbours from the tree of landmarks L
      K_neighbor_search search(L, w, D+1, FT(0), true,
                               CGAL::Distance_adapter<std::ptrdiff_t,Point_d*,Euclidean_distance>(&(landmarks_ext[0])) );
      for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
        {
          if (std::find(WL[i].begin(), WL[i].end(), (it->first)%nbL) == WL[i].end())
            WL[i].push_back((it->first)%nbL);
        }
      if (i == landmarks_ind[WL[i][0]])
        {
          FT dist = ed.transformed_distance(W[i], landmarks[WL[i][1]]);
          if (dist < lambda)
            lambda = dist;
        }
    }
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
  //std::cout << "Bad links around ";
  std::vector< int > count_bad(D);
  std::vector< int > count_good(D);
  for (auto u: witnessComplex.complex_vertex_range())
    {
      if (!witnessComplex.has_good_link(u, count_bad, count_good, D))
         {
           count_badlinks++;
           Point_d& l = landmarks[u];
           Fuzzy_sphere fs(l, sqrt(lambda)*3, 0, traits);
           std::vector<int> curr_perturb;
           L.search(std::insert_iterator<std::vector<int>>(curr_perturb,curr_perturb.begin()),fs);
           for (int i: curr_perturb)
             perturbL.insert(i%nbL);
       }
    }
  for (unsigned int i = 0; i != count_good.size(); i++)
    if (count_good[i] != 0)
      std::cout << "count_good[" << i << "] = " << count_good[i] << std::endl;
  for (unsigned int i = 0; i != count_bad.size(); i++)
    if (count_bad[i] != 0)
      std::cout << "count_bad[" << i << "] = " << count_bad[i] << std::endl;
  std::cout << "\nBad links total: " << count_badlinks << " Points to perturb: " << perturbL.size() << std::endl;
  
  //*********************** Perturb bad link landmarks  
  for (auto u: perturbL)
    {
      Random_point_iterator rp(D,sqrt(lambda)/8);
      std::vector<FT> point;
      for (int i = 0; i < D; i++)
        {
          while (K().squared_distance_d_object()(*rp,origin) < lambda/256)
            rp++;
          FT coord = landmarks[u][i] + (*rp)[i];
          if (coord > 1)
            point.push_back(coord-1);
          else if (coord < -1)
            point.push_back(coord+1);
          else
            point.push_back(coord);
        }
      landmarks[u] = Point_d(point);
    }
  std::cout << "lambda=" << lambda << std::endl;
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
  return count_badlinks;
}

void exaustive_search(Point_Vector& W, int width)
{
  int D = W[0].size()+1;
  int nb_points = pow(2,D);
  std::vector<bool> face_centers(D, false);
  int bl = 0; //Bad links
  std::vector<std::vector<bool>> good_patterns;
  for (int i = 0; i < nb_points; ++i)
    {
      int cell_i = i;
      for (int l = 0; l < D; ++l)
        {
          if (cell_i%2 == 0)
            face_centers[l] = false;
          else
            face_centers[l] = true;
          cell_i /= 2;
        }
      std::cout << "**Current pattern ";
      for (bool b: face_centers)
        if (b) std::cout << '1';
        else   std::cout << '0';
      std::cout << "\n";
      Point_Vector landmarks;
      std::vector<int> landmarks_ind;
      Point_Vector W_copy(W);
      landmark_choice_cs(W_copy, width, landmarks, landmarks_ind, face_centers);
      if (landmarks.size() != 0)
        {
          bl = landmark_perturbation(W_copy, landmarks, landmarks_ind);
          if ((1.0*bl)/landmarks.size() < 0.5)
            good_patterns.push_back(face_centers);
        }
    }
  std::cout << "The following patterns worked: ";
  for (std::vector<bool> pattern : good_patterns)
    {
      std::cout << "[";
      for (bool b: pattern)
        if (b) std::cout << '1';
        else   std::cout << '0';
      std::cout << "] ";
    }
  std::cout << "\n";        
}

int main (int argc, char * const argv[])
{
  unsigned nbP       = atoi(argv[1]);
  unsigned width     = atoi(argv[2]);
  unsigned dim       = atoi(argv[3]);
  std::string code   = (std::string) argv[4];
  bool e_option = false;
  int c;
  if (argc != 5)
    {
      std::cerr << "Usage: " << argv[0]
                << "witness_complex_cubic_systems nbP width dim code || witness_complex_systems -e nbP width dim\n"
                << "where nbP stands for the number of witnesses, width for the width of the grid, dim for dimension "
                << "and code is a sequence of (dim+1) symbols 0 and 1 representing if we take the centers of k-dimensional faces of the cubic system depending if it is 0 or 1."
                << "-e stands for the 'exaustive' option";
      return 0;
    }
  while ((c = getopt (argc, argv, "e::")) != -1)
    switch(c)
      {
      case 'e' :
        e_option  = true;
        nbP       = atoi(argv[2]);
        width     = atoi(argv[3]);
        dim       = atoi(argv[4]);
        break;
      default  :
        nbP       = atoi(argv[1]);
        width     = atoi(argv[2]);
        dim       = atoi(argv[3]);
        code   = (std::string) argv[4];
      }
  Point_Vector point_vector;
  generate_points_random_box(point_vector, nbP, dim);
  
  // Exaustive search
  if (e_option)
    {
      std::cout << "Start exaustive search!\n";
      exaustive_search(point_vector, width);
      return 0;
    }
  // Search with a specific cubic system
  std::vector<bool> face_centers;
  if (code.size() != dim+1)
    {
      std::cerr << "The code should contain (dim+1) symbols";
      return 1;
    }
  for (char c: code)
    if (c == '0')
      face_centers.push_back(false);
    else
      face_centers.push_back(true);
  std::cout << "Let the carnage begin!\n";
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  
  landmark_choice_cs(point_vector, width, L, chosen_landmarks, face_centers);
  
  int nbL = width; //!!!!!!!!!!!!!
  int bl = nbL, curr_min = bl;
  write_points("landmarks/initial_pointset",point_vector);
  //write_points("landmarks/initial_landmarks",L);
  //for (int i = 0; i < 1; i++)
  for (int i = 0; bl > 0; i++)
    {
      std::cout << "========== Start iteration " << i << "== curr_min(" << curr_min << ")========\n";
      bl=landmark_perturbation(point_vector, L, chosen_landmarks);
      if (bl < curr_min)
        curr_min=bl;
      write_points("landmarks/landmarks0",L);
    }
  
}
