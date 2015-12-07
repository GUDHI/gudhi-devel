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
#include <chrono>

#include <sys/types.h>
#include <sys/stat.h>
//#include <stdlib.h>

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

void generate_points_grid(Point_Vector& W, int width, int D)
{
  int nb_points = 1;
  for (int i = 0; i < D; ++i)
    nb_points *= width;
  for (int i = 0; i < nb_points; ++i)
    {
      std::vector<double> point;
      int cell_i = i;
      for (int l = 0; l < D; ++l)
        {
          point.push_back(0.01*(cell_i%width));
          cell_i /= width;
        }
      W.push_back(point);
    }
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

void insert_delaunay_landmark_with_copies(Point_Vector& W, int chosen_landmark, std::vector<int>& landmarks_ind, Delaunay_triangulation& delaunay, int& landmark_count)
{
  int D = W[0].size();
  int nb_cells = pow(3, D);
  for (int i = 0; i < nb_cells; ++i)
    {
      std::vector<FT> point;
      int cell_i = i;
      for (int l = 0; l < D; ++l)
        {
                point.push_back(W[chosen_landmark][l] + 2.0*(cell_i%3-1));
                cell_i /= 3;
        }
      delaunay.insert(point);
    }
  landmarks_ind.push_back(chosen_landmark);
  landmark_count++;
}




////////////////////////////////////////////////////////////////////////
// OLD CODE VVVVVVVV
////////////////////////////////////////////////////////////////////////


/*
bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta)
{
  Euclidean_distance ed;
  Delaunay_triangulation::Vertex_handle v;
  Delaunay_triangulation::Face f(t.current_dimension()); 
  Delaunay_triangulation::Facet ft; 
  Delaunay_triangulation::Full_cell_handle c; 
  Delaunay_triangulation::Locate_type lt;
  c = t.locate(p, lt, f, ft, v);
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    if (!t.is_infinite(fc_it))
      {
        std::vector<Point_d> vertices;
        for (auto v_it = fc_it->vertices_begin(); v_it != fc_it->vertices_end(); ++v_it)
          vertices.push_back((*v_it)->point());
        Sphere_d cs(D, vertices.begin(), vertices.end());
        Point_d center_cs = cs.center();
        FT r = sqrt(ed.transformed_distance(center_cs, fc_it->vertex(1)->point()));
        FT dist2 = ed.transformed_distance(center_cs, p);
        //if the new point is inside the protection ball of a non conflicting simplex
        if (dist2 >= r*r && dist2 <= (r+delta)*(r+delta)) 
          return true;
      }
  return false;
}

bool triangulation_is_protected(Delaunay_triangulation& t, FT delta)
{
  Euclidean_distance ed;
  int D = t.current_dimension();
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    if (!t.is_infinite(fc_it))
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        {
          //check if vertex belongs to the face
          bool belongs = false;
          for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
            if (v_it == *fc_v_it)
              {
                belongs = true;
                break;
              }
          if (!belongs)
            {
              std::vector<Point_d> vertices;
              for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
                vertices.push_back((*fc_v_it)->point());
              Sphere_d cs(D, vertices.begin(), vertices.end());
              Point_d center_cs = cs.center();
              FT r = sqrt(ed.transformed_distance(center_cs, fc_it->vertex(1)->point()));
              FT dist2 = ed.transformed_distance(center_cs, v_it->point());
              //if the new point is inside the protection ball of a non conflicting simplex
              if (dist2 <= (r+delta)*(r+delta)) 
                return false;
            }
        }
  return true;
}

void fill_landmark_copies(Point_Vector& W, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  int D = W[0].size();
  int nb_cells = pow(3, D);
  int nbL = landmarks_ind.size();
  // Fill landmarks
  for (int i = 0; i < nb_cells-1; ++i)
    for (int j = 0; j < nbL; ++j)
      {
        int cell_i = i;
        Point_d point;
        for (int l = 0; l < D; ++l)
          {
            point.push_back(W[landmarks_ind[j]][l] + 2.0*(cell_i-1));
            cell_i /= 3;
          }
        landmarks.push_back(point);
      }
}

void landmark_choice_by_delaunay(Point_Vector& W, int nbP, int nbL, Point_Vector& landmarks, std::vector<int>& landmarks_ind, FT delta)
{
  int D = W[0].size();
  Delaunay_triangulation t(D);
  CGAL::Random rand;
  int chosen_landmark;
  int landmark_count = 0;
  for (int i = 0; i <= D+1; ++i)
    {
      do chosen_landmark = rand.get_int(0,nbP);
      while (std::count(landmarks_ind.begin(),landmarks_ind.end(),chosen_landmark)!=0);
      insert_delaunay_landmark_with_copies(W, chosen_landmark, landmarks_ind, t, landmark_count);
    }
  while (landmark_count < nbL)
    {      
      do chosen_landmark = rand.get_int(0,nbP);
      while (std::count(landmarks_ind.begin(),landmarks_ind.end(),chosen_landmark)!=0);      
      // If no conflicts then insert in every copy of T^3
      if (!is_violating_protection(W[chosen_landmark], t, D, delta))
        insert_delaunay_landmark_with_copies(W, chosen_landmark, landmarks_ind, t, landmark_count);
    }
  fill_landmark_copies(W, landmarks, landmarks_ind);
}


void landmark_choice_protected_delaunay(Point_Vector& W, int nbP, Point_Vector& landmarks, std::vector<int>& landmarks_ind, FT delta)
{
  int D = W[0].size();
  Torus_distance td;
  Euclidean_distance ed;
  Delaunay_triangulation t(D);
  CGAL::Random rand;
  int landmark_count = 0;
  std::list<int> index_list;
  //      shuffle the list of indexes (via a vector)
  {
    std::vector<int> temp_vector;
    for (int i = 0; i < nbP; ++i)
      temp_vector.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(temp_vector.begin(), temp_vector.end(), std::default_random_engine(seed));
    for (std::vector<int>::iterator it = temp_vector.begin(); it != temp_vector.end(); ++it)
      index_list.push_front(*it);
  }
  //      add the first D+1 vertices to form one non-empty cell 
  for (int i = 0; i <= D+1; ++i)
    {
      insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count);
      index_list.pop_front();
    }
  //      add other vertices if they don't violate protection
  std::list<int>::iterator list_it = index_list.begin();
  while (list_it != index_list.end())
    if (!is_violating_protection(W[*list_it], t, D, delta))
      {
          // If no conflicts then insert in every copy of T^3      
          insert_delaunay_landmark_with_copies(W, *list_it, landmarks_ind, t, landmark_count);
          index_list.erase(list_it);
          list_it = index_list.begin();
        }
    else
      list_it++;
  fill_landmark_copies(W, landmarks, landmarks_ind);
}


int landmark_perturbation(Point_Vector &W, int nbL, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  //******************** Preface: origin point
  int D = W[0].size();
  std::vector<FT> orig_vector;
  for (int i=0; i<D; i++)
    orig_vector.push_back(0);
  Point_d origin(orig_vector);

  //******************** Constructing a WL matrix
  int nbP = W.size();
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
       if (!witnessComplex.has_good_link(u, count_bad, count_good))
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


int main (int argc, char * const argv[])
{
  if (argc != 5)
    {
      std::cerr << "Usage: " << argv[0]
                << " nbP nbL dim delta\n";
      return 0;
    }
  int nbP       = atoi(argv[1]);
  int nbL       = atoi(argv[2]);
  int dim       = atoi(argv[3]);
  FT delta      = atof(argv[4]);
  
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  generate_points_random_box(point_vector, nbP, dim);
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  bool ok=false;
  while (!ok)
    {
      ok = true;
      L = {};
      chosen_landmarks = {};
      //landmark_choice_by_delaunay(point_vector, nbP, nbL, L, chosen_landmarks, delta);
      landmark_choice_protected_delaunay(point_vector, nbP, L, chosen_landmarks, delta);
      nbL = chosen_landmarks.size();
      std::cout << "Number of landmarks is " << nbL << std::endl;
      //int width = (int)pow(nbL, 1.0/dim); landmark_choice_bcc(point_vector, nbP, width, L, chosen_landmarks);
      for (auto i: chosen_landmarks)
        {
          ok = ok && (std::count(chosen_landmarks.begin(),chosen_landmarks.end(),i) == 1);
          if (!ok) break;
        }
      
    }
  int bl = nbL, curr_min = bl;
  write_points("landmarks/initial_pointset",point_vector);
  //write_points("landmarks/initial_landmarks",L);
  //for (int i = 0; i < 1; i++)
  for (int i = 0; bl > 0; i++)
    {
      std::cout << "========== Start iteration " << i << "== curr_min(" << curr_min << ")========\n";
      bl=landmark_perturbation(point_vector, nbL, L, chosen_landmarks);
      if (bl < curr_min)
        curr_min=bl;
      write_points("landmarks/landmarks0",L);
    }
  
}
*/
