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

// Avoiding the max arity issue with CGAL
#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif

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
#include "generators.h"
#include "output.h"
//#include "protected_sets/protected_sets.h"
#include "protected_sets/protected_sets_paper2.h"

#include <CGAL/Cartesian_d.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/enum.h>

#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/Delaunay_triangulation.h>


#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

using namespace Gudhi;
//using namespace boost::filesystem;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef K::Vector_d Vector_d;
typedef K::Oriented_side_d Oriented_side_d;
typedef K::Has_on_positive_side_d Has_on_positive_side_d;

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

typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;
typedef Delaunay_triangulation::Facet Facet;
typedef Delaunay_triangulation::Vertex_handle Delaunay_vertex;
typedef Delaunay_triangulation::Full_cell_handle Full_cell_handle;
//typedef CGAL::Sphere_d<K> Sphere_d;
typedef K::Sphere_d Sphere_d;
typedef K::Hyperplane_d Hyperplane_d;

/*//////////////////////////////////////
 * GLOBAL VARIABLES ********************
 *//////////////////////////////////////

//NA bool toric=false;
bool power_protection = true;
bool grid_points = true;
bool is2d = true;
//FT  _sfty = pow(10,-14);
bool torus = false;


bool triangulation_is_protected(Delaunay_triangulation& t, FT delta)
{
  std::cout << "Start protection verification\n";
  Euclidean_distance ed;
  // Fill the map Vertices -> Numbers
  std::map<Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
  int ind = 0;
  for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
    {
      if (t.is_infinite(v_it))
        continue;
      index_of_vertex[v_it] = ind++;
    }
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    if (!t.is_infinite(fc_it))
      {
        std::vector<Point_d> vertices;
        for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
          vertices.push_back((*fc_v_it)->point());                
        Sphere_d cs( vertices.begin(), vertices.end());
        Point_d center_cs = cs.center();
        FT r = sqrt(ed.transformed_distance(center_cs, fc_it->vertex(0)->point()));
        for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
          if (!t.is_infinite(v_it))
          //check if vertex belongs to the face
            if (!vertex_is_in_full_cell(v_it, fc_it))
              {
                FT dist2 = ed.transformed_distance(center_cs, v_it->point());
                //if the new point is inside the protection ball of a non conflicting simplex
                //std::cout << "Dist^2 = " << dist2 << " (r+delta)*(r+delta) = " << (r+delta)*(r+delta) << " r^2 = " << r*r <<"\n";
                if (!power_protection)
                  if (dist2 <= (r+delta)*(r+delta) && dist2 >= r*r)
                    {
                      write_delaunay_mesh(t, v_it->point(), is2d);
                      // Output the problems
                      std::cout << "Problematic vertex " << index_of_vertex[v_it] << " ";
                      std::cout << "Problematic cell ";
                      for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
                        if (!t.is_infinite(*vh_it))
                          std::cout << index_of_vertex[*vh_it] << " ";
                      std::cout << "\n";
                      std::cout << "r^2 = " << r*r << ", d^2 = " << dist2 << ", (r+delta)^2 = " << (r+delta)*(r+delta) << "\n";
                      return false;
                    }
                if (power_protection)
                  if (dist2 <= r*r+delta*delta && dist2 >= r*r)
                    {
                      write_delaunay_mesh(t, v_it->point(), is2d);
                      std::cout << "Problematic vertex " << *v_it << " ";
                      std::cout << "Problematic cell " << *fc_it << "\n";
                      std::cout << "r^2 = " << r*r << ", d^2 = " << dist2 << ", r^2+delta^2 = " << r*r+delta*delta << "\n";
                      return false;
                    }
              }
      } 
  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// SAMPLING RADIUS 
//////////////////////////////////////////////////////////////////////////////////////////////////////////

FT sampling_radius(Delaunay_triangulation& t, FT epsilon0)
{
  FT epsilon2 = 0;
  Point_d control_point;
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    {
      if (t.is_infinite(fc_it))
        continue;
      Point_Vector vertices;
      for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
        vertices.push_back((*fc_v_it)->point());
      Sphere_d cs( vertices.begin(), vertices.end());
      Point_d csc = cs.center();
      bool in_cube = true; 
      for (auto xi = csc.cartesian_begin(); xi != csc.cartesian_end(); ++xi)
        if (*xi > 1.0 || *xi < -1.0)
          {
            in_cube = false; break;
          }
      if (!in_cube)
        continue;
      FT r2 = Euclidean_distance().transformed_distance(cs.center(), *(vertices.begin()));
      if (epsilon2 < r2)
        {
          epsilon2 = r2;
          control_point = (*vertices.begin());
        }
    }
  if (epsilon2 < epsilon0*epsilon0)
    {
      std::cout << "ACHTUNG! E' < E\n";
      std::cout << "eps = " << epsilon0 << " eps' = " << sqrt(epsilon2) << "\n";
      write_delaunay_mesh(t, control_point, is2d);
    }
  return sqrt(epsilon2);
}

FT point_sampling_radius_by_delaunay(Point_Vector& points, FT epsilon0)
{
  Delaunay_triangulation t(points[0].size());
  t.insert(points.begin(), points.end());
  return sampling_radius(t, epsilon0);
}

// A little script to make a tikz histogram of epsilon distribution
// Returns the average epsilon
FT epsilon_histogram(Delaunay_triangulation& t, int n)
{
  FT epsilon_max = 0; //sampling_radius(t,0);
  FT sum_epsilon = 0;
  int count_simplices = 0;
  std::vector<int> histo(n+1, 0);
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    {
      if (t.is_infinite(fc_it))
        continue;
      Point_Vector vertices;
      for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
        vertices.push_back((*fc_v_it)->point());
      Sphere_d cs( vertices.begin(), vertices.end());
      Point_d csc = cs.center();
      bool in_cube = true; 
      for (auto xi = csc.cartesian_begin(); xi != csc.cartesian_end(); ++xi)
        if (*xi > 1.0 || *xi < -1.0)
          {
            in_cube = false; break;
          }
      if (!in_cube)
        continue;
      FT r = sqrt(Euclidean_distance().transformed_distance(cs.center(), *(vertices.begin())));
      if (r > epsilon_max)
        epsilon_max = r;
      sum_epsilon += r;
      count_simplices++;
      histo[floor(r/epsilon_max*n)]++;
    }
  std::ofstream ofs ("histogram.tikz", std::ofstream::out);
  FT barwidth = 20.0/n;
  int max_value = *(std::max_element(histo.begin(), histo.end()));
  std::cout << max_value << std::endl;
  FT ten_power = pow(10, ceil(log10(max_value)));
  FT max_histo = ten_power;
  if (max_value/ten_power < 2)
    max_histo = 0.2*ten_power;
  if (max_value/ten_power < 5)
    max_histo = 0.5*ten_power;  
  std::cout << ceil(log10(max_value)) << std::endl << max_histo << std::endl;
  FT unitht = max_histo/10.0;

  ofs << "\\draw[->] (0,0) -- (0,11);\n" <<
    "\\draw[->] (0,0) -- (21,0);\n" <<
    "\\foreach \\i in {1,...,10}\n" <<
    "\\draw (0,\\i) -- (-0.1,\\i);\n" <<
    "\\foreach \\i in {1,...,20}\n" <<
    "\\draw (\\i,0) -- (\\i,-0.1);\n" <<
 
    "\\node at (-1,11) {$\\epsilon$};\n" << 
    "\\node at (22,-1) {$\\epsilon/\\epsilon_{max}$};\n" << 
    "\\node at (-0.5,-0.5) {0};\n" << 
    "\\node at (-0.5,10) {" << max_histo << "};\n" << 
    "\\node at (20,-0.5) {1};\n";
    

  for (int i = 0; i < n; ++i)
    ofs << "\\draw (" << barwidth*i << "," << histo[i]/unitht << ") -- ("
        << barwidth*(i+1) << "," << histo[i]/unitht << ") -- ("
  << barwidth*(i+1) << ",0) -- (" << barwidth*i << ",0) -- cycle;\n";
  
  ofs.close();

  //return sum_epsilon/count_simplices;
  return epsilon_max;
}

FT epsilon_histogram_by_delaunay(Point_Vector& points, int n)
{
  Delaunay_triangulation t(points[0].size());
  t.insert(points.begin(), points.end());
  return epsilon_histogram(t, n);
}


int landmark_perturbation(Point_Vector &W, int nbL, Point_Vector& landmarks, std::vector<int>& landmarks_ind, std::vector<std::vector<int>>& full_cells)
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
  //write_wl(out_file,WL);
  
  //******************** Constructng a witness complex
  std::cout << "Entered witness complex construction\n";
  Witness_complex<> witnessComplex;
  witnessComplex.setNbL(nbL);
  witnessComplex.witness_complex(WL);

  //******************** Verifying if all full cells are in the complex

  int in=0, not_in=0;
  for (auto cell : full_cells)
    {
      //print_vector(cell);
      if (witnessComplex.find(cell) != witnessComplex.null_simplex())
        in++;
      else
        not_in++;
    }
  std::cout << "Out of all the cells in Delaunay triangulation:\n" << in << " are in the witness complex\n" <<
    not_in << " are not.\n";
  
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
  /*
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
  
  //write_edges("landmarks/edges", witnessComplex, landmarks);
  /*
  return count_badlinks;
  */
  return 0;
}

int main (int argc, char * const argv[])
{
  power_protection = true;//false;
  grid_points = false;//true;
  torus = true;
  
  if (argc != 4)
    {
      std::cerr << "Usage: " << argv[0]
                << " nbP dim delta\n";
      return 0;
    }
  int nbP       = atoi(argv[1]);
  int dim       = atoi(argv[2]);
  double theta0 = atof(argv[3]);
  //double delta  = atof(argv[3]);

  is2d = (dim == 2);
  
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  if (grid_points)
    {
      generate_points_grid(point_vector, (int)pow(nbP, 1.0/dim), dim, torus);
      nbP = (int)pow((int)pow(nbP, 1.0/dim), dim);
    }
  else
    generate_points_random_box(point_vector, nbP, dim);
  FT epsilon = point_sampling_radius_by_delaunay(point_vector, 0);
  //FT epsilon = epsilon_histogram_by_delaunay(point_vector,50);
  std::cout << "Initial epsilon = " << epsilon << std::endl;
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  //write_points("landmarks/initial_pointset",point_vector);
  //write_points("landmarks/initial_landmarks",L);
  CGAL::Timer timer;

  int n = 1;
  std::vector<FT> values(n,0);
  std::vector<FT> time(n,0);

  //FT step = 0.001;
  //FT delta = 0.01*epsilon;
  //FT alpha = 0.5;
  //FT step = atof(argv[3]);

  start_experiments(point_vector, theta0, chosen_landmarks, epsilon);
    
  // for (int i = 0; i < n; i++)
  // //for (int i = 0; bl > 0; i++)
  //   {
  //     //std::cout << "========== Start iteration " << i << "== curr_min(" << curr_min << ")========\n";
  //     //double delta = pow(10, -(1.0*i)/2);
  //     //delta = step*i*epsilon;
  //     //theta0 = step*i;
  //     std::cout << "delta/epsilon = " << delta/epsilon << std::endl;
  //     std::cout << "theta0 = " << theta0 << std::endl;
  //     // Averaging the result
  //     int sum_values = 0;
  //     int nb_iterations = 1;
  //     std::vector<std::vector<int>> full_cells;
  //     for (int i = 0; i < nb_iterations; ++i)
  //       {
  //         //L = {};
  //         chosen_landmarks = {};
  //         //full_cells = {};
  //         //timer.start();
  //         //protected_delaunay(point_vector, nbP, L, chosen_landmarks, delta, epsilon, alpha, theta0, full_cells, torus, power_protection);
  //         protected_delaunay(point_vector, chosen_landmarks, delta, epsilon, alpha, theta0, torus, power_protection);
  //         //timer.stop();
  //         sum_values += chosen_landmarks.size();
  //       }
  //     //FT epsilon2 = point_sampling_radius_by_delaunay(L, epsilon);
  //     //std::cout << "Final epsilon = " << epsilon2  << ". Ratio = " << epsilon2/epsilon << std::endl;
  //     //write_points("landmarks/initial_landmarks",L);
  //     //std::cout << "delta/epsilon' = " << delta/epsilon2 << std::endl;
  //     FT nbL = (sum_values*1.0)/nb_iterations;
  //     //values[i] = pow((1.0*nbL)/nbP, -1.0/dim);
  //     values[i] = (1.0*nbL)/nbP;
  //     std::cout << "Number of landmarks = " << nbL << ", time= " << timer.time() << "s"<< std::endl;
  //     //landmark_perturbation(point_vector, nbL, L, chosen_landmarks, full_cells);
  //     time[i] = timer.time();
  //     timer.reset();
  //     //write_points("landmarks/landmarks0",L);
  //   }
  
  // // OUTPUT A PLOT
  // FT hstep = 20.0/(n-1);
  // FT wstep = 10.0;

  // std::ofstream ofs("N'Nplot.tikz", std::ofstream::out);
  // ofs << "\\draw[red] (0," << wstep*values[0] << ")";
  // for (int i = 1; i < n; ++i)
  //   ofs << " -- (" << hstep*i << "," << wstep*values[i] << ")";
  // ofs << ";\n";
  // ofs.close();
  /*
  wstep = 0.1;
  ofs = std::ofstream("time.tikz", std::ofstream::out);
  ofs << "\\draw[red] (0," << wstep*time[0] << ")";
  for (int i = 1; i < n; ++i)
    ofs << " -- (" << hstep*i << "," << wstep*time[i] << ")";
  ofs << ";\n";
  ofs.close();
  
  
  std::vector<std::vector<int>> full_cells;
  timer.start();
  landmark_choice_protected_delaunay(point_vector, nbP, L, chosen_landmarks, delta, full_cells);
  timer.stop();
  FT epsilon2 = point_sampling_radius_by_delaunay(L);
  std::cout << "Final epsilon = " << epsilon2  << ". Ratio = " << epsilon/epsilon2 << std::endl;
  write_points("landmarks/initial_landmarks",L);
  int nbL = chosen_landmarks.size();
  std::cout << "Number of landmarks = " << nbL << ", time= " << timer.time() << "s"<< std::endl;
  //landmark_perturbation(point_vector, nbL, L, chosen_landmarks, full_cells);
  timer.reset();
  */
}
