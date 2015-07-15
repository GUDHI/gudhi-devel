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
typedef CGAL::Random_points_in_cube_d<Point_d> Random_cube_iterator;
typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;

typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;
typedef Delaunay_triangulation::Facet Facet;
typedef Delaunay_triangulation::Vertex_handle Delaunay_vertex;
typedef Delaunay_triangulation::Full_cell_handle Full_cell_handle;
//typedef CGAL::Sphere_d<K> Sphere_d;
typedef K::Sphere_d Sphere_d;
typedef K::Hyperplane_d Hyperplane_d;


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


void insert_delaunay_landmark_with_copies(Point_Vector& W, int chosen_landmark, std::vector<int>& landmarks_ind, Delaunay_triangulation& delaunay, int& landmark_count)
{
  delaunay.insert(W[chosen_landmark]);
  landmarks_ind.push_back(chosen_landmark);
  landmark_count++;
}

bool vertex_is_in_full_cell(Delaunay_triangulation::Vertex_handle v, Full_cell_handle fc)
{
  for (auto v_it = fc->vertices_begin(); v_it != fc->vertices_end(); ++v_it)
    if (*v_it == v)
      return true;
  return false;
}

bool new_cell_is_violated(Delaunay_triangulation& t, std::vector<Point_d>& vertices, bool is_infinite, const Point_d& p, FT delta)
{
  if (!is_infinite)
  // FINITE CASE
    {
      Sphere_d cs(vertices.begin(), vertices.end());
      Point_d center_cs = cs.center();
      FT r = sqrt(Euclidean_distance().transformed_distance(center_cs, vertices[0]));
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it)) 
          {
            //CGAL::Oriented_side side = Oriented_side_d()(cs, (v_it)->point()); 
            if (std::find(vertices.begin(), vertices.end(), v_it->point()) == vertices.end())
              {
                FT dist2 = Euclidean_distance().transformed_distance(center_cs, (v_it)->point());
                //if (dist2 >= r*r && dist2 <= (r+delta)*(r+delta))
                if (dist2 >= r*r && dist2 <= r*r+delta*delta) 
                  return true;                 
              }
          }
    }
  else
  // INFINITE CASE
    {
      Delaunay_triangulation::Vertex_iterator v = t.vertices_begin();
      while (t.is_infinite(v) || std::find(vertices.begin(), vertices.end(), v->point()) == vertices.end())
        v++;
      Hyperplane_d facet_plane(vertices.begin(), vertices.end(), v->point(), CGAL::ON_POSITIVE_SIDE);
      Vector_d orth_v = facet_plane.orthogonal_vector();
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it))
          if (std::find(vertices.begin(), vertices.end(), v_it->point()) == vertices.end())
            {
              std::vector<FT> coords;
              Point_d p = v_it->point();
              auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
              for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
                coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
              Point_d p_delta = Point_d(coords);
              bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
              bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
              if (!p_is_inside && p_delta_is_inside)
                return true;
            }
    }
  return false;
}
  

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, Full_cell_handle c, Full_cell_handle parent_cell, int index, int D, FT delta, std::vector<Full_cell_handle>& marked_cells)
{
  Euclidean_distance ed;
  std::vector<Point_d> vertices;
  if (!t.is_infinite(c))
    {
      // if the cell is finite, we look if the protection is violated
      for (auto v_it = c->vertices_begin(); v_it != c->vertices_end(); ++v_it)
        vertices.push_back((*v_it)->point());
      Sphere_d cs( vertices.begin(), vertices.end());
      Point_d center_cs = cs.center();
      FT r = sqrt(ed.transformed_distance(center_cs, vertices[0]));
      FT dist2 = ed.transformed_distance(center_cs, p);
      // if the new point is inside the protection ball of a non conflicting simplex
      //if (dist2 >= r*r && dist2 <= (r+delta)*(r+delta))
      if (dist2 >= r*r && dist2 <= r*r+delta*delta) 
        return true;    
      c->tds_data().mark_visited();
      marked_cells.push_back(c);
      // if the new point is inside the circumscribing ball : continue violation searching on neughbours
      if (dist2 < r*r)
        for (int i = 0; i < D+1; ++i)
          {
            Full_cell_handle next_c = c->neighbor(i);
            if (next_c->tds_data().is_clear() &&
                is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells))
              return true;
          }
      // if the new point is outside the protection sphere
      else
        {
          // facet f is on the border of the conflict zone : check protection of simplex {p,f}
          // the new simplex is guaranteed to be finite
          vertices.clear(); vertices.push_back(p);
          for (int i = 0; i < D+1; ++i)
            if (i != index)
              vertices.push_back(parent_cell->vertex(i)->point());
          Delaunay_vertex vertex_to_check;
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, false, vertex_to_check->point(), delta)) 
            return true;            
        }
    } 
  else
    {
      // Inside of the convex hull is + side. Outside is - side.
      for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
        if (!t.is_infinite(*vh_it))
          vertices.push_back((*vh_it)->point());
      Delaunay_triangulation::Vertex_iterator v_it = t.vertices_begin();
      while (t.is_infinite(v_it) || vertex_is_in_full_cell(v_it, c))
        v_it++;
      Hyperplane_d facet_plane(vertices.begin(), vertices.end(), v_it->point(), CGAL::ON_POSITIVE_SIDE);
      //CGAL::Oriented_side outside = Oriented_side_d()(facet_plane, v_it->point());
      Vector_d orth_v = facet_plane.orthogonal_vector();
      std::vector<FT> coords;
      auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
      for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
        coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
      Point_d p_delta = Point_d(coords);
      bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
      bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);

      if (!p_is_inside && p_delta_is_inside)
        return true;
      //if the cell is infinite we look at the neighbours regardless
      c->tds_data().mark_visited();
      marked_cells.push_back(c);
      if (p_is_inside)
        for (int i = 0; i < D+1; ++i)
          {
            Full_cell_handle next_c = c->neighbor(i);
            if (next_c->tds_data().is_clear() &&
                is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells))
              return true;
          }
      else
        {
          // facet f is on the border of the conflict zone : check protection of simplex {p,f}
          // the new simplex is finite if the parent cell is finite
          vertices.clear(); vertices.push_back(p);
          bool new_simplex_is_finite = false;
          for (int i = 0; i < D+1; ++i)
            if (i != index)
              {
                if (t.is_infinite(parent_cell->vertex(i)))
                  new_simplex_is_finite = true;
                else
                  vertices.push_back(parent_cell->vertex(i)->point());
              }
          Delaunay_vertex vertex_to_check;
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, new_simplex_is_finite, vertex_to_check->point(), delta)) 
            return true;
        }
    }
  return false;
}

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta)
{
  Euclidean_distance ed;
  Delaunay_triangulation::Vertex_handle v;
  Delaunay_triangulation::Face f(t.current_dimension()); 
  Delaunay_triangulation::Facet ft; 
  Delaunay_triangulation::Full_cell_handle c; 
  Delaunay_triangulation::Locate_type lt;
  std::vector<Full_cell_handle> marked_cells;
  c = t.locate(p, lt, f, ft, v);
  bool violation_existing_cells = is_violating_protection(p, t, c, c, 0, D, delta, marked_cells);
  for (Full_cell_handle fc : marked_cells)
    fc->tds_data().clear();
  return violation_existing_cells;
}

bool old_is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta)
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
        Sphere_d cs( vertices.begin(), vertices.end());
        Point_d center_cs = cs.center();
        FT r = sqrt(ed.transformed_distance(center_cs, fc_it->vertex(1)->point()));
        FT dist2 = ed.transformed_distance(center_cs, p);
        //if the new point is inside the protection ball of a non conflicting simplex
        if (dist2 >= r*r && dist2 <= (r+delta)*(r+delta)) 
          return true;
      }
  t.insert(p, c);
  return false;
}

void write_delaunay_mesh(Delaunay_triangulation& t, const Point_d& p)
{
  std::ofstream ofs ("delaunay.mesh", std::ofstream::out);
  int nbV = t.number_of_vertices()+1;
  ofs << "MeshVersionFormatted 1\nDimension 2\n";
  ofs << "Vertices\n" << nbV << "\n";
  int ind = 1; //index of a vertex
  std::map<Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
  for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
    {
      if (t.is_infinite(v_it))
        continue;
      for (auto coord = v_it->point().cartesian_begin(); coord != v_it->point().cartesian_end(); ++coord)
        ofs << *coord << " ";
      ofs << "508\n";
      index_of_vertex[v_it] = ind++;
    }
  for (auto coord = p.cartesian_begin(); coord != p.cartesian_end(); ++coord)
    ofs << *coord << " ";
  ofs << "208\n";
  /*
  int nbFacets = 0;
  for (auto ft_it = t.finite_facets_begin(); ft_it != t.finite_facets_end(); ++ft_it)
    nbFacets++;
  ofs << "\nEdges\n" << nbFacets << "\n\n";
  for (auto ft_it = t.facets_begin(); ft_it != t.facets_end(); ++ft_it)
    {
      if (t.is_infinite(ft_it))
        continue;
      for (auto vh_it = ft_it->vertices_begin(); vh_it != ft_it->vertices_end(); ++vh_it)
        ofs << index_of_vertex[*vh_it] << " "; 
    }
  */
  ofs << "Triangles " << t.number_of_finite_full_cells()+1 << "\n";
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    {
      if (t.is_infinite(fc_it))
        continue;
      for (auto vh_it = fc_it->vertices_begin(); vh_it != fc_it->vertices_end(); ++vh_it)
        ofs << index_of_vertex[*vh_it] << " ";
      ofs << "508\n";
    }
  ofs << nbV << " " << nbV << " " << nbV << " " << 208 << "\n";
  ofs << "End\n";
  ofs.close();
}

bool triangulation_is_protected(Delaunay_triangulation& t, FT delta)
{
  // Verification part
  Euclidean_distance ed;
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    if (!t.is_infinite(fc_it))
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it))
          //check if vertex belongs to the face
          if (!vertex_is_in_full_cell(v_it, fc_it))
            {
              std::vector<Point_d> vertices;
              for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
                vertices.push_back((*fc_v_it)->point());
              Sphere_d cs( vertices.begin(), vertices.end());
              Point_d center_cs = cs.center();
              FT r = sqrt(ed.transformed_distance(center_cs, fc_it->vertex(0)->point()));
              FT dist2 = ed.transformed_distance(center_cs, v_it->point());
              //if the new point is inside the protection ball of a non conflicting simplex
              //std::cout << "Dist^2 = " << dist2 << " (r+delta)*(r+delta) = " << (r+delta)*(r+delta) << " r^2 = " << r*r <<"\n";
              //if (dist2 <= (r+delta)*(r+delta) && dist2 >= r*r)
              if (dist2 <= r*r+delta*delta && dist2 >= r*r)
                {
                  write_delaunay_mesh(t, v_it->point());
                  std::cout << "Problematic vertex " << *v_it << " ";
                  std::cout << "Problematic cell " << *fc_it << "\n";
                  std::cout << "r^2 = " << r*r << ", d^2 = " << dist2 << ", r^2+delta^2 = " << r*r+delta*delta << "\n";
                  return false;
                }
            }
        
  return true;
}

void fill_landmarks(Point_Vector& W, Point_Vector& landmarks, std::vector<int>& landmarks_ind)
{
  for (unsigned j = 0; j < landmarks_ind.size(); ++j)
    landmarks.push_back(W[landmarks_ind[j]]);
}

void fill_full_cell_vector(Delaunay_triangulation& t, std::vector<std::vector<int>>& full_cells)
{
  // Store vertex indices in a map
  int ind = 0; //index of a vertex
  std::map<Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
  for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
    if (t.is_infinite(v_it))
      continue;
    else
      index_of_vertex[v_it] = ind++;
  // Write full cells as vectors in full_cells
  for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it)
    {
      if (t.is_infinite(fc_it))
        continue;
      std::vector<int> cell;
      for (auto v_it = fc_it->vertices_begin(); v_it != fc_it->vertices_end(); ++v_it)
        cell.push_back(index_of_vertex[*v_it]);
      full_cells.push_back(cell);
    }
}

FT sampling_radius(Delaunay_triangulation& t)
{
  FT epsilon2 = 4.0;
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
      if (epsilon2 > r2)
        epsilon2 = r2;
    }
  return sqrt(epsilon2);
}

FT point_sampling_radius_by_delaunay(Point_Vector& points)
{
  Delaunay_triangulation t(points[0].size());
  t.insert(points.begin(), points.end());
  return sampling_radius(t);
}

void landmark_choice_protected_delaunay(Point_Vector& W, int nbP, Point_Vector& landmarks, std::vector<int>& landmarks_ind, FT delta, std::vector<std::vector<int>>& full_cells)
{
  unsigned D = W[0].size();
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
    //CGAL::spatial_sort(temp_vector.begin(), temp_vector.end());
    for (std::vector<int>::iterator it = temp_vector.begin(); it != temp_vector.end(); ++it)
      index_list.push_front(*it);
  }
  for (unsigned pos1 = 0; pos1 < D+1; ++pos1)
    {
      std::vector<FT> point;
      for (unsigned i = 0; i < pos1; ++i)
        point.push_back(-1);
      if (pos1 != D)
        point.push_back(1);
      for (unsigned i = pos1+1; i < D; ++i)
        point.push_back(0);
      assert(point.size() == D);
      W[index_list.front()] = Point_d(point);
      insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count);
      index_list.pop_front();
    }
  //      add the first D+1 vertices to form one finite cell
  /*
  for (int i = 0; i <= D+1; ++i)
    {
      t.insert(W[index_list.front()]);
      insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count);
      index_list.pop_front();
    }
  */
  /*
  {   
    std::vector<FT> coords;
    for (int i = 0; i < D; ++i)
      coords.push_back(-1);
    W[index_list.front()] = Point_d(coords);
    insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count);
    index_list.pop_front();
    for (int i = 0; i < D; ++i)
      {
        coords.clear();
        for (int j = 0; j < D; ++j)
          if (i == j)
            coords.push_back(1);
          else
            coords.push_back(-1);
        W[index_list.front()] = Point_d(coords);
        insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count);
        index_list.pop_front();
      }
  }
  */
  //std::cout << t;
  //assert(t.number_of_vertices() == D+1);
  //assert(landmarks_ind.size() == D+1);
  //assert(W[landmarks_ind[0]][0] == 0);
  //      add other vertices if they don't violate protection
  std::list<int>::iterator list_it = index_list.begin();
  while (list_it != index_list.end())
    {
      if (!is_violating_protection(W[*list_it], t, D, delta))
        {
          // If no conflicts then insert in every copy of T^3      
          is_violating_protection(W[*list_it], t, D, delta);
          insert_delaunay_landmark_with_copies(W, *list_it, landmarks_ind, t, landmark_count);
          index_list.erase(list_it);
          list_it = index_list.begin();
          //std::cout << "index_list_size() = " << index_list.size() << "\n";
        }
      else
        {
          list_it++;
          //std::cout << "!!!!!WARNING!!!!! A POINT HAS BEEN OMITTED!!!\n";
        }
      //write_delaunay_mesh(t, W[*list_it]);
    }
  fill_landmarks(W, landmarks, landmarks_ind);
  fill_full_cell_vector(t, full_cells);
  if (triangulation_is_protected(t, delta))
    std::cout << "Triangulation is ok\n";
  else
    std::cout << "Triangulation is BAD!! T_T しくしく!\n";
  write_delaunay_mesh(t, Point_d(std::vector<FT>({0,0})));
  //std::cout << t << std::endl;
}

template <typename T>
void print_vector(std::vector<T> v)
{
  std::cout << "[";
  if (!v.empty())
    {
      std::cout << *(v.begin());
      for (auto it = v.begin()+1; it != v.end(); ++it)
        {
          std::cout << ",";
          std::cout << *it;
        }
    }
  std::cout << "]";
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
  write_wl(out_file,WL);
  
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
  /*
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
  */
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
  
  write_edges("landmarks/edges", witnessComplex, landmarks);
  /*
  return count_badlinks;
  */
  return 0;
}


int main (int argc, char * const argv[])
{
  if (argc != 4)
    {
      std::cerr << "Usage: " << argv[0]
                << " nbP dim delta\n";
      return 0;
    }
  int nbP       = atoi(argv[1]);
  int dim       = atoi(argv[2]);
  double delta  = atof(argv[3]);
  
  std::cout << "Let the carnage begin!\n";
  Point_Vector point_vector;
  generate_points_random_box(point_vector, nbP, dim);
  FT epsilon = point_sampling_radius_by_delaunay(point_vector);
  std::cout << "Initial epsilon = " << epsilon << std::endl;
  Point_Vector L;
  std::vector<int> chosen_landmarks;
  //write_points("landmarks/initial_pointset",point_vector);
  //write_points("landmarks/initial_landmarks",L);
  CGAL::Timer timer;
  /*
  for (int i = 0; i < 11; i++)
  //for (int i = 0; bl > 0; i++)
    {
      //std::cout << "========== Start iteration " << i << "== curr_min(" << curr_min << ")========\n";
      double delta = pow(10, -(1.0*i)/2);
      std::cout << "delta = " << delta << std::endl;
      L = {}; chosen_landmarks = {};
      std::vector<std::vector<int>> full_cells;
      timer.start();
      landmark_choice_protected_delaunay(point_vector, nbP, L, chosen_landmarks, delta, full_cells);
      timer.stop();
      FT epsilon2 = point_sampling_radius_by_delaunay(L);
      std::cout << "Final epsilon = " << epsilon2  << ". Ratio = " << epsilon/epsilon2 << std::endl;
      write_points("landmarks/initial_landmarks",L);
      int nbL = chosen_landmarks.size();
      std::cout << "Number of landmarks = " << nbL << ", time= " << timer.time() << "s"<< std::endl;
      landmark_perturbation(point_vector, nbL, L, chosen_landmarks, full_cells);
      timer.reset();
      //write_points("landmarks/landmarks0",L);
    }
  */
  std::vector<std::vector<int>> full_cells;
  timer.start();
  landmark_choice_protected_delaunay(point_vector, nbP, L, chosen_landmarks, delta, full_cells);
  timer.stop();
  FT epsilon2 = point_sampling_radius_by_delaunay(L);
  std::cout << "Final epsilon = " << epsilon2  << ". Ratio = " << epsilon/epsilon2 << std::endl;
  write_points("landmarks/initial_landmarks",L);
  int nbL = chosen_landmarks.size();
  std::cout << "Number of landmarks = " << nbL << ", time= " << timer.time() << "s"<< std::endl;
  landmark_perturbation(point_vector, nbL, L, chosen_landmarks, full_cells);
  timer.reset();
}
