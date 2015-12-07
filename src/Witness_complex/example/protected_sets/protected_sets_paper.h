#ifndef PROTECTED_SETS_H
#define PROTECTED_SETS_H

#include <algorithm>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include <CGAL/Kernel_d/Vector_d.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/policies.hpp>

#include "output_tikz.h"

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point_d;
typedef K::Line_d Line_d;
typedef K::Vector_d Vector_d;
typedef K::Oriented_side_d Oriented_side_d;
typedef K::Has_on_positive_side_d Has_on_positive_side_d;
typedef K::Sphere_d Sphere_d;
typedef K::Hyperplane_d Hyperplane_d;

typedef CGAL::Delaunay_triangulation<K> Delaunay_triangulation;
typedef Delaunay_triangulation::Facet Facet;
typedef Delaunay_triangulation::Vertex_handle Delaunay_vertex;
typedef Delaunay_triangulation::Full_cell_handle Full_cell_handle;

typedef std::vector<Point_d> Point_Vector;
typedef CGAL::Euclidean_distance<Traits_base> Euclidean_distance;

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
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_sphere;


FT  _sfty = pow(10,-14);

bool experiment1, experiment2 = false;

/* Experiment 1: epsilon as function on time **********************/
std::vector<FT> eps_vector;

/* Experiment 2: R/epsilon on delta *******************************/
std::vector<FT> epsratio_vector;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// AUXILLARY FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Insert a point in Delaunay triangulation. If you are working in a flat torus, the procedure adds all the 3^d copies in adjacent cubes as well 
 *  
 *  W is the initial point vector
 *  chosen_landmark is the index of the chosen point in W
 *  landmarks_ind is the vector of indices of already chosen points in W
 *  delaunay is the Delaunay triangulation
 *  landmark_count is the current number of chosen vertices
 *  torus is true iff you are working on a flat torus [-1,1]^d
 *  OUT: Vertex handle to the newly inserted point
 */
Delaunay_vertex insert_delaunay_landmark_with_copies(Point_Vector& W, int chosen_landmark, std::vector<int>& landmarks_ind, Delaunay_triangulation& delaunay, int& landmark_count, bool torus)
{
  if (!torus)
    {
      Delaunay_vertex v =delaunay.insert(W[chosen_landmark]);
      landmarks_ind.push_back(chosen_landmark);
      landmark_count++;
      return v;
    }
  else
    {
      int D = W[0].size();
      int nb_cells = pow(3, D);
      Delaunay_vertex v;
      for (int i = 0; i < nb_cells; ++i)
        {
          std::vector<FT> point;
          int cell_i = i;
          for (int l = 0; l < D; ++l)
            {
              point.push_back(W[chosen_landmark][l] + 2.0*(cell_i%3-1));
              cell_i /= 3;
            }
          if (i == nb_cells/2)
            v = delaunay.insert(point); //v = center point
          else
            delaunay.insert(point);
        }
      landmarks_ind.push_back(chosen_landmark);
      landmark_count++;
      return v;
    }
}

/** Small check if the vertex v is in the full cell fc
 */

bool vertex_is_in_full_cell(Delaunay_triangulation::Vertex_handle v, Full_cell_handle fc)
{
  for (auto v_it = fc->vertices_begin(); v_it != fc->vertices_end(); ++v_it)
    if (*v_it == v)
      return true;
  return false;
}

/** Fill chosen point vector from indices with copies if you are working on a flat torus
 *  
 *  IN:  W is the point vector
 *  OUT: landmarks is the output vector
 *  IN:  landmarks_ind is the vector of indices
 *  IN:  torus is true iff you are working on a flat torus [-1,1]^d
 */

void fill_landmarks(Point_Vector& W, Point_Vector& landmarks, std::vector<int>& landmarks_ind, bool torus)
{
  if (!torus)
    for (unsigned j = 0; j < landmarks_ind.size(); ++j)
      landmarks.push_back(W[landmarks_ind[j]]);
  else
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
}

/** Fill a vector of all simplices in the Delaunay triangulation giving integer indices to vertices
 *
 *  IN: t is the Delaunay triangulation
 *  OUT: full_cells is the output vector
 */

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
      std::vector<int> cell;
      for (auto v_it = fc_it->vertices_begin(); v_it != fc_it->vertices_end(); ++v_it)
        cell.push_back(index_of_vertex[*v_it]);
      full_cells.push_back(cell);
    }
}

bool sphere_intersects_cube(Point_d& c, FT r)
{
  bool in_cube = true;
  //  int i = 0, D = p.size();
  for (auto xi = c.cartesian_begin(); xi != c.cartesian_end(); ++xi)
    // if ((*xi < 1.0 || *xi > -1.0) &&
    //     (*xi-r < 1.0 || *xi-r > -1.0) &&
    //     (*xi+r < 1.0 || *xi+r > -1.0))

    if ((*xi-r < -1.0 && *xi+r < -1.0) ||
        (*xi-r > 1.0  && *xi+r >  1.0 ))
      {
        in_cube = false; break;
      }
  return in_cube;
}

/** Recursive function for checking if the simplex is good,
 *  meaning it does not contain a k-face, which is not theta0^(k-1) thick 
 */

bool is_theta0_good(std::vector<Point_d>& vertices, FT theta0)
{
  if (theta0 > 1)
    {
      std::cout << "Warning! theta0 is set > 1\n";
      return false;
    }
  int D = vertices.size()-1;
  if (D <= 1)
    return true; // Edges are always good
  //******** Circumscribed sphere
  Euclidean_distance ed;
  Sphere_d cs(vertices.begin(), vertices.end());
  FT r = sqrt(cs.squared_radius());
  for (std::vector<Point_d>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it)
    {
      std::vector<Point_d> facet;
      for (std::vector<Point_d>::iterator f_it = vertices.begin(); f_it != vertices.end(); ++f_it)
        if (f_it != v_it)
          facet.push_back(*f_it);
      // Compute the altitude

      if (vertices[0].size() == 3 && D == 2)
        {
          //Vector_d l = facet[0] - facet[1];
          FT orth_length2 = ed.transformed_distance(facet[0],facet[1]);
          K::Cartesian_const_iterator_d l_it, p_it, s_it, c_it;
          FT h = 0;
          // Scalar product = <sp,l>
          FT scalar = 0;
          for (p_it = v_it->cartesian_begin(),
               s_it = facet[0].cartesian_begin(),
               l_it = facet[1].cartesian_begin();
               p_it != v_it->cartesian_end();
               ++l_it, ++p_it, ++s_it)
            scalar += (*l_it - *s_it)*(*p_it - *s_it);
          // Gram-Schmidt for one vector
          for (p_it = v_it->cartesian_begin(),
               s_it = facet[0].cartesian_begin(),
               l_it = facet[1].cartesian_begin();
               p_it != v_it->cartesian_end();
               ++l_it, ++p_it, ++s_it)
            {
              FT hx = (*p_it - *s_it) - scalar*(*l_it - *s_it)/orth_length2;
              h += hx*hx;
            }
          h = sqrt(h);
      
          if (h/(2*r) < pow(theta0, D-1))
            return false;
          if (!is_theta0_good(facet, theta0))
            return false;
        }
      else
        {
      Hyperplane_d tau_h(facet.begin(), facet.end(), *v_it);
      Vector_d orth_tau = tau_h.orthogonal_vector();
      FT orth_length = sqrt(orth_tau.squared_length());
      K::Cartesian_const_iterator_d o_it, p_it, s_it, c_it;
      FT h = 0;
      for (o_it = orth_tau.cartesian_begin(),
           p_it = v_it->cartesian_begin(),
           s_it = (facet.begin())->cartesian_begin();
           o_it != orth_tau.cartesian_end();
           ++o_it, ++p_it, ++s_it)
        h += (*o_it)*(*p_it - *s_it)/orth_length;
      h = fabs(h);
      if (h/(2*r) < pow(theta0, D-1))
        return false;
      if (!is_theta0_good(facet, theta0))
        return false;
        }
    }
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IS VIOLATED TEST
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Check if a newly created cell is protected from old vertices
 *  
 *  t is the Delaunay triangulation
 *  vertices is the vector containing the point to insert and a facet f in t
 *  v1 is the vertex of t, such that f and v1 form a simplex
 *  v2 is the vertex of t, such that f and v2 form another simplex
 *  delta is the protection constant
 *  power_protection is true iff the delta-power protection is used
 */

bool new_cell_is_violated(Delaunay_triangulation& t, std::vector<Point_d>& vertices, const Delaunay_vertex& v1, const Delaunay_vertex v2, FT delta, bool power_protection, FT theta0)
{
  assert(vertices.size() == vertices[0].size() ||
         vertices.size() == vertices[0].size() + 1); //simplex size = d | d+1
  assert(v1 != v2);
  if (vertices.size() == vertices[0].size() + 1)
  // FINITE CASE
    {
      Sphere_d cs(vertices.begin(), vertices.end());
      Point_d center_cs = cs.center();
      FT r = sqrt(Euclidean_distance().transformed_distance(center_cs, vertices[0]));
      /*
      for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        if (!t.is_infinite(v_it)) 
          {
            //CGAL::Oriented_side side = Oriented_side_d()(cs, (v_it)->point()); 
            if (std::find(vertices.begin(), vertices.end(), v_it->point()) == vertices.end())
              {
                FT dist2 = Euclidean_distance().transformed_distance(center_cs, (v_it)->point());
                if (!power_protection)
                  if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
                    return true;
                if (power_protection)
                  if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
                    return true;                 
              }
          }
      */
      // Check if the simplex is theta0-good
      if (!is_theta0_good(vertices, theta0))
        return true;
      // Is the center inside the box? (only Euclidean case)
      // if (!torus)
      //   {
      //     bool inside_the_box = true;
      //     for (c_it = center_cs.cartesian_begin(); c_it != center_cs.cartesian_end(); ++c_it)
      //       if (*c_it > 1.0 || *c_it < -1.0)
      //         {
      //           inside_the_box = false; break;
      //         }
      //     if (inside_the_box && h/r < theta0)
      //       return true;
      //   }
      // Check the two vertices (if not infinite)
      if (!t.is_infinite(v1))
        {
          FT dist2 = Euclidean_distance().transformed_distance(center_cs, v1->point());
          if (!power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
              return true;
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
              return true;
        }
      if (!t.is_infinite(v2))
        {
          FT dist2 = Euclidean_distance().transformed_distance(center_cs, v2->point());
          if (!power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
              return true;
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
              return true;
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
      /*
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
      */
      if (!t.is_infinite(v1))
        {
          std::vector<FT> coords;
          Point_d p = v1->point();
          auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
          for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
            coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
          Point_d p_delta = Point_d(coords);
          bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
          bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
          if (!power_protection && !p_is_inside && p_delta_is_inside)
            return true;
        }
      if (!t.is_infinite(v2))
        {
          std::vector<FT> coords;
          Point_d p = v2->point();
          auto orth_i = orth_v.cartesian_begin(), p_i = p.cartesian_begin();
          for (; orth_i != orth_v.cartesian_end(); ++orth_i, ++p_i)
            coords.push_back((*p_i) - (*orth_i) * delta / sqrt(orth_v.squared_length()));
          Point_d p_delta = Point_d(coords);
          bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p);
          bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);
          if (!power_protection && !p_is_inside && p_delta_is_inside)
            return true;
        }
    }
  return false;
}
  
/** Auxillary recursive function to check if the point p violates the protection of the cell c and
 *  if there is a violation of an eventual new cell
 *
 *  p is the point to insert
 *  t is the current triangulation
 *  c is the current cell (simplex)
 *  parent_cell is the parent cell (simplex)
 *  index is the index of the facet between c and parent_cell from parent_cell's point of view
 *  D is the dimension of the triangulation
 *  delta is the protection constant
 *  marked_cells is the vector of all visited cells containing p in their circumscribed ball
 *  power_protection is true iff you are working with delta-power protection
 *
 *  OUT: true iff inserting p hasn't produced any violation so far
 */

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, Full_cell_handle c, Full_cell_handle parent_cell, int index, int D, FT delta, std::vector<Full_cell_handle>& marked_cells, bool power_protection, FT theta0)
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
      if (!power_protection)
        if (dist2 >= r*r-_sfty && dist2 <= (r+delta)*(r+delta))
          return true;
      if (power_protection)
        if (dist2 >= r*r-_sfty && dist2 <= r*r+delta*delta) 
          return true;    
      // if the new point is inside the circumscribing ball : continue violation searching on neighbours
      //if (dist2 < r*r)
      //if (dist2 < (5*r+delta)*(5*r+delta))
      if (dist2 < r*r)
        {
          c->tds_data().mark_visited();
          marked_cells.push_back(c);
          for (int i = 0; i < D+1; ++i)
            {
              Full_cell_handle next_c = c->neighbor(i);
              if (next_c->tds_data().is_clear() &&
                  is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells, power_protection, theta0))
                return true;
            }
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
          Delaunay_vertex vertex_to_check = t.infinite_vertex();
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta, power_protection, theta0)) 
          //if (new_cell_is_violated(t, vertices, vertex_to_check->point(), delta)) 
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
      bool p_is_inside = !Has_on_positive_side_d()(facet_plane, p) && (Oriented_side_d()(facet_plane, p) != CGAL::ZERO);
      bool p_delta_is_inside = !Has_on_positive_side_d()(facet_plane, p_delta);

      // If we work with power protection, we just ignore any conflicts
      if (!power_protection && !p_is_inside && p_delta_is_inside)
        return true;
      //if the cell is infinite we look at the neighbours regardless
      if (p_is_inside)
        {
          c->tds_data().mark_visited();
          marked_cells.push_back(c);    
          for (int i = 0; i < D+1; ++i)
            {
              Full_cell_handle next_c = c->neighbor(i);
              if (next_c->tds_data().is_clear() &&
                  is_violating_protection(p, t, next_c, c, i, D, delta, marked_cells, power_protection, theta0))
                return true;
            }
        }
      else
        {
          // facet f is on the border of the conflict zone : check protection of simplex {p,f}
          // the new simplex is finite if the parent cell is finite
          vertices.clear(); vertices.push_back(p);
          for (int i = 0; i < D+1; ++i)
            if (i != index)
              if (!t.is_infinite(parent_cell->vertex(i)))
                vertices.push_back(parent_cell->vertex(i)->point());
          Delaunay_vertex vertex_to_check = t.infinite_vertex();
          for (auto vh_it = c->vertices_begin(); vh_it != c->vertices_end(); ++vh_it)
            if (!vertex_is_in_full_cell(*vh_it, parent_cell))
              {
                vertex_to_check = *vh_it; break;
              }
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta, power_protection, theta0)) 
          //if (new_cell_is_violated(t, vertices, vertex_to_check->point(), delta)) 
            return true;
        }
    }
  //c->tds_data().clear_visited();
  //marked_cells.pop_back();
  return false;
}

/** Checks if inserting the point p in t will make conflicts
 *
 *  p is the point to insert
 *  t is the current triangulation
 *  D is the dimension of triangulation
 *  delta is the protection constant
 *  power_protection is true iff you are working with delta-power protection
 *  OUT: true iff inserting p produces a violation of delta-protection.
 */

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta, bool power_protection, FT theta0)
{
  Euclidean_distance ed;
  Delaunay_triangulation::Vertex_handle v;
  Delaunay_triangulation::Face f(t.current_dimension()); 
  Delaunay_triangulation::Facet ft; 
  Delaunay_triangulation::Full_cell_handle c; 
  Delaunay_triangulation::Locate_type lt;
  std::vector<Full_cell_handle> marked_cells;
  //c = t.locate(p, lt, f, ft, v);
  c = t.locate(p);
  bool violation_existing_cells = is_violating_protection(p, t, c, c, 0, D, delta, marked_cells, power_protection, theta0);
  for (Full_cell_handle fc : marked_cells)
    fc->tds_data().clear();
  return violation_existing_cells;
}


////////////////////////////////////////////////////////////////////////
// INITIALIZATION
////////////////////////////////////////////////////////////////////////

// Query for a sphere near a cite in all copies of a torus
// OUT points_inside
void torus_search(Tree& treeW, int D, Point_d cite, FT r, std::vector<int>& points_inside)
{
  int nb_cells = pow(3, D);
  Delaunay_vertex v;
  for (int i = 0; i < nb_cells; ++i)
    {
      std::vector<FT> cite_copy;
      int cell_i = i;
      for (int l = 0; l < D; ++l)
        {
          cite_copy.push_back(cite[l] + 2.0*(cell_i%3-1));
          cell_i /= 3;
        }
      Fuzzy_sphere fs(cite_copy, r, 0, treeW.traits());
      treeW.search(std::insert_iterator<std::vector<int>>(points_inside, points_inside.end()), fs);
    }
}

  
void initialize_torus(Point_Vector& W, Tree& treeW, Delaunay_triangulation& t, FT epsilon, std::vector<int>& landmarks_ind, int& landmark_count)
{
  int D = W[0].size();
  if (D == 2)
    {
      int xw = 6, yw = 4;
      // Triangular lattice close to regular triangles h=0.866a ~ 0.875a   : 48p
      for (int i = 0; i < xw; ++i)
        for (int j = 0; j < yw; ++j)
          {
            Point_d cite1(std::vector<FT>{2.0/xw*i, 1.0/yw*j});
            std::vector<int> points_inside;
            torus_search(treeW, D, cite1, epsilon, points_inside);
            assert(points_inside.size() > 0);
            insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                 landmarks_ind, t, landmark_count, true);
            Point_d cite2(std::vector<FT>{2.0/xw*(i+0.5), 1.0/yw*(j+0.5)});
            points_inside.clear();
            torus_search(treeW, D, cite2, epsilon, points_inside);
            assert(points_inside.size() > 0);
            insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                 landmarks_ind, t, landmark_count, true);            
          }
    }
  else if (D == 3)
    {
      int wd = 3;
      // Body-centered cubic lattice : 54p
      for (int i = 0; i < wd; ++i)
        for (int j = 0; j < wd; ++j)
          for (int k = 0; k < wd; ++k)
            {
              Point_d cite1(std::vector<FT>{2.0/wd*i, 2.0/wd*j, 2.0/wd*k});
              std::vector<int> points_inside;
              torus_search(treeW, D, cite1, epsilon, points_inside);
              assert(points_inside.size() > 0);
              insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                   landmarks_ind, t, landmark_count, true);
              Point_d cite2(std::vector<FT>{2.0/wd*(i+0.5), 2.0/wd*(j+0.5), 2.0/wd*(k+0.5)});
              points_inside.clear();
              torus_search(treeW, D, cite2, epsilon, points_inside);
              assert(points_inside.size() > 0);
              insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                   landmarks_ind, t, landmark_count, true);
            }
    }
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!! THE INTERFACE FOR LANDMARK CHOICE IS BELOW !!!!!!!!!!//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// Struct for R_max_heap elements

struct R_max_handle
{
  FT value;
  Point_d center;
  
  R_max_handle(FT value_, Point_d c): value(value_), center(c)
  {}
};

struct R_max_compare
{
  bool operator()(const R_max_handle& rmh1, const R_max_handle& rmh2) const
  {
    return rmh1.value < rmh2.value;
  }
};

// typedef boost::heap::fibonacci_heap<R_max_handle, boost::heap::compare<R_max_compare>> Heap;

// void make_heap(Delaunay_triangulation& t, Heap& R_max_heap)
// {
//   R_max_heap.clear();
//   for (auto fc_it = t.full_cells_begin(); fc_it != t.full_cells_end(); ++fc_it) 
//     {
//       if (t.is_infinite(fc_it))
//         continue;
//       Point_Vector vertices;
//       for (auto fc_v_it = fc_it->vertices_begin(); fc_v_it != fc_it->vertices_end(); ++fc_v_it)
//         vertices.push_back((*fc_v_it)->point());
//       Sphere_d cs( vertices.begin(), vertices.end());
//       Point_d csc = cs.center();
//       FT r = sqrt(cs.squared_radius());
//       // A ball is in the heap, if it intersects the cube
//       bool accepted = sphere_intersects_cube(csc, sqrt(r)); 
//       if (!accepted)
//         continue;
//       R_max_heap.push(R_max_handle(r, fc_it, csc));
//     }
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// SAMPLING RADIUS 
//////////////////////////////////////////////////////////////////////////////////////////////////////////

R_max_handle sampling_radius(Delaunay_triangulation& t)
{
  FT epsilon2 = 0;
  Point_d final_center;
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
          final_center = csc;         
          control_point = (*vertices.begin());
        }
    }
  return R_max_handle(sqrt(epsilon2), final_center);
}

///////////////////////////////////////////////////////////////////////
// LANDMARK CHOICE PROCEDURE
///////////////////////////////////////////////////////////////////////

/** Procedure to compute a maximal protected subset from a point cloud. All OUTs should be empty at call.
 *
 *  IN:  W is the initial point cloud having type Epick_d<Dynamic_dimension_tag>::Point_d
 *  IN:  nbP is the size of W
 *  OUT: landmarks is the output vector for the points
 *  OUT: landmarks_ind is the output vector for the indices of the selected points in W
 *  IN:  delta is the constant of protection
 *  OUT: full_cells is the output vector of the simplices in the final Delaunay triangulation
 *  IN:  torus is true iff you are working on a flat torus [-1,1]^d
 */ 

void protected_delaunay(Point_Vector& W,
                        //Point_Vector& landmarks,
                        std::vector<int>& landmarks_ind,
                        FT delta,
                        FT epsilon,
                        FT alpha,
                        FT theta0,
                        //std::vector<std::vector<int>>& full_cells,
                        bool torus,
                        bool power_protection
                        )
{
  //bool return_ = true;
  unsigned D = W[0].size();
  int nbP = W.size();
  Torus_distance td;
  Euclidean_distance ed;
  Delaunay_triangulation t(D);
  CGAL::Random rand;
  int landmark_count = 0;
  std::list<int> index_list;
  //****************** Kd Tree W
  STraits traits(&(W[0]));
  Tree treeW(boost::counting_iterator<std::ptrdiff_t>(0),
         boost::counting_iterator<std::ptrdiff_t>(nbP),
         typename Tree::Splitter(),
         traits);
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
  //******************** Initialize point set
  if (!torus)
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
        insert_delaunay_landmark_with_copies(W, index_list.front(), landmarks_ind, t, landmark_count, torus);
        index_list.pop_front();
      }
  else
    initialize_torus(W, treeW, t, epsilon, landmarks_ind, landmark_count);
  //std::cout << "Size of treeW: " << treeW.size() << "\n";
  //std::cout << "Size of t: " << t.number_of_vertices() << "\n";
  //******************* Initialize heap for R_max
  //Heap R_max_heap;
  //make_heap(t, R_max_heap);

  
  R_max_handle rh = sampling_radius(t);
  FT epsilon0 = rh.value; 
  if (experiment1) eps_vector.push_back(pow(1/rh.value,D));
  //******************** Iterative algorithm
  std::vector<int> candidate_points;
  torus_search(treeW, D,
               rh.center,
               alpha*rh.value,
               candidate_points);
  std::list<int>::iterator list_it;
  std::vector<int>::iterator cp_it = candidate_points.begin();
  while (cp_it != candidate_points.end())
    {
      if (!is_violating_protection(W[*cp_it], t, D, delta, power_protection, theta0))
        {
          insert_delaunay_landmark_with_copies(W, *cp_it, landmarks_ind, t, landmark_count, torus);
          //make_heap(t, R_max_heap);
          rh = sampling_radius(t);
          if (experiment1) eps_vector.push_back(pow(1/rh.value,D));
          //std::cout << "rhvalue = " << rh.value << "\n";
          //std::cout << "D = " << 
          candidate_points.clear();
          torus_search(treeW, D,
                       rh.center,
                       alpha*rh.value,
                       candidate_points);
          /*
          // PIECE OF CODE FOR DEBUGGING PURPOSES
          
          Delaunay_vertex inserted_v = insert_delaunay_landmark_with_copies(W, *list_it, landmarks_ind, t, landmark_count);
          if (triangulation_is_protected(t, delta))
            {
              index_list.erase(list_it);
              list_it = index_list.begin();
            }
          else
            { //THAT'S WHERE SOMETHING'S WRONG
              t.remove(inserted_v);
              landmarks_ind.pop_back();
              landmark_count--;
              write_delaunay_mesh(t, W[*list_it], is2d);
              is_violating_protection(W[*list_it], t_old, D, delta); //Called for encore
            }
          */
          //std::cout << "index_list_size() = " << index_list.size() << "\n";
        }
      else
        {
          cp_it++;
          //std::cout << "!!!!!WARNING!!!!! A POINT HAS BEEN OMITTED!!!\n";
        }
      //if (list_it != index_list.end())
      //  write_delaunay_mesh(t, W[*list_it], is2d);
    }
  if (experiment2) epsratio_vector.push_back(rh.value/epsilon0);
  std::cout << "The iteration ended when cp_count = " << candidate_points.size() << "\n";
  std::cout << "alphaRmax = " << alpha*rh.value << "\n";
  std::cout << "epsilon' = " << rh.value << "\n";
  std::cout << "nbL = " << landmarks_ind.size() << "\n";
  //fill_landmarks(W, landmarks, landmarks_ind, torus);
  //fill_full_cell_vector(t, full_cells);
  /*
  if (triangulation_is_protected(t, delta))
    std::cout << "Triangulation is ok\n";
  else
    {
      std::cout << "Triangulation is BAD!! T_T しくしく!\n";
    }
  */
  //write_delaunay_mesh(t, W[0], is2d);
  //std::cout << t << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Series of experiments
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void start_experiments(Point_Vector& W, FT theta0, std::vector<int>& landmarks_ind, FT epsilon)
{
  // Experiment 1
  experiment1 = true;
  protected_delaunay(W, landmarks_ind, 0.1*epsilon, epsilon, 0.5, 0, true, true);
  write_tikz_plot(eps_vector,"epstime.tikz");
  experiment1 = false;

  // Experiment 2
  // experiment2 = true;
  // for (FT delta = 0; delta < epsilon; delta += 0.1*epsilon)
  //   protected_delaunay(W, landmarks_ind, delta, epsilon, 0.5, 0, true, true);
  // write_tikz_plot(epsratio_vector,"epsratio_delta.tikz");
  // experiment2 = false;

}
  
#endif
