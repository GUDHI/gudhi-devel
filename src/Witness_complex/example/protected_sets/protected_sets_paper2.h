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
#include "../output.h"
#include "../generators.h"

#include <CGAL/point_generators_d.h>


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

typedef CGAL::Random_points_in_ball_d<Point_d> Random_point_iterator;


FT  _sfty = pow(10,-14);

bool experiment1, experiment2, experiment3, experiment5 = false;

/* Experiment 1: epsilon as function on time **********************/
std::vector<FT> eps_vector;

/* Experiment 2: R/epsilon on alpha *******************************/
std::vector<FT> epsratio_vector;
std::vector<FT> epsslope_vector;

/* Experiment 3: theta on delta ***********************************/
std::vector<FT> thetamin_vector; FT curr_theta;
std::vector<FT> gammamin_vector;

/* Statistical data ***********************************************/
int refused_case1, refused_case2, refused_bad, refused_centers1, refused_centers2;

void initialize_statistics()
{
  refused_case1 = 0;
  refused_case2 = 0;
  refused_bad   = 0;
  refused_centers1 = 0;
  refused_centers2 = 0;
}

void print_statistics()
{
  std::cout << " * Old simplex not protected: " << refused_case1 << "\n";
  std::cout << " * New simplex not protected: " << refused_case2 << "\n";
  std::cout << " * New simplex not good:      " << refused_bad << "\n";
  std::cout << " * New-old centers too close:      " << refused_centers1 << "\n";
  std::cout << " * New-new centers too close:      " << refused_centers2 << "\n";
}

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
          if (experiment3 && thetamin_vector[thetamin_vector.size()-1] > pow(h/(2*r), 1.0/(D-1)))
            {
              thetamin_vector[thetamin_vector.size()-1] = pow(h/(2*r), 1.0/(D-1));
              //std::cout << "theta=" << h/(2*r) << ", ";
            }
          if (h/(2*r) < pow(theta0, D-1))
              return false;
          if (!is_theta0_good(facet, theta0))
            return false;
        }
    }
  return true;
}

/** Recursive function for checking the goodness of a simplex,
 *  meaning it does not contain a k-face, which is not theta0^(k-1) thick 
 */

FT theta(std::vector<Point_d>& vertices)
{
  FT curr_value = 1.0;
  int D = vertices.size()-1;
  if (D <= 1)
    return 1; // Edges are always good
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
      curr_value = std::min(curr_value, theta(facet)); // Check the corresponding facet
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
          curr_value = std::min(curr_value, std::pow(h/(2*r), 1.0/(D-1)));
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
          curr_value = std::min(curr_value, pow(h/(2*r), 1.0/(D-1)));
        }
    }
  return curr_value;
}

// Doubling in a way 1->2->5->10
void double_round(int& i)
{
  FT order10 = pow(10,std::floor(std::log10(i)));
  int digit = std::floor( i / order10);
  std::cout << digit;
  if (digit == 1)
    i *= 2;
  else if (digit == 2)
    i = 5*i/2;
  else if (digit == 5)
    i *= 2;
  else
    std::cout << "digit not correct. digit = " << digit << std::endl;
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

bool new_cell_is_violated(Delaunay_triangulation& t, std::vector<Point_d>& vertices, const Delaunay_vertex& v1, const Delaunay_vertex v2, FT delta0, bool power_protection, FT theta0, FT gamma0)
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
            if (dist2 >= r*r-_sfty && dist2 <= (r+r*delta0)*(r+r*delta0))
              { refused_case2++; return true;}
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+r*r*delta0*delta0) 
              { refused_case2++; return true;}
          // Check if the centers are not too close
          std::vector<Point_d> sigma(vertices);
          sigma[0] = v1->point();
          Sphere_d cs_sigma(sigma.begin(), sigma.end());
          Point_d csc_sigma = cs_sigma.center();
          FT r_sigma = sqrt(cs_sigma.squared_radius());
          FT dcc = sqrt(Euclidean_distance().transformed_distance(center_cs, csc_sigma));
          if (experiment3 && dcc/r < gammamin_vector[gammamin_vector.size()-1])
            gammamin_vector[gammamin_vector.size()-1] = dcc/r;
          if (experiment3 && dcc/r_sigma < gammamin_vector[gammamin_vector.size()-1])
            gammamin_vector[gammamin_vector.size()-1] = dcc/r_sigma;
          if (dcc < r*gamma0 || dcc < r_sigma*gamma0)
            { refused_centers1++; return true; }
        }
      if (!t.is_infinite(v2))
        {
          FT dist2 = Euclidean_distance().transformed_distance(center_cs, v2->point());
          if (!power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= (r+r*delta0)*(r+r*delta0))
              { refused_case2++; return true;}
          if (power_protection)
            if (dist2 >= r*r-_sfty && dist2 <= r*r+r*r*delta0*delta0) 
              { refused_case2++; return true;}
          // Check if the centers are not too close
          std::vector<Point_d> sigma(vertices);
          sigma[0] = v2->point();
          Sphere_d cs_sigma(sigma.begin(), sigma.end());
          Point_d csc_sigma = cs_sigma.center();
          FT r_sigma = sqrt(cs_sigma.squared_radius());
          FT dcc = sqrt(Euclidean_distance().transformed_distance(center_cs, csc_sigma));
          if (experiment3 && dcc/r < gammamin_vector[gammamin_vector.size()-1])
            gammamin_vector[gammamin_vector.size()-1] = dcc/r;
          if (experiment3 && dcc/r_sigma < gammamin_vector[gammamin_vector.size()-1])
            gammamin_vector[gammamin_vector.size()-1] = dcc/r_sigma;
          if (dcc < r*gamma0 || dcc < r_sigma*gamma0)
            { refused_centers1++; return true; }
        }
      // Check if the simplex is theta0-good
      if (!is_theta0_good(vertices, theta0))
        { refused_bad++; return true;}
      
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
            coords.push_back((*p_i) - (*orth_i) * delta0 / sqrt(orth_v.squared_length()));
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
            coords.push_back((*p_i) - (*orth_i) * delta0 / sqrt(orth_v.squared_length()));
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

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, Full_cell_handle c, Full_cell_handle parent_cell, int index, int D, FT delta0, std::vector<Full_cell_handle>& marked_cells, bool power_protection, FT theta0, FT gamma0)
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
        if (dist2 >= r*r-_sfty && dist2 <= (r+r*delta0)*(r+r*delta0))
          { refused_case1++; return true;}
      if (power_protection)
        if (dist2 >= r*r-_sfty && dist2 <= r*r+delta0*delta0*r*r) 
          { refused_case1++; return true;}    
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
                  is_violating_protection(p, t, next_c, c, i, D, delta0, marked_cells, power_protection, theta0, gamma0))
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
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta0, power_protection, theta0, gamma0)) 
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
        coords.push_back((*p_i) - (*orth_i) * delta0 / sqrt(orth_v.squared_length()));
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
                  is_violating_protection(p, t, next_c, c, i, D, delta0, marked_cells, power_protection, theta0, gamma0))
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
          if (new_cell_is_violated(t, vertices, vertex_to_check, parent_cell->vertex(index), delta0, power_protection, theta0, gamma0)) 
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

bool is_violating_protection(Point_d& p, Delaunay_triangulation& t, int D, FT delta0, bool power_protection, FT theta0, FT gamma0)
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
  bool violation_existing_cells = is_violating_protection(p, t, c, c, 0, D, delta0, marked_cells, power_protection, theta0, gamma0);
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

  
void initialize_torus(Point_Vector& W, Tree& treeW, Delaunay_triangulation& t, FT epsilon, std::vector<int>& landmarks_ind, int& landmark_count, std::vector<bool>& point_taken)
{
  initialize_statistics();
  int D = W[0].size();
  if (D == 2)
    {
      int xw = 6, yw = 4;
      // Triangular lattice close to regular triangles h=0.866a ~ 0.875a   : 48p
      for (int i = 0; i < xw; ++i)
        for (int j = 0; j < yw; ++j)
          {
            Point_d cite1(std::vector<FT>{2.0/xw*i, 2.0/yw*j});
            std::vector<int> points_inside;
            torus_search(treeW, D, cite1, epsilon, points_inside);
            //std::cout << "i=" << i << ", j=" << j << " "; print_vector(points_inside); std::cout << "\n";
            std::vector<int>::iterator p_it = points_inside.begin();
            while (p_it != points_inside.end() && point_taken[*p_it])
              ++p_it;
            assert(p_it != points_inside.end());
            //W[*p_it] = cite1;  // debug purpose
            insert_delaunay_landmark_with_copies(W, *p_it,
                                                 landmarks_ind, t, landmark_count, true);
            point_taken[*p_it] = true;

            Point_d cite2(std::vector<FT>{2.0/xw*(i+0.5), 2.0/yw*(j+0.5)});
            points_inside.clear();
            torus_search(treeW, D, cite2, epsilon, points_inside);
            //std::cout << "i=" << i << ", j=" << j << " "; print_vector(points_inside); std::cout << "\n";
            p_it = points_inside.begin();
            while (p_it != points_inside.end() && point_taken[*p_it])
              ++p_it;
            assert(p_it != points_inside.end());
            //W[*p_it] = cite2;  // debug purpose
            insert_delaunay_landmark_with_copies(W, *p_it,
                                                 landmarks_ind, t, landmark_count, true);
            point_taken[*p_it] = true;            
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
              std::vector<int>::iterator p_it = points_inside.begin();
              while (p_it != points_inside.end() && point_taken[*p_it])
                ++p_it;
              assert(p_it != points_inside.end());
              insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                   landmarks_ind, t, landmark_count, true);
              point_taken[*p_it] = true;
              
              Point_d cite2(std::vector<FT>{2.0/wd*(i+0.5), 2.0/wd*(j+0.5), 2.0/wd*(k+0.5)});
              points_inside.clear();
              torus_search(treeW, D, cite2, epsilon, points_inside);
              p_it = points_inside.begin();
              while (p_it != points_inside.end() && point_taken[*p_it])
                ++p_it;
              assert(p_it != points_inside.end());
              insert_delaunay_landmark_with_copies(W, *(points_inside.begin()),
                                                   landmarks_ind, t, landmark_count, true);
              point_taken[*p_it] = true;
            }
    }
  //write_mesh
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

FT sampling_fatness(Delaunay_triangulation& t)
{
  FT curr_theta = 1.0;
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
      FT theta_f = theta(vertices);
      curr_theta = std::min(curr_theta, theta_f);
      //std::cout << "theta(sigma) = " << theta_f << "\n";
    }
  return curr_theta;
}

// Generate an epsilon sample for a given epsilon
void generate_epsilon_sample_torus(Point_Vector& W, FT epsilon, int dim, Delaunay_triangulation& t)
{
  W.clear();
  t.clear();
  int point_count = 0;
  std::vector<int> point_ind;
  // std::vector<FT> coords;
  FT curr_eps = 2*dim;
  // Initialize
  // for (int i = 0; i < dim; ++i)
  //   coords.push_back(-1);
  // R_max_handle rmh(2*sqrt(dim), Point_d(coords));
  // int N = dim; std::floor(std::pow(1/epsilon,dim));
  // std::cout << N << "\n";
  typedef CGAL::Random_points_in_cube_d<Point_d> Random_cube_iterator;
  Random_cube_iterator rp(dim, 1.0);
  W.push_back(*rp++);
  insert_delaunay_landmark_with_copies(W, W.size()-1, point_ind, t, point_count, true);
  curr_eps = sampling_radius(t).value;      
  while (curr_eps > epsilon)
    {
      
      W.push_back(*rp++);
      insert_delaunay_landmark_with_copies(W, W.size()-1, point_ind, t, point_count, true);
      
      Point_d c = sampling_radius(t).center;
      W.push_back(c);
      insert_delaunay_landmark_with_copies(W, W.size()-1, point_ind, t, point_count, true);
      curr_eps = sampling_radius(t).value;
      
      std::cout << "curr_eps = " << curr_eps << "\n";
    }
  // Iterate and insert in a torus
  // while (rmh.value > epsilon)
  //   {
  //     W.push_back(rmh.center);
  //     insert_delaunay_landmark_with_copies(W, W.size()-1, point_ind, t, point_count, true);
  //     rmh = sampling_radius(t);
  //     //std::cout << rmh.value;
  //   }
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
                        FT alpha,
                        FT epsilon,
                        FT delta0,
                        FT theta0,
                        FT gamma0,
                        //std::vector<std::vector<int>>& full_cells,
                        bool torus,
                        bool power_protection
                        )
{
  //bool return_ = true;
  unsigned D = W[0].size();
  int nbP = W.size();
  //FT beta = 1/(1-alpha);
  //FT Ad = pow((4*alpha + 8*beta)/alpha, D);
  //FT theta0 = 1/Ad;
  //FT delta0 = pow(1/Ad,D);
  Torus_distance td;
  Euclidean_distance ed;
  Delaunay_triangulation t(D);
  std::vector<bool> point_taken(nbP,false);
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
    initialize_torus(W, treeW, t, epsilon, landmarks_ind, landmark_count, point_taken);
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
      if (!point_taken[*cp_it] && !is_violating_protection(W[*cp_it], t, D, delta0, power_protection, theta0, gamma0))
        {
          Delaunay_vertex v = insert_delaunay_landmark_with_copies(W, *cp_it, landmarks_ind, t, landmark_count, torus);
          {
            // Simple check if the new cells don't have centers too close one to another
            std::vector<Full_cell_handle> inc_cells;
            std::back_insert_iterator<std::vector<Full_cell_handle>> out(inc_cells);
            t.tds().incident_full_cells(v, out);

            std::vector<Sphere_d> spheres;
            for (auto i_it = inc_cells.begin(); i_it != inc_cells.end(); ++i_it)
              {
                std::vector<Point_d> vertices;
                for (auto v_it = (*i_it)->vertices_begin(); v_it != (*i_it)->vertices_end(); ++v_it)
                  vertices.push_back((*v_it)->point());
                spheres.push_back(Sphere_d(vertices.begin(), vertices.end()));
              }
            for (auto s_it = spheres.begin(); s_it != spheres.end(); ++s_it)
              for (auto t_it = s_it+1; t_it != spheres.end(); ++t_it)
                {
                  FT ddc2 = ed.transformed_distance(s_it->center(),t_it->center());
                  if (ddc2 < gamma0*gamma0*s_it->squared_radius() ||
                      ddc2 < gamma0*gamma0*t_it->squared_radius())
                    { refused_centers2++; }
                }
          }
            
          //std::cout << *cp_it << ",\n";
          //make_heap(t, R_max_heap);
          point_taken[*cp_it] = true;
          rh = sampling_radius(t);
          if (experiment1) eps_vector.push_back(pow(1/rh.value,D));
          //std::cout << "rhvalue = " << rh.value << "\n";
          //std::cout << "D = " << 
          candidate_points.clear();
          torus_search(treeW, D,
                       rh.center,
                       alpha*rh.value,
                       candidate_points);
          cp_it = candidate_points.begin();
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
  if (experiment2) epsslope_vector.push_back( (pow(1/rh.value,D)-pow(1/epsilon0,D))/(landmarks_ind.size() - 48) );
  std::cout << "The iteration ended when cp_count = " << candidate_points.size() << "\n";
  std::cout << "alphaRmax = " << alpha*rh.value << "\n";
  std::cout << "epsilon' = " << rh.value << "\n";
  std::cout << "nbL = " << landmarks_ind.size() << "\n";
  print_statistics();
  //print_vector(landmarks_ind); std::cout << std::endl;
  //std::sort(landmarks_ind.begin(), landmarks_ind.end());
  print_vector(landmarks_ind); std::cout << std::endl;
  if (experiment3) thetamin_vector[thetamin_vector.size()-1] = sampling_fatness(t);
  std::cout << "theta = " << sampling_fatness(t) << "\n";
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
  write_delaunay_mesh(t, W[0], true);
  //std::cout << t << std::endl;
}

void run_experiment5(Point_Vector& W,
                     int D,
                     FT alpha,
                     FT epsilon,
                     FT delta0,
                     FT theta0,
                     FT gamma0,
                     //std::vector<std::vector<int>>& full_cells,
                     bool torus,
                     bool power_protection
                     )
{
  // INITIALIZATION
  Delaunay_triangulation t(D);
  std::vector<int> landmarks_ind;
  int landmark_count = 0;
  initialize_statistics();
  if (D == 2)
    {
      int xw = 6, yw = 4;
      // Triangular lattice close to regular triangles h=0.866a ~ 0.875a   : 48p
      for (int i = 0; i < xw; ++i)
        for (int j = 0; j < yw; ++j)
          {
            Point_d cite1(std::vector<FT>{2.0/xw*i, 2.0/yw*j});
            W.push_back(cite1);  // debug purpose
            insert_delaunay_landmark_with_copies(W, W.size()-1,
                                                 landmarks_ind, t, landmark_count, true);

            Point_d cite2(std::vector<FT>{2.0/xw*(i+0.5), 2.0/yw*(j+0.5)});
            W.push_back(cite2);  // debug purpose
            insert_delaunay_landmark_with_copies(W, W.size()-1,
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
              W.push_back(cite1);  // debug purpose
              insert_delaunay_landmark_with_copies(W, W.size()-1,
                                                   landmarks_ind, t, landmark_count, true);

              Point_d cite2(std::vector<FT>{2.0/wd*(i+0.5), 2.0/wd*(j+0.5), 2.0/wd*(k+0.5)});
              W.push_back(cite2);  // debug purpose
              insert_delaunay_landmark_with_copies(W, W.size()-1,
                                                   landmarks_ind, t, landmark_count, true);
            }
    }

  // ITERATIONS
  R_max_handle rh = sampling_radius(t);
  Point_d rp = *(Random_point_iterator(D, alpha*rh.value));
  int death_count = 0;
  std::cout << "death count " << death_count << " rp = " << rp << "\n";
  while (death_count < 100)
    {
      std::vector<FT> coords;
      for (auto c_it = rh.center.cartesian_begin(),
           r_it = rp.cartesian_begin();
           c_it != rh.center.cartesian_end();
           ++c_it, ++r_it)
        coords.push_back(*c_it + *r_it);
      Point_d new_p(coords);
      if (!is_violating_protection(new_p, t, D, delta0, power_protection, theta0, gamma0))
        {
          W.push_back(new_p);
          insert_delaunay_landmark_with_copies(W, W.size()-1, landmarks_ind, t, landmark_count, torus);
          rh = sampling_radius(t);
          rp = *(Random_point_iterator(D, alpha*rh.value));
          death_count = 0;
          std::cout << "death count " << death_count << " rp = " << rp << "\n";
        }
      else
        {
          rp = *(Random_point_iterator(D, alpha*rh.value));
          death_count++;
          std::cout << "death count " << death_count << " rp = " << rp << "\n";
        }
      //Point_d new_p = (*rp++) + Vector_d;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Series of experiments
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void start_experiments(Point_Vector& W, FT alpha, std::vector<int>& landmarks_ind, FT epsilon)
{
  int experiment_no = 1;
  FT delta0 = 0.1;
  FT theta0 = 0.1;
  FT gamma0 = 0.01;
  std::string suffix;
  //std::cout << "ようこそジプシー我が神秘の部屋へ:\n";
  while (experiment_no != 0)
    {
      std::cout << "Enter experiment no (0 to exit): ";
      std::cin >> experiment_no;
      switch (experiment_no)
        {
        case 1:
          // Experiment 1
          experiment1 = true;         
          eps_vector = {};
          std::cout << "Enter delta0: "; std::cin >> delta0;
          std::cout << "Enter theta0: "; std::cin >> theta0;
          std::cout << "Enter gamma0: "; std::cin >> gamma0;
          protected_delaunay(W, landmarks_ind, alpha, epsilon, delta0, theta0, gamma0, true, true);
          write_tikz_plot(eps_vector,"epstime.tikz");
          experiment1 = false;
          break;
          
        case 2:
          // Experiment 2
          suffix = "";
          experiment2 = true;
          epsratio_vector = {0};
          epsslope_vector = {0};
          std::cout << "File name suffix: ";
          std::cin >> suffix;
          for (FT alpha = 0.01; alpha < 0.999; alpha += 0.01)
            {
              landmarks_ind.clear();
              std::cout << "Test for alpha = " << alpha << "\n";
              protected_delaunay(W, landmarks_ind, alpha, epsilon, delta0, theta0, gamma0, true, true);
            }
          write_tikz_plot(epsratio_vector,"epsratio_alpha." + suffix + ".tex");
          write_tikz_plot(epsslope_vector,"epsslope_alpha." + suffix + ".tex");
          experiment2 = false;
          break;

        case 3:
          // Experiment 3
          experiment3 = true;
          thetamin_vector = {};
          gammamin_vector = {};
          theta0 = 0;
          gamma0 = 0;
          for (FT delta0 = 0; delta0 < 0.999; delta0 += 0.05)
            {
              landmarks_ind.clear();
              thetamin_vector.push_back(1.0); //0.7489 fatness of the initialization
              gammamin_vector.push_back(10);
              std::cout << "Test for delta0 = " << delta0 << "\n";
              protected_delaunay(W, landmarks_ind, alpha, epsilon, delta0, theta0, gamma0, true, true);
            }
          write_tikz_plot(thetamin_vector,"thetamin_delta.tex");
          write_tikz_plot(gammamin_vector,"gammamin_delta.tex");
          experiment3 = false;
          break;

        // case 4:
        //   // Experiment 4
        //   {
        //     int dim;
        //     std::cout << "Enter dimension: ";
        //     std::cin >> dim;
        //     Delaunay_triangulation t(dim);
        //     // for (FT eps = 0.7; eps < 1.1; eps += 0.1)
        //   //     {
        //   //       generate_epsilon_sample_torus(W, eps, dim, t);
        //   //       for (auto v_it = t.vertices_begin(); v_it != t.vertices_end(); ++v_it)
        //   //         {
        //   //           if (t.is_infinite(v_it))
        //   //             continue;
        //   //           bool in_cube = true; 
        //   //           for (auto xi = v_it->cartesian_begin(); xi != v_it->cartesian_end(); ++xi)
        //   //             if (*xi > 1.0 || *xi < -1.0)
        //   //               {
        //   //                 in_cube = false; break;
        //   //               }
        //   //           if (!in_cube)
        //   //             continue;
        //   //           for (auto t.tds().incident_full_cells())
        //   //         }
        //   //       std::cout << "eps = " << eps << ", real epsilon = " << sampling_radius(t).value << "\n";
        //   //     }
        //   // }
        //   break;


        case 5:
          // Experiment 5
          experiment5 = true;
          //     std::cout << "Enter dimension: ";
          //     std::cin >> dim;

          landmarks_ind.clear();
          W.clear();
          run_experiment5(W,  alpha, epsilon, delta0, theta0, gamma0, true, true);
          experiment5 = false;
          break;
        }

    }

}
  
#endif
