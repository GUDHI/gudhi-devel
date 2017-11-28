#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#include <gudhi/Points_off_io.h>
#include <gudhi/Ad_simplex.h>

#include <CGAL/Epick_d.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cxx-prettyprint/prettyprint.hpp"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;

using Matrix = Eigen::SparseMatrix<FT>;
using Triplet = Eigen::Triplet<FT>;
using Simplex_id = std::vector<int>;
using Vertex_id = Simplex_id;
using Pointer_range = typename std::vector<Point_vector::iterator>;


struct Lexicographic {
  bool operator() (const Simplex_id &lhs, const Simplex_id &rhs) {
    assert (lhs.size() == rhs.size());
    auto l_it = lhs.begin();
    auto r_it = rhs.begin();
    for (; l_it != lhs.end(); ++l_it, ++r_it)
      if (*l_it < *r_it)
        return true;
      else if (*l_it > *r_it)
        return false;
    return false;
  }
};

using SPMap = std::map<Simplex_id, Pointer_range, Lexicographic>;
using SPointer_range = typename std::vector<SPMap::iterator>;
using VSMap = std::map<Vertex_id, std::vector<int>, Lexicographic>;

struct Lexicographic_ptr {
  bool operator() (const SPMap::iterator &lhs, const SPMap::iterator &rhs) {
    Lexicographic lx;
    return lx(lhs->first,rhs->first);
  }
};

using SiMap = std::map<SPMap::iterator, int, Lexicographic_ptr>;


int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

/** Common gcd simplification */
template <class Id>
Id reduced_id(Id& id) {
  int common_gcd = 0;
  for (auto i: id) {
    common_gcd = gcd(i, common_gcd);
    if (common_gcd == 1)
      return id;
  }
  Id id_red(id);
  for (auto i_it = id_red.begin(); i_it != id_red.end(); ++i_it) {
    *i_it = *i_it / common_gcd;
  }
  return id_red;
}

/** A conversion from Cartesian coordinates to root coordinates.
 *  The matrix' rows are root vectors (or normal vectors of a simplex in general).
 */
template <class Point,
          class Matrix>
Point root_coordinates(Point p, Matrix& root_t, short d)
{
  // short d = p.size();
  std::vector<double> p_r;
  for (int i = 0; i < d+1; i++) {
    FT sc_prod = 0;
    /* for now no root normalization takes place */
    // FT root_norm_sq = 0;
    // for (int j = 0; j < d; j++)
    //   root_norm_sq += root_t.coeff(i,j)*root_t.coeff(i,j);
    // FT root_norm = sqrt()
    for (int j = 0; j < d; j++) {
      sc_prod += root_t.coeff(i,j) * p[j];
    }
    p_r.push_back(sc_prod);
  }
  return Point(p_r);
}

/** A conversion from Cartesian coordinates to root coordinates in a point range.
 *  The matrix' rows are root vectors (or normal vectors of a simplex in general).
 *  The input point range is rewritten.
 */
template <class Point_list,
          class Matrix>
Point_list root_coordinates_range(Point_list& points, Matrix& root_t)
{
  short d = points[0].size();
  Point_list points_r;
  for (auto p: points) {
    points_r.push_back(root_coordinates(p,root_t,d));
  }
  return points_r;
}


/** Current state of the algorithm.
 *  Input: a point cloud 'point_vector'
 *  Output: a reconstruction (a simplicial complex?, a Czech-like complex?)
 */

int main(int argc, char * const argv[]) {
  std::cout << "Marching cube adaptation for Coxeter triangulations\n";
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0]
        << " path_to_off_point_file\n";
    return 0;
  }
  Point_vector point_vector;
  Gudhi::Points_off_reader<Point_d> off_reader(argv[1]);
  if (!off_reader.is_valid()) {
      std::cerr << "Coxeter triangulations - Unable to read file " << argv[1] << "\n";
      exit(-1);  // ----- >>
    }
  point_vector = Point_vector(off_reader.get_point_cloud());
  int N = point_vector.size();
  short d = point_vector[0].size();
  // short d = 2;
  std::cout << "Successfully read " << N << " points in dimension " << d << std::endl;

  // The A root vectors, computed as a matrix

  std::vector<Triplet> cartan_triplets;
  cartan_triplets.reserve(3*d-2);
  for (int i = 0; i < d; i++) {
    cartan_triplets.push_back(Triplet(i,i,2.0));
  }
  for (int i = 1; i < d; i++) {
    cartan_triplets.push_back(Triplet(i-1,i,-1.0));
    cartan_triplets.push_back(Triplet(i,i-1,-1.0));
  }
  Matrix cartan(d,d);
  cartan.setFromTriplets(cartan_triplets.begin(), cartan_triplets.end());
  std::cout << "cartan =" << std::endl << cartan << std::endl;

  Eigen::SimplicialLLT<Matrix, Eigen::Lower> chol(cartan);
  Matrix do_vect = chol.matrixL();
  std::cout << "do_vect^t =" << std::endl << do_vect << std::endl;
  
  std::vector<Triplet> r_t_triplets;
  r_t_triplets.reserve(2*d);
  for (int i = 0; i < d; i++) {
    r_t_triplets.push_back(Triplet(i,i,1.0));
  }
  for (int i = 0; i < d; i++) {
    r_t_triplets.push_back(Triplet(d,i,-1.0));
  }
  Matrix r_t(d+1,d);
  r_t.setFromTriplets(r_t_triplets.begin(), r_t_triplets.end()); 
  std::cout << "r_t =" << std::endl << r_t << std::endl;

  Matrix root_t = r_t * do_vect;
  std::cout << "norm_t =" << std::endl << root_t << std::endl;
  
  // Compute the root coordinates the root matrix
  // std::cout << "First point is:";
  // for (auto x: point_vector[0])
  //   std::cout << " " << x;
  point_vector = root_coordinates_range(point_vector, root_t);
  // std::cout << ", the root coordinates are";
  // Point_d p_r = root_coordinates(point_vector[0], root_t, d);
  // for (auto x: p_r)
  //   std::cout << " " << x;
  // std::cout << ".\n";

  // The first fill of a map: simplex coordinates -> points
  std::cout << point_vector[0] << std::endl;
  for (auto x: point_vector[0])
    std::cout << std::floor(x) << " ";
  std::cout << std::endl;
  SPMap sp_map;
  for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
    Simplex_id s_id(1,1);
    for (auto x: *p_it)
      s_id.push_back(std::floor(x));
    auto find_it = sp_map.find(s_id);
    if (find_it == sp_map.end())
      sp_map.emplace(s_id, Pointer_range(1, p_it));
    else
      find_it->second.push_back(p_it);
  }
  // std::cout << "SPMap composition:\n";
  // for (auto m: sp_map) {
  //   std::cout << m.first << ": " << m.second.size() << " elements.\n";
  // }

  // small test
  Simplex_id p1 = {4,12,10,-8};
  std::cout << "Non-reduced: " << p1 << ", reduced: " << reduced_id(p1) << ".\n";

  SiMap si_map;
  int si_index = 0;
  for (auto m_it = sp_map.begin(); m_it != sp_map.end(); ++m_it, si_index++)
    si_map.emplace(m_it, si_index);

  std::cout << "SIMap composition:\n";
  for (auto m: si_map) {
    std::cout << m.first->first << ": index " << m.second << ".\n";
  }
  
  // map : vertex coordinates -> simplex coordinates
  VSMap vs_map;
  for (auto m_it = sp_map.begin(); m_it != sp_map.end(); ++m_it) {
    for (int i = 1; i <= d; i++) {
      Vertex_id v_id(m_it->first);
      v_id[i] += 1;
      v_id = reduced_id(v_id);
      auto find_it = vs_map.find(v_id);
      if (find_it == vs_map.end())
        vs_map.emplace(v_id, std::vector<int>(1, si_map[m_it]));
      else
        find_it->second.push_back(si_map[m_it]);
    }
  }
  std::cout << "VSMap composition:\n";
  for (auto m: vs_map) {
    std::cout << m.first << ": " << m.second << ".\n";
  }
  
}
