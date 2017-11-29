#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#include <gudhi/Points_off_io.h>
#include <gudhi/Ad_simplex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/output_tikz.h>

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


struct Simplex_tree_options_no_persistence {
  typedef Gudhi::linear_indexing_tag Indexing_tag;
  typedef int Vertex_handle;
  typedef FT Filtration_value;
  typedef std::uint32_t Simplex_key;
  static const bool store_key = true;
  static const bool store_filtration = false;
  static constexpr bool contiguous_vertices = true;
};
using Simplex_tree = Gudhi::Simplex_tree<Simplex_tree_options_no_persistence>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

std::vector<FT> bounding_box_dimensions(Point_vector& points) {
  std::vector<FT> lower, upper, difference;
  for (auto x: points[0]) {
    lower.push_back(x);
    upper.push_back(x);
  }
  for (auto p: points)
    for (int i = 0; i < p.size(); i++) {
      if (p[i] < lower[i])
        lower[i] = p[i];
      if (p[i] > upper[i])
        upper[i] = p[i];
    }
  for (int i = 0; i < lower.size(); i++)
    difference.push_back(upper[i]-lower[i]);
  return difference;
}

/** \brief Write triangles (tetrahedra in 3d) of a simplicial complex in a file, compatible with medit.
 *  `landmarks_ind` represents the set of landmark indices in W  
 *  `st` is the Simplex_tree to be visualized,
 *  `shr` is the Simplex_handle_range of simplices in `st` to be visualized
 *  `is2d` should be true if the simplicial complex is 2d, false if 3d
 *  `l_is_v` = landmark is vertex
 */
void write_coxeter_mesh(Point_vector& W, SiMap& map, Matrix& root_t, std::string file_name = "coxeter.mesh")
{
  short d = W[0].size();
  assert(d < 4);

  std::vector<FT> bbox_dimensions = bounding_box_dimensions(W);
  std::cout << bbox_dimensions << "\n";
  
  std::ofstream ofs (file_name, std::ofstream::out);
  if (d <= 2)
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
  else
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
  int num_vertices = W.size(), num_edges = 0, num_triangles = 0, num_tetrahedra = 0;
  // if (d <= 2) {
  //   num_triangles = W.size();
  //   num_edges = W.size()*3;
  //   num_vertices += W.size()*3;
  // }
  // else {
  //   num_tetrahedra = W.size();
  //   num_triangles = W.size()*4;
  //   num_edges = W.size()*6;
  //   num_vertices += W.size()*4;
  // }

  ofs << "Vertices\n" << num_vertices << "\n";
  
  for (auto p: W) {
    for (auto coord = p.cartesian_begin(); coord != p.cartesian_end() && coord != p.cartesian_begin()+3 ; ++coord)
      ofs << *coord << " ";
    ofs << "508\n";
  }
  // for (auto m: map) {
  //   FT denom = m.first->first[0];
  //   Eigen::VectorXf x(d+1);
  //   for (int i = 0; i < d; i++) {
  //     x(i) = m.first->first[i+1]/denom;
  //   }
    
  // }
  
  ofs << "Edges " << num_edges << "\n";
  for (int i = 1; i <= W.size(); i++) 
    ofs << i << " " << i << " " << "100\n";
  ofs << "Triangles " << num_triangles << "\n";
  for (int i = 1; i <= W.size(); i++) 
    ofs << i << " " << i << " " << i << " " << "100\n";
  if (d == 3) {
    ofs << "Tetrahedra " << num_tetrahedra << "\n";  
    for (int i = 1; i <= W.size(); i++) 
      ofs << i << " " << i << " " << i << " " << i << " " << "100\n";
    ofs << "End\n";
    ofs.close();
  }
}

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
  Point_vector rc_point_vector = root_coordinates_range(point_vector, root_t);
  // std::cout << ", the root coordinates are";
  // Point_d p_r = root_coordinates(point_vector[0], root_t, d);
  // for (auto x: p_r)
  //   std::cout << " " << x;
  // std::cout << ".\n";

  // The first fill of a map: simplex coordinates -> points
  std::cout << rc_point_vector[0] << std::endl;
  for (auto x: rc_point_vector[0])
    std::cout << std::floor(x) << " ";
  std::cout << std::endl;
  SPMap sp_map;
  for (auto p_it = rc_point_vector.begin(); p_it != rc_point_vector.end(); ++p_it) {
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

  // simplex tree construction
  Simplex_tree st;
  for (auto m: vs_map) {
    st.insert_simplex_and_subfaces(m.second);
  }
  std::cout << st;

  Persistent_cohomology pcoh(st);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(11);

  pcoh.compute_persistent_cohomology(-0.1);
  pcoh.output_diagram();

  write_coxeter_mesh(point_vector, si_map, root_t, "coxeter.mesh");
}
