#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/output_tikz.h>

#include <CGAL/Epick_d.h>

//#include <Eigen/Dense>
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

using SPMap = std::map<Simplex_id, Pointer_range>;
using SPointer_range = typename std::vector<SPMap::iterator>;
using VSMap = std::map<Vertex_id, std::vector<int>>;

struct Lexicographic_ptr {
  bool operator() (const SPMap::iterator &lhs, const SPMap::iterator &rhs) {
    return std::lexicographical_compare(lhs->first.begin(), lhs->first.end(), rhs->first.begin(), rhs->first.end());
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
using Simplex_tree = Gudhi::Simplex_tree<>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

std::vector<FT> bounding_box_dimensions(Point_vector& points) {
  std::vector<FT> lower, upper, difference;
  for (auto x: points[0]) {
    lower.push_back(x);
    upper.push_back(x);
  }
  for (auto p: points)
    for (unsigned i = 0; i < p.size(); i++) {
      if (p[i] < lower[i])
        lower[i] = p[i];
      if (p[i] > upper[i])
        upper[i] = p[i];
    }
  for (unsigned i = 0; i < lower.size(); i++)
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
void write_coxeter_mesh(Point_vector& W, VSMap& vs_map, Matrix& root_t, std::string file_name = "coxeter.mesh")
{
  short d = W[0].size();
  if (d > 3);
  
  std::ofstream ofs (file_name, std::ofstream::out);
  if (d <= 2)
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
  else
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
  
  // int num_vertices = (d+1)*W.size() + vs_map.size(), num_edges = 0, num_triangles = 0, num_tetrahedra = 0;
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

  ofs << "Vertices\n" << vs_map.size() << "\n";

  W.clear();
  for (auto m: vs_map) {
    FT denom = m.first[0];
    Eigen::VectorXd b(d);
    for (int i = 0; i < d; i++) {
      b(i) = m.first[i+1]/denom;
    }
    Eigen::SimplicialLDLT<Matrix, Eigen::Upper> chol(root_t);
    Eigen::VectorXd x = chol.solve(b);
    if(chol.info()!=Eigen::Success) {
      std::cout << "solving failed\n";
    }
    std::vector<FT> np;
    for (int i = 0; i < d; i++)
      np.push_back(x(i));
    W.push_back(Point_d(np));
  }
  std::map<int,std::vector<int>> sv_map;
  int j = 1;
  for (auto m: vs_map) {
    for (auto s: m.second) {
      auto find_it = sv_map.find(s);
      if (find_it == sv_map.end())
        sv_map.emplace(s, std::vector<int>(1,j));
      else
        find_it->second.push_back(j);
    }
    j++;
  }

  
  // FT p_prop = 0.001;
  int p_col = 208;
  std::vector<FT> bbox_dimensions = bounding_box_dimensions(W);
  std::cout << bbox_dimensions << "\n";
  for (auto p: W) {
    for (auto coord = p.cartesian_begin(); coord != p.cartesian_end() && coord != p.cartesian_begin()+3 ; ++coord) 
      ofs << *coord << " "; 
    ofs << "508\n";
    // for (int i = 0; i < d; i++) {
    //   int j = 0;
    //   for (auto coord = p.cartesian_begin(); coord != p.cartesian_end() && coord != p.cartesian_begin()+3 ; ++coord, j++)
    //     if (j == i)
    //       ofs << *coord + bbox_dimensions[i]*p_prop << " ";
    //     else
    //       ofs << *coord << " ";
    //   ofs << "108\n";
    // }
  }
  // num_edges = ((d+1)*d/2)*W.size()/2;
  // num_triangles = ((d+1)*d*(d-1)/6)*W.size();
  // num_tetrahedra = W.size();
  
  // ofs << "Edges " << num_edges << "\n";
  // for (unsigned i = 0; i < W.size(); i++) {
  //   ofs << (d+1)*i+1 << " " << (d+1)*i+2 << " " << p_col << "\n";
  //   ofs << (d+1)*i+1 << " " << (d+1)*i+3 << " " << p_col << "\n";
  //   ofs << (d+1)*i+2 << " " << (d+1)*i+3 << " " << p_col << "\n";
  //   if (d == 3) {
  //     ofs << (d+1)*i+1 << " " << (d+1)*i+4 << " " << p_col << "\n";
  //     ofs << (d+1)*i+2 << " " << (d+1)*i+4 << " " << p_col << "\n";
  //     ofs << (d+1)*i+3 << " " << (d+1)*i+4 << " " << p_col << "\n";
  //   }
  // }
  if (d == 2) {
    ofs << "Triangles " << sv_map.size() << "\n";
    for (auto m: sv_map) {
      for (auto i: m.second)
        ofs << i << " ";
      ofs << p_col << "\n";
    }
  }
  if (d == 3) {
    ofs << "Tetrahedra " << sv_map.size() << "\n";  
    for (auto m: sv_map) {
      for (auto i: m.second)
        ofs << i << " ";
      ofs << p_col << "\n";
    }
  }
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
  for (int i = 0; i < d; i++) {
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

  Coxeter_system cs('A', d-1);
  cs.emplace('A', 1);

  Matrix root_t;// = cs.root_t_;
  // Point_vector rc_point_vector = root_coordinates_range(point_vector, root_t);

  // The first fill of a map: simplex coordinates -> points
  std::cout << "Point cartesian coordinates: " << point_vector[0] << std::endl
    // << "Simple root root coordinates: " << root_coordinates(point_vector[0], root_t, d) << std::endl
            << "Alcove coordinates: " << cs.alcove_coordinates(point_vector[0], 1) << std::endl;
  // std::cout << std::endl;
  SPMap sp_map;
  for (auto p_it = point_vector.begin(); p_it != point_vector.end(); ++p_it) {
    Simplex_id s_id = cs.alcove_coordinates(*p_it, 1); 
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

  // // small test
  // Simplex_id p1 = {4,12,10,-8};
  // std::cout << "Non-reduced: " << p1 << ", reduced: " << reduced_id(p1) << ".\n";

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
  // add_vertices_to_map(si_map.begin(), vs_map);
  for (auto si_it = si_map.begin(); si_it != si_map.end(); ++si_it) {
    std::vector<Vertex_id> vertices = cs.vertices_of_alcove(si_it->first->first);
    for (Vertex_id v: vertices) {
      auto find_it = vs_map.find(v);
      if (find_it == vs_map.end())
        vs_map.emplace(v, std::vector<int>(1, si_it->second));
      else
        find_it->second.push_back(si_it->second);    
    }
  }
  std::cout << "VSMap composition:\n";
  for (auto m: vs_map) {
    std::cout << m.first << ": " << m.second << ".\n";
  }

  // simplex tree construction
  // Simplex_tree st;
  // for (auto m: vs_map) {
  //   st.insert_simplex_and_subfaces(m.second, m.first[0]);
  // }
  // std::cout << st;

  // Persistent_cohomology pcoh(st);
  // // initializes the coefficient field for homology
  // pcoh.init_coefficients(11);

  // pcoh.compute_persistent_cohomology(-0.1);
  // pcoh.output_diagram();

  write_coxeter_mesh(point_vector, vs_map, root_t, "coxeter.mesh");
}
