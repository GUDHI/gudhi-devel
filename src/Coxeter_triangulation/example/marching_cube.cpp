#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>

#include <gudhi/Points_off_io.h>
#include <gudhi/Coxeter_system.h>
#include <gudhi/Coxeter_complex.h>

#include <CGAL/Epick_d.h>

//#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cxx-prettyprint/prettyprint.hpp"

using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using FT = K::FT;
using Point_d = K::Point_d;
using Point_vector = std::vector< Point_d >;
using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

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
template <class VSMap,
          class Matrix>
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

void rec_test(std::vector<unsigned>& decomposition, Coxeter_system& cs, Point_vector& point_vector) {
  if (decomposition[0] == point_vector[0].size()) {
    // test if the number of vertices of 1 alcove is more than 10^6
    unsigned num_vertices = 1;
    for (auto d_it = decomposition.begin()+1; d_it != decomposition.end(); ++d_it) {
      num_vertices *= *d_it + 1;
      if (num_vertices > 1000000) {
        std::cout << "Configuration " << decomposition << ": too many vertices. Abandon.\n";
        return;
      }
    }
    std::cout << std::vector<unsigned>(decomposition.begin()+1, decomposition.end()) << std::endl;
    Coxeter_complex(point_vector, cs);
    return;
  }
  unsigned i = decomposition.back();
  if (decomposition.back() == 0)
    i = 1;
  for (; i <= point_vector[0].size() - decomposition[0]; ++i) {
    decomposition.push_back(i);
    decomposition[0] += i;
    cs.emplace_back('A', i);
    rec_test(decomposition, cs, point_vector);
    cs.pop_back();
    decomposition[0] -= i;
    decomposition.pop_back();
  }
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
  
  std::vector<unsigned> decomposition(1,0); // first coordinate is the sum
  decomposition.reserve(d+1);
  Coxeter_system cs;
  rec_test(decomposition, cs, point_vector);
}
