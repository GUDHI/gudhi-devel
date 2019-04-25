#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <fstream>
#include <cstdlib>
#include <unordered_set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/filesystem.hpp>

#include <gudhi_patches/CGAL/Epick_d.h>
#include <gudhi/random_point_generators.h> // construct_point
#include <gudhi/Coxeter_triangulation_ds.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include <gudhi/Kd_tree_search.h>
#include <gudhi/Points_off_io.h>
#include "output_max_cells_to_medit.h"
#include "output_allgowerschmidt.h"
#include "output_transition_graph_to_medit.h"

#include <CGAL/Search_traits_d.h>
#include <CGAL/Fuzzy_sphere.h>

#include <Eigen/Eigenvalues> 
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include "functions/sphere_S1_in_R2.h"
#include "functions/sphere_S2_in_R3.h"
#include "functions/chair_in_R3.h"


using Cell_id = Gudhi::Cell_id;
// using Point_d = Eigen::VectorXd;
using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Point_d = Kernel::Point_d;
using Point_range = std::vector<Point_d>;
using CT = Gudhi::Coxeter_triangulation_ds;

using Traits = CGAL::Search_traits_d<Kernel, CGAL::Dynamic_dimension_tag>;
using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;

std::size_t GUDHI_CT_NB_POINTS_FOR_PCA = 20;
Point_range points;
std::size_t snapshot_num = 0;

std::vector<Cell_id> ambient_edges_transv;
std::vector<Cell_id> ambient_edges_n_transv;

// template <>
void output_curve_to_medit() {
  std::vector<std::vector<double> > vertex_points;
  std::vector<std::vector<std::size_t> > edges, triangles, tetrahedra;
  std::vector<std::size_t> trn_mask, edge_mask;
  std::vector<bool> neighborhood_mask(points.size(), false);

  std::ofstream ofs ("curve" + std::to_string(snapshot_num) + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb ("curve" + std::to_string(snapshot_num++) + ".bb", std::ofstream::out);

  double r = 0.005/ambient_edges_transv.begin()->level();
  std::vector<double>
    p0 = {r/std::sqrt(2), -r/2},
    p1 = {-r/std::sqrt(2), -r/2},
    p2 = {0, r};  
  for (std::size_t i = 0; i < points.size(); ++i) {
    Point_d& p = points[i];
    std::size_t ind = vertex_points.size()+1;
    triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
    trn_mask.push_back(1); trn_mask.push_back(1); trn_mask.push_back(1);
    vertex_points.push_back(std::vector<double>({p[0]+p0[0], p[1] + p0[1]}));
    vertex_points.push_back(std::vector<double>({p[0]+p1[0], p[1] + p1[1]}));
    vertex_points.push_back(std::vector<double>({p[0]+p2[0], p[1] + p2[1]}));
  }
  CT ct(2);
  for (auto c: ambient_edges_n_transv) {
    std::size_t ind = vertex_points.size()+1;
    edges.push_back(std::vector<std::size_t>({ind, ind+1}));
    edge_mask.push_back(2);
    for (auto v: ct.vertex_range(c)) {
      Eigen::VectorXd v_ = ct.cartesian_coordinates(v);
      vertex_points.push_back(std::vector<double>({v_(0), v_(1)}));
    }
  }
  for (auto c: ambient_edges_transv) {
    std::size_t ind = vertex_points.size()+1;
    edges.push_back(std::vector<std::size_t>({ind, ind+1}));
    edge_mask.push_back(3);
    for (auto v: ct.vertex_range(c)) {
      Eigen::VectorXd v_ = ct.cartesian_coordinates(v);
      vertex_points.push_back(std::vector<double>({v_(0), v_(1)}));
    }
    edges.push_back(std::vector<std::size_t>({ind+2, ind+3}));
    edge_mask.push_back(4);
    Coface_it c_it(c, ct, 2), c_end;
    for (; c_it != c_end; ++c_it) {
      Eigen::VectorXd v_ = ct.barycenter(*c_it);
      vertex_points.push_back(std::vector<double>({v_(0), v_(1)}));
    }
  }
  // std::size_t ind = vertex_points.size()+1;
  // if (vertices.size() == 2) {
  //   edges.push_back(std::vector<std::size_t>({ind, ind+1}));
  //   for (auto v: vertices)
  //     vertex_points.push_back(std::vector<double>({v(0), v(1), v(2)}));
  //   nn_it = kns_range.begin();
  //   Point_d& p = points[nn_it->first];
  //   Eigen::VectorXd
  //     u0 = eigenvectors.col(1),
  //     u1 = eigenvectors.col(2);
  //   ind = vertex_points.size()+1;
  //   triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
  //   triangles.push_back(std::vector<std::size_t>({ind, ind+2, ind+3}));
  //   triangles.push_back(std::vector<std::size_t>({ind, ind+3, ind+4}));
  //   triangles.push_back(std::vector<std::size_t>({ind, ind+4, ind+1}));
  //   trn_mask.push_back(3); trn_mask.push_back(3); trn_mask.push_back(3); trn_mask.push_back(3);
  //   vertex_points.push_back(std::vector<double>({p[0], p[1], p[2]}));
  //   vertex_points.push_back(std::vector<double>({p[0] + u0(0) + u1(0),
  // 						 p[1] + u0(1) + u1(1),
  // 						 p[2] + u0(2) + u1(2)}));
  //   vertex_points.push_back(std::vector<double>({p[0] + u0(0) - u1(0),
  // 						 p[1] + u0(1) - u1(1),
  // 						 p[2] + u0(2) - u1(2)}));
  //   vertex_points.push_back(std::vector<double>({p[0] - u0(0) - u1(0),
  // 						 p[1] - u0(1) - u1(1),
  // 						 p[2] - u0(2) - u1(2)}));
  //   vertex_points.push_back(std::vector<double>({p[0] - u0(0) + u1(0),
  // 						 p[1] - u0(1) + u1(1),
  // 						 p[2] - u0(2) + u1(2)}));
  // }
  // if (vertices.size() == 3) {
  //   triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
  //   trn_mask.push_back(4);
  //   for (auto v: vertices)
  //     vertex_points.push_back(std::vector<double>({v(0), v(1), v(2)}));
  //   nn_it = kns_range.begin();
  //   Point_d& p = points[nn_it->first];
  //   Eigen::VectorXd u0 = eigenvectors.col(2);
  //   ind = vertex_points.size()+1;
  //   edges.push_back(std::vector<std::size_t>({ind, ind+1}));
  //   vertex_points.push_back(std::vector<double>({p[0] - u0(0),
  // 						 p[1] - u0(1),
  // 						 p[2] - u0(2)}));
  //   vertex_points.push_back(std::vector<double>({p[0] + u0(0),
  // 						 p[1] + u0(1),
  // 						 p[2] + u0(2)}));
  // }  
  ofs << "MeshVersionFormatted 1\nDimension 2\n";
  ofs_bb << "2 1 ";
  ofs << "Vertices\n" << vertex_points.size() << "\n";
  for (auto p: vertex_points) {
    ofs << p[0] << " " << p[1] << " 215\n";
  }
  ofs << "Edges " << edges.size() << "\n";
  auto em_it = edge_mask.begin();
  for (auto s: edges) {
    for (auto v: s) {
      ofs << v << " ";
    }
    ofs << *em_it++ << std::endl;
  }
  ofs << "Triangles " << triangles.size() << "\n";
  ofs_bb << triangles.size()+tetrahedra.size() << " 1\n";
  auto m_it = trn_mask.begin();
  // auto f_it = filtrations.begin();
  for (auto s: triangles) {
    for (auto v: s) {
      ofs << v << " ";
    }
    ofs << *m_it++ << std::endl;
    // ofs_bb << *f_it++ << "\n";
  }
  ofs.close();
  ofs_bb.close();
}


template <class KNSRange,
	  class VertexRange,
	  class Eigenvectors>
void output_config_to_medit(KNSRange& kns_range,
			    VertexRange& vertices,
			    Eigenvectors& eigenvectors) {
  std::vector<std::vector<double> > vertex_points;
  std::vector<std::vector<std::size_t> > edges, triangles, tetrahedra;
  std::vector<std::size_t> trn_mask, tet_mask;
  std::vector<bool> neighborhood_mask(points.size(), false);

  std::ofstream ofs ("non_insertion" + std::to_string(snapshot_num) + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb ("non_insertion" + std::to_string(snapshot_num++) + ".bb", std::ofstream::out);

  double r = 0.005;
  std::vector<double>
    p0 = {r, 0, -r/std::sqrt(2)},
    p1 = {-r, 0, -r/std::sqrt(2)},
    p2 = {0, r, r/std::sqrt(2)},
    p3 = {0, -r, r/std::sqrt(2)};

  std::size_t nb_pts_pca = std::min(static_cast<std::size_t>(std::pow(GUDHI_CT_NB_POINTS_FOR_PCA, 2)), points.size());
  // std::size_t nb_pts_pca = std::min(static_cast<std::size_t>(std::pow(90, 2)), points.size());
  auto nn_it = kns_range.begin();
  for (std::size_t i = 0; i < nb_pts_pca && nn_it != kns_range.end(); ++i, ++nn_it)
    neighborhood_mask[nn_it->first] = true;
  
  for (std::size_t i = 0; i < points.size(); ++i) {
    Point_d& p = points[i];
    std::size_t ind = vertex_points.size()+1;
    triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
    triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+3}));
    triangles.push_back(std::vector<std::size_t>({ind, ind+2, ind+3}));
    triangles.push_back(std::vector<std::size_t>({ind+1, ind+2, ind+3}));
    if (neighborhood_mask[i]) {
      trn_mask.push_back(1); trn_mask.push_back(1); trn_mask.push_back(1); trn_mask.push_back(1);
    }
    else {
      trn_mask.push_back(2); trn_mask.push_back(2); trn_mask.push_back(2); trn_mask.push_back(2);
    }
    vertex_points.push_back(std::vector<double>({p[0]+p0[0], p[1] + p0[1], p[2] + p0[2]}));
    vertex_points.push_back(std::vector<double>({p[0]+p1[0], p[1] + p1[1], p[2] + p1[2]}));
    vertex_points.push_back(std::vector<double>({p[0]+p2[0], p[1] + p2[1], p[2] + p2[2]}));
    vertex_points.push_back(std::vector<double>({p[0]+p3[0], p[1] + p3[1], p[2] + p3[2]}));
  }
  std::size_t ind = vertex_points.size()+1;
  if (vertices.size() == 2) {
    edges.push_back(std::vector<std::size_t>({ind, ind+1}));
    for (auto v: vertices)
      vertex_points.push_back(std::vector<double>({v(0), v(1), v(2)}));
    nn_it = kns_range.begin();
    Point_d& p = points[nn_it->first];
    Eigen::VectorXd
      u0 = eigenvectors.col(1),
      u1 = eigenvectors.col(2);
    ind = vertex_points.size()+1;
    triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
    triangles.push_back(std::vector<std::size_t>({ind, ind+2, ind+3}));
    triangles.push_back(std::vector<std::size_t>({ind, ind+3, ind+4}));
    triangles.push_back(std::vector<std::size_t>({ind, ind+4, ind+1}));
    trn_mask.push_back(3); trn_mask.push_back(3); trn_mask.push_back(3); trn_mask.push_back(3);
    vertex_points.push_back(std::vector<double>({p[0], p[1], p[2]}));
    vertex_points.push_back(std::vector<double>({p[0] + u0(0) + u1(0),
						 p[1] + u0(1) + u1(1),
						 p[2] + u0(2) + u1(2)}));
    vertex_points.push_back(std::vector<double>({p[0] + u0(0) - u1(0),
						 p[1] + u0(1) - u1(1),
						 p[2] + u0(2) - u1(2)}));
    vertex_points.push_back(std::vector<double>({p[0] - u0(0) - u1(0),
						 p[1] - u0(1) - u1(1),
						 p[2] - u0(2) - u1(2)}));
    vertex_points.push_back(std::vector<double>({p[0] - u0(0) + u1(0),
						 p[1] - u0(1) + u1(1),
						 p[2] - u0(2) + u1(2)}));
  }
  if (vertices.size() == 3) {
    triangles.push_back(std::vector<std::size_t>({ind, ind+1, ind+2}));
    trn_mask.push_back(4);
    for (auto v: vertices)
      vertex_points.push_back(std::vector<double>({v(0), v(1), v(2)}));
    nn_it = kns_range.begin();
    Point_d& p = points[nn_it->first];
    Eigen::VectorXd u0 = eigenvectors.col(2);
    ind = vertex_points.size()+1;
    edges.push_back(std::vector<std::size_t>({ind, ind+1}));
    vertex_points.push_back(std::vector<double>({p[0] - u0(0),
						 p[1] - u0(1),
						 p[2] - u0(2)}));
    vertex_points.push_back(std::vector<double>({p[0] + u0(0),
						 p[1] + u0(1),
						 p[2] + u0(2)}));
  }  
  ofs << "MeshVersionFormatted 1\nDimension 3\n";
  ofs_bb << "3 1 ";
  ofs << "Vertices\n" << vertex_points.size() << "\n";
  for (auto p: vertex_points) {
    ofs << p[0] << " " << p[1] << " " << p[2] << " 215\n";
  }
  ofs << "Edges " << edges.size() << "\n";
  for (auto s: edges) {
    for (auto v: s) {
      ofs << v << " ";
    }
    ofs << "515" << std::endl;
  }
  ofs << "Triangles " << triangles.size() << "\n";
  ofs_bb << triangles.size()+tetrahedra.size() << " 1\n";
  auto m_it = trn_mask.begin();
  // auto f_it = filtrations.begin();
  for (auto s: triangles) {
    for (auto v: s) {
      ofs << v << " ";
    }
    ofs << *m_it++ << std::endl;
    // ofs_bb << *f_it++ << "\n";
  }
  ofs << "Tetrahedra " << tetrahedra.size() << "\n";
  for (auto s: tetrahedra) {
    for (auto v: s) {
      ofs << v << " ";
    }
    // ofs << "545\n";
    // ofs_bb << *f_it++ << "\n";
  }
  ofs.close();
  ofs_bb.close();
}


template <class Kd_tree_search>
bool intersects(const Cell_id& c,
		const Kd_tree_search& kd_tree,
		const CT& ct) {
  std::size_t amb_d = ct.dimension();
  std::size_t intr_d = amb_d - c.dimension();
  std::size_t cod_d = amb_d - intr_d;
  Kernel k;
  auto pt = k.construct_point_d_object();
  auto coord = k.compute_coordinate_d_object();
  auto scalar = k.scalar_product_d_object();
  auto sq_dist = k.squared_distance_d_object();
  std::size_t nb_pts_pca = std::min(static_cast<std::size_t>(std::pow(GUDHI_CT_NB_POINTS_FOR_PCA, intr_d)), points.size());  
  Eigen::MatrixXd pt_matrix(nb_pts_pca, amb_d);
  Eigen::VectorXd b = ct.barycenter(c);
  Point_d b_point = pt(amb_d, b.data(), b.data()+amb_d);

  // KNS with Fuzzy_sphere
  Eigen::VectorXd v0 = ct.cartesian_coordinates(*ct.vertex_range(c).begin());
  double rad = std::sqrt(sq_dist(b_point, pt(amb_d, v0.data(), v0.data()+amb_d)));
  double fs_err = 1e-5;
  Fuzzy_sphere fs(b_point, rad, fs_err);
  std::vector<Point_d> kns_points;
  kd_tree.search(std::back_inserter<std::vector<Point_d> >(kns_points), fs);
  
  auto kns_range = kd_tree.k_nearest_neighbors(b_point, nb_pts_pca, false);
  auto nn_it = kns_range.begin();
  Point_d closest_point = points[nn_it->first];
  std::size_t curr_row = 0;
  for (; nn_it != kns_range.end(), curr_row < (std::size_t)pt_matrix.rows(); ++nn_it, ++curr_row)
    for (std::size_t i = 0; i < amb_d; ++i)
      pt_matrix(curr_row, i) = CGAL::to_double(coord(points[nn_it->first], i));
  Eigen::MatrixXd centered = pt_matrix.rowwise() - pt_matrix.colwise().mean();
  Eigen::MatrixXd cov = centered.adjoint() * centered;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

  Eigen::VectorXd free_term(cod_d);
  for (std::size_t i = 0; i < cod_d; ++i) {
    Point_d ui = pt(amb_d,
		    eig.eigenvectors().col(i).data(),
		    eig.eigenvectors().col(i).data() + amb_d);
    // free_term(i) = scalar(ui, closest_point);
    // free_term(i) = scalar(ui, b_point);
    double cheat_lambda = 0, cheat_colambda = 1 - cheat_lambda;
    free_term(i) = cheat_colambda*scalar(ui, closest_point) + cheat_lambda*scalar(ui, b_point);
  }
    
  Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
  for (std::size_t i = 0; i < cod_d + 1; ++i)
    matrix(0, i) = 1;
  std::size_t j = 0;
  for (auto v: ct.vertex_range(c)) {
    Eigen::VectorXd v_vec = ct.cartesian_coordinates(v);
    for (std::size_t i = 0; i < cod_d; ++i)
      matrix(i+1, j) = v_vec.dot(eig.eigenvectors().col(i)) - free_term(i);
    j++;
  }
  Eigen::VectorXd z(cod_d + 1);
  z(0) = 1;
  for (std::size_t i = 1; i < cod_d + 1; ++i)
    z(i) = 0;
  Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
  double toler = 1e-4;
  for (std::size_t i = 0; i < cod_d + 1; ++i)
    if (lambda(i) < -toler || lambda(i) > 1+toler) {
      // std::vector<Eigen::VectorXd> vertices;
      // for (auto v: ct.vertex_range(c))
      // 	vertices.push_back(ct.cartesian_coordinates(v));
      // output_config_to_medit(kns_range, vertices, eig.eigenvectors());
      ambient_edges_n_transv.push_back(c);
      return false;
    }
  ambient_edges_transv.push_back(c);  
  return true;
}


template <class Point,
	  class Kd_tree_search>
void compute_complex(const Point& seed_point,
		     std::size_t intr_d,
		     double level,
		     std::unordered_set<Cell_id>& max_cells,
		     const Kd_tree_search& kd_tree,
		     bool output_to_medit = true,
		     std::string file_name_prefix = "reconstruction") {
  std::size_t amb_d = seed_point.size();
  std::size_t cod_d = amb_d - intr_d;
  CT ct(amb_d);
  std::unordered_set<Cell_id> facet_cells;
  std::unordered_map<Cell_id, std::size_t> order_map;
    
  std::queue<Cell_id> queue;
  Cell_id c = ct.locate_point(seed_point, level);
  for (auto f: ct.face_range(c, cod_d))
    if (intersects(f, kd_tree, ct) && max_cells.emplace(f).second && order_map.emplace(std::make_pair(f, order_map.size())).second)
      queue.emplace(f);

  // std::size_t snapshot_num = 0;
  while (!queue.empty()) {
    Cell_id s = queue.front();
    queue.pop();
    // Graph_node node = boost::add_vertex(graph);
    // nc_map.emplace(std::make_pair(node, s));
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(s, ct, cod_d+1);
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
    for (; cof_it != cof_end; ++cof_it)
      if (facet_cells.emplace(*cof_it).second)
	for (auto f: ct.face_range(*cof_it, cod_d))
	  if (intersects(f, kd_tree, ct) && max_cells.emplace(f).second && order_map.emplace(std::make_pair(f, order_map.size())).second)
	    queue.emplace(f);
    // std::cout << "queue.size() = " << queue.size() << "\n";
    // output_max_cells_to_medit(max_cells, ct, file_name_prefix + std::to_string(snapshot_num++));
  }
  std::cout << "#max_cells = " << max_cells.size() << "\n";
  std::cout << "#facet_cells = " << facet_cells.size() << "\n";

  // std::unordered_set<Cell_id> half_set;
  // for (auto c: max_cells)
  //   if (c.value(0) == 0)
  //     half_set.emplace(c);
  // output_max_cells_to_medit(half_set, ct, order_map, file_name_prefix);
  // output_transition_graph_to_medit(half_set, ct, file_name_prefix + "_tg");
  if (output_to_medit)
    output_max_cells_to_medit(max_cells, ct, order_map, cod_d, file_name_prefix);
  // output_allgowerschmidt_to_medit(max_cells, ct, order_map, file_name_prefix + "_as");
  output_transition_graph_to_medit(max_cells, ct, file_name_prefix + "_tg");
  output_curve_to_medit();
}
  

int main(int argc, char * const argv[]) {
  Kernel k;
  
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0]
	      << " path_to_point_file intrinsic_dimension level\n";
    return 0;
  }

  boost::filesystem::path path(argv[1]);
  std::size_t intr_d = atoi(argv[2]);
  double level = atof(argv[3]);

  std::cout << "File: " << path.string() << "\nIntrinsic dimension: " << intr_d << "\nLevel: " << level << "\n";
  Gudhi::Points_off_reader<Point_d> off_reader(path.string());
  if (!off_reader.is_valid()) {
    std::cerr << "Coxeter triangulation - Unable to read file " << path.string() << "\n";
    exit(-1);  // ----- >>
  }
  points = off_reader.get_point_cloud();
  std::cout << "Number of points: " << points.size() << "\n";
  if (points.empty()) {
    std::cerr << "Coxeter triangulation - No points found\n";
    exit(-1);  // ----- >>
  }
  Gudhi::spatial_searching::Kd_tree_search<Kernel, Point_range> kd_tree(points);
  std::cout << "Depth of the kd tree: " << kd_tree.tree_depth() << "\n";

  std::unordered_set<Cell_id> max_cells;
  std::string file_name = path.filename().stem().string();
  std::cout << "File name: " << file_name << "\n";
  compute_complex(*points.begin(), intr_d, level, max_cells, kd_tree, true, file_name);
}
