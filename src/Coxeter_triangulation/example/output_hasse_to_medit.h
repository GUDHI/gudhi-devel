#ifndef OUTPUT_HASSE_TO_MEDIT_H_
#define OUTPUT_HASSE_TO_MEDIT_H_

#include <string>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include <gudhi/Coxeter_triangulation/Cell_id.h>
#include <gudhi/Coxeter_triangulation_ds.h>
#include <gudhi/Hasse_diagram.h>

bool allgow_switch = false;

using Cell_id = Gudhi::Cell_id;
using Point_d_eigen = Eigen::VectorXd;
using Coface_it = Gudhi::Coxeter_triangulation_ds::Coface_iterator;
using Face_it = Gudhi::Coxeter_triangulation_ds::Face_iterator;

using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
using Hasse_boundary = std::vector<std::pair<Hasse_cell*, int> >;
// The vertices are always smaller than the non-vertices
struct Hasse_cell_comparator {
  bool operator() (Hasse_cell* l_it, Hasse_cell* r_it) {
    if (l_it->get_dimension() == 0)
      if (r_it->get_dimension() == 0)
	return l_it < r_it;
      else
	return true;
    else if (r_it->get_dimension() == 0)
      return false;
    else
      return l_it->get_boundary() < r_it->get_boundary();
  }
};
using Hasse_diagram = std::set<Hasse_cell*, Hasse_cell_comparator>;

struct Cell_id_comparator {
  bool operator()(const Cell_id& lhs, const Cell_id& rhs) {
    if (lhs.size() < rhs.size())
      return true;
    else if (lhs.size() > rhs.size())
      return false;
    for (std::size_t k = 0; k < lhs.size(); ++k) {
      if (lhs.value(k) < rhs.value(k))
	return true;
      else if (lhs.value(k) > rhs.value(k))
	return false;
      else if (!lhs.mask(k) && rhs.mask(k))
	return true;
      else if (lhs.mask(k) && !rhs.mask(k))
	return false;
    }
    return false;
  }
};
using VC_map = std::map<Cell_id, Hasse_cell*, Cell_id_comparator>;
using VP_map = std::map<Hasse_cell*, Point_d_eigen>;


template <class Point>
void perturb_vertex(Point& vertex, double rad) {
  int d = vertex.size();
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef typename Kernel::Point_d Point_d;
  CGAL::Random_points_on_sphere_d<Point_d> rp(d+1, rad);
  for (int i = 0; i < d; ++i)
    vertex[i] += (*rp)[i];
}


template <class VPMap,
	  class HasseDiagram>
void output_hasse_to_medit(const HasseDiagram& hasse_diagram,
			   const VPMap& vp_map,
			   bool barycentric,
			   std::string file_name = "hasse_diagram_output")
{
  using Hasse_cell_ptr = typename VPMap::key_type;
  using Point_t = typename VPMap::mapped_type;

  if (vp_map.empty())
    return;
  // determine dimension from VP_map
  unsigned d = vp_map.begin()->second.size();

  std::ofstream ofs (file_name + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb (file_name + ".bb", std::ofstream::out);

  std::vector<Hasse_cell_ptr> edges, polygons, polytopes;
  std::vector<std::vector<std::size_t> > triangles, tetrahedra; 
  std::vector<double> filtrations;
  std::vector<std::size_t> mask;
  std::vector<Point_t> vertex_points;
  std::map<Hasse_cell_ptr, std::size_t> ci_map;
  int index = 1;
  for (auto vp_pair: vp_map) {
    ci_map.emplace(vp_pair.first, index++);
    vertex_points.push_back(vp_pair.second);
    perturb_vertex(vertex_points.back(), 0.00001);
    // std::cout << "Inserted vertex " << vp_pair.first << ", index = " << ci_map.at(vp_pair.first) << "\n";
  }
  for (auto s: hasse_diagram)
    if (s->get_dimension() == 1)
      edges.push_back(s);
    else if (s->get_dimension() == 2)
      polygons.push_back(s);
    else if (s->get_dimension() == 3)
      polytopes.push_back(s);

  if (barycentric)
    for (auto e: edges) {
      Point_t edge_barycenter(d);
      for (std::size_t i = 0; i < d; ++i)
	edge_barycenter(i) = 0;
      for (auto v_pair: e->get_boundary())
	for (std::size_t i = 0; i < d; ++i)
	  edge_barycenter(i) += vp_map.at(v_pair.first)(i);
      vertex_points.emplace_back(edge_barycenter / 2);
      // perturb_vertex(vertex_points.back(), 0.00001);
      std::size_t e_index = ci_map.size()+1;
      ci_map.emplace(std::make_pair(e, e_index));
      // std::cout << "Inserted edge " << e << ", index = " << ci_map.at(e) << "\n";
    }
  std::size_t ref_size = vertex_points.size() - 1;
  for (auto h: polygons) {
    std::set<Hasse_cell_ptr> v_cells;
    for (auto e_pair: h->get_boundary())
      for (auto v_pair: e_pair.first->get_boundary())
	v_cells.emplace(v_pair.first);
    Point_t barycenter(d);
    for (std::size_t i = 0; i < d; ++i)
      barycenter(i) = 0;
    for (auto vc: v_cells) {
      for (std::size_t i = 0; i < d; ++i)
	barycenter(i) += vp_map.at(vc)(i);
    }
    vertex_points.emplace_back(barycenter / v_cells.size());
    // perturb_vertex(vertex_points.back(), 0.00001);
    std::size_t h_index = ci_map.size()+1;
    ci_map.emplace(std::make_pair(h, h_index));
    // std::cout << "Inserted polygon " << h << ", index = " << ci_map.at(h) << "\n";
    for (auto e_pair: h->get_boundary()) {
      if (barycentric) {	  
	for (auto v_pair: e_pair.first->get_boundary()) {
	  triangles.push_back(std::vector<std::size_t>({h_index, ci_map.at(e_pair.first), ci_map.at(v_pair.first)}));
	  filtrations.push_back(v_pair.first->get_filtration());
	  if (allgow_switch)
	    mask.push_back(std::round(v_pair.first->get_filtration()));
	  else
	    mask.push_back(std::round(h->get_filtration()));
	}
      }
      else { // not barycentric
	std::vector<std::size_t> triangle(1, h_index); 
	for (auto v_pair: e_pair.first->get_boundary())
	  triangle.push_back(ci_map.at(v_pair.first));
	triangles.push_back(triangle);
	filtrations.push_back(e_pair.first->get_filtration());
	mask.push_back(h_index - ref_size);
      }
    }
  }
  for (auto p: polytopes) {
    std::set<Hasse_cell_ptr> v_cells;
    for (auto h_pair: p->get_boundary())
      for (auto e_pair: h_pair.first->get_boundary())
	for (auto v_pair: e_pair.first->get_boundary())
	  v_cells.emplace(v_pair.first);
    Point_t barycenter(d);
    for (std::size_t i = 0; i < d; ++i)
      barycenter(i) = 0;
    for (auto vc: v_cells) {
      for (std::size_t i = 0; i < d; ++i)
	barycenter(i) += vp_map.at(vc)(i);
    }
    vertex_points.emplace_back(barycenter / v_cells.size());
    // perturb_vertex(vertex_points.back(), 0.00001);
    std::size_t p_index = vertex_points.size();
    ci_map.emplace(std::make_pair(p, p_index));
    for (auto h_pair: p->get_boundary())
      if (barycentric) {
	for (auto e_pair: h_pair.first->get_boundary())
	  for (auto v_pair: e_pair.first->get_boundary()) {
	    tetrahedra.push_back(std::vector<std::size_t>({p_index, ci_map.at(h_pair.first), ci_map.at(e_pair.first), ci_map.at(v_pair.first)}));
	    filtrations.push_back(v_pair.first->get_filtration());
	  }
      }
      else { // not barycentric // code below should fit
      }
  }
  
  if (d == 2) {
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 215\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto s: edges) {
      for (auto v: s->get_boundary()) {
	ofs << ci_map.at(v.first) << " ";
      }
      ofs << "515" << std::endl;
    }
    ofs << "Triangles " << triangles.size() << "\n";
    ofs_bb << triangles.size() << " 1\n";
    auto f_it = filtrations.begin();
    for (auto s: triangles) {
      for (auto v: s)
	ofs << v << " ";
      ofs << "516" << std::endl;
      ofs_bb << *f_it++ << "\n";
    }
  }
  else {
    ofs << "MeshVersionFormatted 1\nDimension 3\n";
    ofs_bb << "3 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " " << p[2] << " 215\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto s: edges) {
      for (auto v: s->get_boundary()) {
	ofs << ci_map.at(v.first) << " ";
      }
      ofs << "515" << std::endl;
    }
    ofs << "Triangles " << triangles.size() << "\n";
    ofs_bb << triangles.size()+tetrahedra.size() << " 1\n";
    auto m_it = mask.begin();
    auto f_it = filtrations.begin();
    for (auto s: triangles) {
      for (auto v: s) {
	ofs << v << " ";
      }
      ofs << *m_it++ << std::endl;
      ofs_bb << *f_it++ << "\n";
    }
    ofs << "Tetrahedra " << tetrahedra.size() << "\n";
    for (auto s: tetrahedra) {
      for (auto v: s) {
	ofs << v << " ";
      }
      ofs << "545\n";
      ofs_bb << *f_it++ << "\n";
    }

    // typedef CGAL::Epick_d<CGAL::Dimension_tag<3> > Kernel_3;
    // typedef typename Kernel_3::Point_d Point_3;
    // typedef CGAL::Delaunay_triangulation<Kernel_3> Delaunay_triangulation;
      
      // // constrain the outer edges, which are found by the convex hull
      // std::vector<Point_3> vertices;
      // std::vector<int> v_indices;
      // std::set<Hasse_cell_ptr> v_cells;
      // // std::cout << "3-cell " << p->get_filtration() << ":\n";
      // for (auto h_pair: p->get_boundary()) {
      // 	// std::cout << "2-cell " << h_pair.first->get_position() << " " << h_pair.first->get_filtration() << ":\n";
      // 	for (auto e_pair: h_pair.first->get_boundary()) {
      // 	  // std::cout << "Edge " << e_pair.first->get_position() << " " << e_pair.first->get_filtration() << ":\n"; 
      // 	  for (auto v_pair: e_pair.first->get_boundary()) {
      // 	    // std::cout << "Vertex " << v_pair.first->get_position() << " " << v_pair.first->get_filtration() << ":\n"; 
      // 	    v_cells.emplace(v_pair.first);
      // 	  }
      // 	}
      // 	// std::cout << "\n";
      // }
      // for (auto vc: v_cells) {
      // 	Point_t& b = vertex_points[ci_map.at(vc)-1];
      // 	vertices.push_back(Point_3(b[0], b[1], b[2]));
      // 	v_indices.push_back(ci_map.at(vc));
      // 	// std::cout << "Vertex " << vc->get_position() << " " << vc->get_filtration() << ":\n"; 
      // }
      // Delaunay_triangulation del(3);
      // index = 0;
      // std::map<typename Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
      // for (auto pt: vertices)
      // 	index_of_vertex.emplace(del.insert(pt), index++);
      // for (auto fc_it = del.full_cells_begin(); fc_it != del.full_cells_end(); ++fc_it) {
      // 	if (del.is_infinite(fc_it))
      // 	  continue;
      // 	std::vector<int> tetrahedron;
      // 	for (auto fv_it = fc_it->vertices_begin(); fv_it != fc_it->vertices_end(); ++fv_it)
      // 	  tetrahedron.push_back(v_indices[index_of_vertex[*fv_it]]);
      // 	tetrahedra.push_back(tetrahedron);
      // 	filtrations.push_back(p->get_filtration());
      // }
  }

  ofs.close();
  ofs_bb.close();
  }


#endif
