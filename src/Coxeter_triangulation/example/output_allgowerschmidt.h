#ifndef OUTPUT_ALLGOWERSCHMIDT_TO_MEDIT_H_
#define OUTPUT_ALLGOWERSCHMIDT_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

template <class CellSet,
	  class Coxeter_triangulation_ds,
	  class OrderMap>
void output_allgowerschmidt_to_medit(const CellSet& max_cells,
				     const Coxeter_triangulation_ds& ct,
				     const OrderMap& order_map,
				     std::string file_name = "reconstruction")
{
  using Cell_id = typename CellSet::value_type;
  using Point_d = Eigen::VectorXd;
  using Coface_it = Gudhi::Coxeter_triangulation_ds::Coface_iterator;
  
  if (max_cells.empty())
    return;
  // determine dimension from VP_map
  unsigned d = ct.dimension();

  std::ofstream ofs (file_name + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb (file_name + ".bb", std::ofstream::out);

  std::vector<Point_d> vertex_points;
  std::unordered_map<Cell_id, std::size_t> ci_map;
  std::unordered_set<Cell_id> edges, polygons;
  for (auto c: max_cells) {
    vertex_points.emplace_back(ct.barycenter(c));
    perturb_vertex(vertex_points.back(), 0.00001);
    ci_map.emplace(std::make_pair(c, vertex_points.size()));
    Coface_it edge_it(c, ct, c.dimension()+1);
    Coface_it cof_end;
    for (; edge_it != cof_end; ++edge_it)
      edges.emplace(*edge_it);
    if (d >= 3 && c.dimension()+2 <= d) {
      Coface_it pol_it(c, ct, c.dimension()+2);
      for (; pol_it != cof_end; ++pol_it)
	polygons.emplace(*pol_it); 
    }
  }
  if (d == 2) {
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 215\n";
    }
    ofs << "Edges " << max_cells.size() << "\n";
    for (auto e: edges) {
      for (auto v: ct.face_range(e, 1)) {
    	ofs << ci_map.at(v) << " ";
      }
      ofs << "515" << std::endl;
    }
  } // d == 3
  else {
    std::unordered_set<Cell_id> edges;
    std::vector<std::vector<std::size_t> > triangles;
    std::vector<double> filtrations;
    std::vector<std::size_t> mask;
    if (max_cells.begin()->dimension() == 2) {
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
      ofs_bb << "3 1 ";
      ofs << "Vertices\n" << vertex_points.size() << "\n";
      for (auto p: vertex_points) {
	ofs << p[0] << " " << p[1] << " 215\n";
      }
      ofs << "Edges " << edges.size() << "\n";
      for (auto e: edges) {
	for (auto v: ct.face_range(e, 1)) {
	  ofs << ci_map.at(v) << " ";
	}
	ofs << "515" << std::endl;
      }
    }
    else { // max_cells.begin()->dimension() == 1
      for (auto c: polygons) {
	// std::size_t mask_val = 515;
        Point_d barycenter(3);
        barycenter[0] = 0;
        barycenter[1] = 0;
        barycenter[2] = 0;
        std::unordered_set<Cell_id> v_cells;
	std::size_t num_vertices = 0;
	for (Cell_id v: ct.face_range(c, 1)) 
	  if (max_cells.find(v) != max_cells.end()) {
	    int ci_value = ci_map.at(v);
	    barycenter[0] += vertex_points[ci_value-1][0];
	    barycenter[1] += vertex_points[ci_value-1][1];
	    barycenter[2] += vertex_points[ci_value-1][2];
	    num_vertices++;
	  }
        vertex_points.push_back(barycenter / num_vertices);
	std::size_t bc_index = vertex_points.size();
	for (Cell_id e: ct.face_range(c, 2)) {
	  Point_d barycenter(3);
	  barycenter[0] = 0;
	  barycenter[1] = 0;
	  barycenter[2] = 0;
	  std::unordered_set<Cell_id> v_cells;
	  for (Cell_id v: ct.face_range(e, 1))
	    if (max_cells.find(v) != max_cells.end()) {
	      int ci_value = ci_map.at(v);
	      barycenter[0] += vertex_points[ci_value-1][0];
	      barycenter[1] += vertex_points[ci_value-1][1];
	      barycenter[2] += vertex_points[ci_value-1][2];
	      std::vector<std::size_t> triangle {bc_index, vertex_points.size()+1, ci_map.at(v)};
	      triangles.push_back(triangle);
	      filtrations.push_back(order_map.at(v));
	      mask.push_back(order_map.at(v));
	    }
	  vertex_points.push_back(barycenter / 2);
	}
      }
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
      ofs_bb << "3 1 ";
      ofs << "Vertices\n" << vertex_points.size() << "\n";
      for (auto p: vertex_points) {
        ofs << p[0] << " " << p[1] << " " << p[2] << " 215\n";
      }
      ofs << "Edges " << edges.size() << "\n";
      for (auto e: edges) {
	for (auto v: ct.vertex_range(e)) {
	  ofs << ci_map.at(v) << " ";
	}
	ofs << "515" << std::endl;
      }
      ofs << "Triangles " << triangles.size() << "\n";
      ofs_bb << triangles.size() << " 1\n";
      auto m_it = mask.begin();
      auto f_it = filtrations.begin();
      for (auto s: triangles) {
        for (auto v: s)
	  ofs << v << " ";
        ofs << *m_it++ << std::endl;
        ofs_bb << *f_it++ << "\n";
      }
    }
  }
  ofs.close();
  ofs_bb.close();
}


#endif
