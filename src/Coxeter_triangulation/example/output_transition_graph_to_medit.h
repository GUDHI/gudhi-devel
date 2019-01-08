#ifndef OUTPUT_TRANSITION_GRAPH_TO_MEDIT_H_
#define OUTPUT_TRANSITION_GRAPH_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

template <class CellSet,
	  class Coxeter_triangulation_ds>
void output_transition_graph_to_medit(const CellSet& max_cells,
				      const Coxeter_triangulation_ds& ct,
				      std::string file_name = "graph")
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
  std::unordered_set<Cell_id> edges;
  for (auto c: max_cells) {
    ci_map.emplace(std::make_pair(c, ci_map.size()+1));
    Point_d barycenter(d);
    for (std::size_t i = 0; i < d; ++i)
      barycenter[i] = 0;
    std::unordered_set<Cell_id> v_cells;
    Coface_it cof_it(c, ct, c.dimension()+1);
    Coface_it cof_end;
    std::size_t num_vertices = 0;
    for (; cof_it != cof_end; ++cof_it) {
      std::size_t num_cof_cells = 0;
      for (auto f: ct.face_range(*cof_it, 1)) 
	if (max_cells.find(f) != max_cells.end())
	  num_cof_cells++;
      if (num_cof_cells > 1)
	edges.emplace(*cof_it);
    }
    Coface_it cof2_it(c, ct, d);
    for (; cof2_it != cof_end; ++cof2_it, ++num_vertices) {
      Eigen::VectorXd v = ct.barycenter(*cof2_it);
      for (std::size_t i = 0; i < d; ++i)
	barycenter[i] += v(i);
    }
    vertex_points.push_back(barycenter / num_vertices);
  }
  if (d == 2) { // TODO: Still same as the max_cells
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 215\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto e: edges) {
      for (auto f: ct.face_range(e, 1))
	if (max_cells.find(f) != max_cells.end()) {
	  ofs << ci_map.at(f) << " ";
	}
      ofs << "515" << std::endl;
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
    for (auto e: edges) {
      for (auto f: ct.face_range(e, 1))
	if (max_cells.find(f) != max_cells.end()) {
	  ofs << ci_map.at(f) << " ";
	}
      ofs << "515" << std::endl;
    }
  }
  ofs.close();
  ofs_bb.close();
}


#endif
