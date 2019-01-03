#ifndef OUTPUT_MAX_CELLS_TO_MEDIT_H_
#define OUTPUT_MAX_CELLS_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

template <class Point>
void perturb_vertex(Point& vertex, double rad) {
  int d = vertex.size();
  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
  typedef typename Kernel::Point_d Point_d;
  CGAL::Random_points_on_sphere_d<Point_d> rp(d+1, rad);
  for (int i = 0; i < d; ++i)
    vertex[i] += (*rp)[i];
}


template <class CellSet,
	  class Coxeter_triangulation_ds>
void output_max_cells_to_medit(const CellSet& max_cells,
			       const Coxeter_triangulation_ds& ct, 
			       std::string file_name = "reconstruction")
{
  using Cell_id = typename CellSet::value_type;
  using Point_d = Eigen::VectorXd;

  if (max_cells.empty())
    return;
  // determine dimension from VP_map
  unsigned d = ct.dimension();

  std::ofstream ofs (file_name + ".mesh", std::ofstream::out);
  std::ofstream ofs_bb (file_name + ".bb", std::ofstream::out);

  std::vector<Point_d> vertex_points;
  int index = 1;
  std::unordered_map<Cell_id, std::size_t> ci_map;
  for (auto c: max_cells) {
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(c, ct, d);
    Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
    for (; cof_it != cof_end; ++cof_it)
      if (ci_map.emplace(std::make_pair(*cof_it, index)).second) {
	vertex_points.emplace_back(ct.barycenter(*cof_it));
	perturb_vertex(vertex_points.back(), 0.00001);
	index++;
      }    
  }
  if (d == 2) {
    std::unordered_set<Cell_id> edges;
    std::vector<std::vector<int> > triangles;
    std::vector<double> filtrations;
    typedef CGAL::Epick_d<CGAL::Dimension_tag<2> > Kernel;
    typedef typename Kernel::Point_d Point_2;
    typedef CGAL::Delaunay_triangulation<Kernel> Delaunay_triangulation;
    for (auto c: max_cells) {
      // constrain the outer edges, which are found by the convex hull
      std::vector<Point_2> vertices;
      std::vector<int> v_indices;
      std::unordered_set<Cell_id> v_cells;
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(c, ct, d);
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
      for (; cof_it != cof_end; ++cof_it)
	v_cells.emplace(*cof_it);
      for (auto vc: v_cells) {
        v_indices.push_back(ci_map.at(vc));
        Point_d& b = vertex_points[v_indices.back()-1];
        vertices.push_back(Point_2(b[0], b[1]));
      }
      Delaunay_triangulation del(2);
      index = 0;
      std::map<typename Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
      for (auto p: vertices)
        index_of_vertex.emplace(del.insert(p), index++);
      for (auto fc_it = del.full_cells_begin(); fc_it != del.full_cells_end(); ++fc_it) {
        if (del.is_infinite(fc_it))
          continue;
        std::vector<int> triangle;
        for (auto fv_it = fc_it->vertices_begin(); fv_it != fc_it->vertices_end(); ++fv_it)
          triangle.push_back(v_indices[index_of_vertex[*fv_it]]);
        triangles.push_back(triangle);
        filtrations.push_back(0);
      }
    }
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 215\n";
    }
    ofs << "Edges " << edges.size() << "\n";
    for (auto e: edges) {
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(e, ct, d);
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
      for (; cof_it != cof_end; ++cof_it) {
    	ofs << ci_map.at(*cof_it) << " ";
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
  } // d == 3
  else {
    // std::vector<Hasse_cell_ptr> edges, polygons, polytopes;
    // std::vector<std::vector<int> > triangles, tetrahedra; 
    // std::vector<double> filtrations;
    // std::vector<int> mask;
    // // for (auto s: hasse_diagram)
    // //   if (s->get_dimension() == 1)
    // //     edges.push_back(s);
    // //   else if (s->get_dimension() == 2)
    // //     polygons.push_back(s);
    // //   else if (s->get_dimension() == 3)
    // //     polytopes.push_back(s);
    // typedef CGAL::Epick_d<CGAL::Dimension_tag<3> > Kernel_3;
    // typedef typename Kernel_3::Point_d Point_3;
    // typedef CGAL::Delaunay_triangulation<Kernel_3> Delaunay_triangulation;
    // for (auto h: polygons) {
    //   // int mask_val = 515;
    //   // Point_t barycenter(3);
    //   // barycenter[0] = 0;
    //   // barycenter[1] = 0;
    //   // barycenter[2] = 0;
    //   // std::set<Hasse_cell_ptr> v_cells;
    //   // for (auto e_pair: h->get_boundary()) {
    //   //   if (e_pair.first->get_coBoundary().size() > 2)
    //   //     mask_val = 517;
    //   //   for (auto v_pair: e_pair.first->get_boundary())
    //   //     v_cells.emplace(v_pair.first);
    //   // }
    //   // for (auto vc: v_cells) {
    //   //   int ci_value = ci_map.at(vc);
    //   //   barycenter[0] += vertex_points[ci_value-1][0] / v_cells.size();
    //   //   barycenter[1] += vertex_points[ci_value-1][1] / v_cells.size();
    //   //   barycenter[2] += vertex_points[ci_value-1][2] / v_cells.size();
    //   // }
    //   // vertex_points.push_back(barycenter);
    //   // for (auto e_pair: h->get_boundary()) {
    //   //   std::vector<int> triangle(1, vertex_points.size());
    //   //   for (auto v_pair: e_pair.first->get_boundary())
    //   //     triangle.push_back(ci_map.at(v_pair.first));
    //   //   triangles.push_back(triangle);
    //   //   filtrations.push_back(h->get_filtration());
    //   //   mask.push_back(mask_val);
    //   // }
    // }
    // for (auto p: polytopes) {
    //   // // constrain the outer edges, which are found by the convex hull
    //   // std::vector<Point_3> vertices;
    //   // std::vector<int> v_indices;
    //   // std::set<Hasse_cell_ptr> v_cells;
    //   // // std::cout << "3-cell " << p->get_filtration() << ":\n";
    //   // for (auto h_pair: p->get_boundary()) {
    //   //   // std::cout << "2-cell " << h_pair.first->get_position() << " " << h_pair.first->get_filtration() << ":\n";
    //   //   for (auto e_pair: h_pair.first->get_boundary()) {
    //   //     // std::cout << "Edge " << e_pair.first->get_position() << " " << e_pair.first->get_filtration() << ":\n"; 
    //   //     for (auto v_pair: e_pair.first->get_boundary()) {
    //   //       // std::cout << "Vertex " << v_pair.first->get_position() << " " << v_pair.first->get_filtration() << ":\n"; 
    //   //       v_cells.emplace(v_pair.first);
    //   //     }
    //   //   }
    //   //   // std::cout << "\n";
    //   // }
    //   // for (auto vc: v_cells) {
    //   //   Point_t& b = vertex_points[ci_map.at(vc)-1];
    //   //   vertices.push_back(Point_3(b[0], b[1], b[2]));
    //   //   v_indices.push_back(ci_map.at(vc));
    //   //   // std::cout << "Vertex " << vc->get_position() << " " << vc->get_filtration() << ":\n"; 
    //   // }
    //   // Delaunay_triangulation del(3);
    //   // index = 0;
    //   // std::map<typename Delaunay_triangulation::Vertex_handle, int> index_of_vertex;
    //   // for (auto pt: vertices)
    //   //   index_of_vertex.emplace(del.insert(pt), index++);
    //   // for (auto fc_it = del.full_cells_begin(); fc_it != del.full_cells_end(); ++fc_it) {
    //   //   if (del.is_infinite(fc_it))
    //   //     continue;
    //   //   std::vector<int> tetrahedron;
    //   //   for (auto fv_it = fc_it->vertices_begin(); fv_it != fc_it->vertices_end(); ++fv_it)
    //   //     tetrahedron.push_back(v_indices[index_of_vertex[*fv_it]]);
    //   //   tetrahedra.push_back(tetrahedron);
    //   //   filtrations.push_back(p->get_filtration());
    //   // }
    // }
    // ofs << "MeshVersionFormatted 1\nDimension 3\n";
    // ofs_bb << "3 1 ";
    // ofs << "Vertices\n" << vertex_points.size() << "\n";
    // for (auto p: vertex_points) {
    //   ofs << p[0] << " " << p[1] << " " << p[2] << " 215\n";
    // }
    // ofs << "Edges " << edges.size() << "\n";
    // for (auto s: edges) {
    //   for (auto v: s->get_boundary()) {
    // 	ofs << ci_map.at(v.first) << " ";
    //   }
    //   ofs << "515" << std::endl;
    // }
    // ofs << "Triangles " << triangles.size() << "\n";
    // ofs_bb << triangles.size()+tetrahedra.size() << " 1\n";
    // auto m_it = mask.begin();
    // auto f_it = filtrations.begin();
    // for (auto s: triangles) {
    //   for (auto v: s) {
    // 	ofs << v << " ";
    //   }
    //   ofs << *m_it++ << std::endl;
    //   ofs_bb << *f_it++ << "\n";
    // }
    // ofs << "Tetrahedra " << tetrahedra.size() << "\n";
    // for (auto s: tetrahedra) {
    //   for (auto v: s) {
    // 	ofs << v << " ";
    //   }
    //   ofs << "545\n";
    //   ofs_bb << *f_it++ << "\n";
    // }
      
  }
  ofs.close();
  ofs_bb.close();
}


#endif
