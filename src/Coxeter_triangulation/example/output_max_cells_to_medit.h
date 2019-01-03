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
    ofs << "MeshVersionFormatted 1\nDimension 2\n";
    ofs_bb << "2 1 ";
    ofs << "Vertices\n" << vertex_points.size() << "\n";
    for (auto p: vertex_points) {
      ofs << p[0] << " " << p[1] << " 215\n";
    }
    ofs << "Edges " << max_cells.size() << "\n";
    for (auto e: max_cells) {
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(e, ct, d);
      Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
      for (; cof_it != cof_end; ++cof_it) {
    	ofs << ci_map.at(*cof_it) << " ";
      }
      ofs << "515" << std::endl;
    }
  } // d == 3
  else {
    std::unordered_set<Cell_id> edges;
    std::vector<std::vector<int> > triangles;
    std::vector<double> filtrations;
    std::vector<std::size_t> mask;
    if (max_cells.begin()->dimension() == 2) {
      ofs << "MeshVersionFormatted 1\nDimension 3\n";
      ofs_bb << "3 1 ";
      ofs << "Vertices\n" << vertex_points.size() << "\n";
      for (auto p: vertex_points) {
	ofs << p[0] << " " << p[1] << " 215\n";
      }
      ofs << "Edges " << max_cells.size() << "\n";
      for (auto e: max_cells) {
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(e, ct, d);
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
	for (; cof_it != cof_end; ++cof_it) {
	  ofs << ci_map.at(*cof_it) << " ";
	}
	ofs << "515" << std::endl;
      }
    }
    else { // max_cells.begin()->dimension() == 1
      for (auto c: max_cells) {
	std::size_t mask_val = 515;
        Point_d barycenter(3);
        barycenter[0] = 0;
        barycenter[1] = 0;
        barycenter[2] = 0;
        std::unordered_set<Cell_id> v_cells;
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(c, ct, 2);
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
	for (; cof_it != cof_end; ++cof_it) {
	  std::size_t num_cof_cells = 0;
	  for (auto f: ct.face_range(*cof_it, 1)) 
	    if (max_cells.find(f) != max_cells.end())
	      num_cof_cells++;
	  if (num_cof_cells > 2)
	    mask_val = 517;
	  edges.emplace(*cof_it);
	}
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof2_it(c, ct, 3);
	for (; cof2_it != cof_end; ++cof2_it) {
          int ci_value = ci_map.at(*cof2_it);
          barycenter[0] += vertex_points[ci_value-1][0] / v_cells.size();
          barycenter[1] += vertex_points[ci_value-1][1] / v_cells.size();
          barycenter[2] += vertex_points[ci_value-1][2] / v_cells.size();
        }
        vertex_points.push_back(barycenter);
        Gudhi::Coxeter_triangulation_ds::Coface_iterator cof3_it(c, ct, 2);
	for (; cof3_it != cof_end; ++cof3_it) {
          std::vector<int> triangle(1, vertex_points.size());
	  Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(*cof3_it, ct, 3);
          for (; cof_it != cof_end; ++cof_it)
            triangle.push_back(ci_map.at(*cof_it));
          triangles.push_back(triangle);
          filtrations.push_back(0);
          mask.push_back(mask_val);
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
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_it(e, ct, d);
	Gudhi::Coxeter_triangulation_ds::Coface_iterator cof_end;
	for (; cof_it != cof_end; ++cof_it) {
	  ofs << ci_map.at(*cof_it) << " ";
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
