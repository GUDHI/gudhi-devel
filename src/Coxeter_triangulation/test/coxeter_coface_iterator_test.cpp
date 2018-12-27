#include "../example/cxx-prettyprint/prettyprint.hpp"
#include <gudhi/Coxeter_triangulation_ds.h>
// #include <gudhi/Coxeter_triangulation/base_n_range.h>
// #include <gudhi_patches/CGAL/Epick_d.h>
// #include <gudhi/random_point_generators.h>

int main() {
  // std::size_t n = 8, k = 2;
  // for (std::vector<size_t> part: base_n_range(n,k)) {
  //   std::cout << part << "\n";
  // }

  Gudhi::Coxeter_triangulation_ds ct_ds_2(3);
  Gudhi::Cell_id c(1.5, 3);
  c.push_back(1, true);
  c.push_back(-1, false);  
  c.push_back(0, false);
  c.push_back(2, false);
  c.push_back(2, false);
  c.push_back(3, false);
  // Gudhi::Cell_id v = *Gudhi::Coxeter_triangulation_ds::Vertex_iterator(c, ct_ds_2);
  Gudhi::Coxeter_triangulation_ds::Coface_iterator(c, ct_ds_2, 2);
  
  
  // std::cout << "Cell: " << c << "\nVertices:\n";
  // for (Gudhi::Cell_id f: ct_ds_2.face_range(c, 1))
  //   std::cout << f << ". Is a face = " << Gudhi::Cell_id::is_face(f,c) << ".\n";

  // // Generate 500 random points in a box, point locate them, compute vertices and check is_face
  // using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  // using Point_d = Kernel::Point_d;

  // std::size_t d = 8, n = 500;
  // Gudhi::Coxeter_triangulation_ds ct_ds_d(d);
  // std::vector<Point_d > points = Gudhi::generate_points_in_cube_d<Kernel>(n, d, 500.0);
  // for (std::vector<double> p: points) {
  //   Gudhi::Cell_id c = ct_ds_d.locate_point(p, 1);
  //   // cell == point located barycenter
  //   if (c != ct_ds_d.locate_point(ct_ds_d.barycenter(c))) {
  // 	std::cerr << "Error: cell = " << c << ", ct_ds_d.locate_point(ct_ds_d.barycenter(c)) = " << ct_ds_d.locate_point(ct_ds_d.barycenter(c)) << ".\n";
  // 	return 1;
  //   }
  //   std::size_t number_of_faces = 0;
  //   for (std::size_t h = 0; h <= c.dimension(); ++h)
  //     for (Gudhi::Cell_id f: ct_ds_d.face_range(c, h)) {
  // 	// f is a face of c
  // 	if (!Gudhi::Cell_id::is_face(f, c)) {
  // 	  std::cerr << "Error: cell = " << c << ", f = " << f << ".\n";
  // 	  return 1;
  // 	}
  // 	// f == point located barycenter of f
  // 	if (f != ct_ds_d.locate_point(ct_ds_d.barycenter(f))) {
  // 	  std::cerr << "Error: face = " << f << ", ct_ds_d.locate_point(ct_ds_d.barycenter(f)) = " << ct_ds_d.locate_point(ct_ds_d.barycenter(f)) << ".\n";
  // 	}
  // 	number_of_faces++;
  //     }
  //   // number of vertices is c.dimension()+1
  //   if (number_of_faces != std::exp2(c.dimension()+1)-1) {
  //     std::cerr << "Error: cell = " << c << ", num_faces = " << number_of_faces << ".\n";
  //     for (std::size_t h = 0; h <= c.dimension(); ++h)
  // 	for (Gudhi::Cell_id f: ct_ds_d.face_range(c, h))
  // 	  std::cout << f << "\n";
  //     return 1;
  //   }
  // }
  // std::cout << "Success!\n";
  return 0;
}

