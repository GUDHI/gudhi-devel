#include "../example/cxx-prettyprint/prettyprint.hpp"
#include <gudhi/Coxeter_triangulation_ds.h>
#include <gudhi/Coxeter_triangulation/restricted_partition_range.h>
#include <gudhi_patches/CGAL/Epick_d.h>
#include <gudhi/random_point_generators.h>

int main() {
  // std::size_t n = 8, k = 2;
  // for (std::vector<size_t> part: base_n_range(n,k)) {
  //   std::cout << part << "\n";
  // }
  // std::size_t n = 8, k = 2;
  // std::vector<std::size_t> restrictions {6, 3};
  // for (std::vector<size_t> part: restricted_partition_range(n, k, restrictions)) {
  //   std::cout << part << "\n";
  // }

  Gudhi::Coxeter_triangulation_ds ct_ds_2(3);
  Gudhi::Cell_id c(1.5, 1);
  // c.push_back(1, true);
  // c.push_back(-1, false);  
  // c.push_back(0, false);
  // c.push_back(2, true);
  // c.push_back(1, false);
  // c.push_back(2, false);
  c.push_back(0, true);
  c.push_back(0, true);  
  c.push_back(0, true);
  c.push_back(0, false);
  c.push_back(0, false);
  c.push_back(0, false);
  std::cout << "c = " << c << "\n";
  // Gudhi::Cell_id v = *Gudhi::Coxeter_triangulation_ds::Vertex_iterator(c, ct_ds_2);
  // for (auto cf: ct_ds_2.coface_range(c, 3)) {}
  auto cof_it = Gudhi::Coxeter_triangulation_ds::Coface_iterator(c, ct_ds_2, 3);
  auto cof_end = Gudhi::Coxeter_triangulation_ds::Coface_iterator();
  for (; cof_it != cof_end; ++cof_it)
    std::cout << *cof_it << "\n";
  
  // Generate 500 random points in a box, point locate them, compute vertices and check is_face
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_d = Kernel::Point_d;

  std::size_t d = 8, n = 10;
  Gudhi::Coxeter_triangulation_ds ct_ds_d(d);
  std::vector<Point_d > points = Gudhi::generate_points_in_cube_d<Kernel>(n, d, 500.0);
  for (std::vector<double> p: points) {
    Gudhi::Cell_id c = ct_ds_d.locate_point(p, 1);
    // cell == point located barycenter
    for (std::size_t h1 = 0; h1 <= c.dimension(); ++h1)
      for (auto f: ct_ds_2.face_range(c, h1))
	for (std::size_t h = h1; h <= d; ++h) {
	  auto cof_it = Gudhi::Coxeter_triangulation_ds::Coface_iterator(f, ct_ds_2, h);
	  auto cof_end = Gudhi::Coxeter_triangulation_ds::Coface_iterator();
	  for (; cof_it != cof_end; ++cof_it) {
	    // f is a face of c
	    if (!Gudhi::Cell_id::is_face(f, *cof_it)) {
	      std::cerr << "Error: coface = " << *cof_it << ", f = " << f << ".\n";
	      return 1;
	    }
	  }
	}
  }
  std::cout << "Success!\n";
  return 0;
}

