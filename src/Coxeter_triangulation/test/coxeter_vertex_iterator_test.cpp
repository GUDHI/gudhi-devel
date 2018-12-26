#include <gudhi/Coxeter_triangulation_ds.h>
#include <gudhi_patches/CGAL/Epick_d.h>
#include <gudhi/random_point_generators.h>

bool is_face(Gudhi::Cell_id& c1, Gudhi::Cell_id& c2) {
  std::size_t k = 0; 
  for (; k < c1.size(); ++k) {
    if (c2.mask(k)) {
      if (!c1.mask(k) || c1.value(k) != c2.value(k))
	return false;
    }
    else { // !c2.mask(k)
      if (!c1.mask(k) && c1.value(k) != c2.value(k))
	return false;
      if (c1.mask(k))
	if (c1.value(k) < c2.value(k) || c1.value(k) > c2.value(k)+1)
	  return false;
    } 
  }
  return true;
}
  
int main() {
  Gudhi::Coxeter_triangulation_ds ct_ds_2(2);
  Gudhi::Cell_id c(1.5, 3);
  c.push_back(1, false);
  c.push_back(-1, false);
  c.push_back(0, false);
  c.push_back(2, false);
  c.push_back(2, false);
  c.push_back(3, false);
  
  std::cout << "Cell: " << c << "\nVertices:\n";
  for (Gudhi::Cell_id v: ct_ds_2.vertex_range(c))
    std::cout << v << ". Is a face = " << is_face(v,c) << ".\n";

  // Generate 500 random points in a box, point locate them, compute vertices and check is_face
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_d = Kernel::Point_d;

  std::size_t d = 20, n = 500;
  Gudhi::Coxeter_triangulation_ds ct_ds_d(d);
  std::vector<Point_d > points = Gudhi::generate_points_in_cube_d<Kernel>(n, d, 500.0);
  for (std::vector<double> p: points) {
    Gudhi::Cell_id c = ct_ds_d.locate_point(p, 1);
    // cell == point located barycenter
    if (c != ct_ds_d.locate_point(ct_ds_d.barycenter(c))) {
	std::cerr << "Error: cell = " << c << ", ct_ds_d.locate_point(ct_ds_d.barycenter(c)) = " << ct_ds_d.locate_point(ct_ds_d.barycenter(c)) << ".\n";
	return 1;
    }
    std::size_t number_of_vertices = 0;
    for (Gudhi::Cell_id v: ct_ds_d.vertex_range(c)) {
      // v is a face of c
      if (!is_face(v,c)) {
	std::cerr << "Error: cell = " << c << ", v = " << v << ".\n";
	return 1;
      }
      // v == point located cartesian coordinates of v
      if (!(v == ct_ds_d.locate_point(ct_ds_d.cartesian_coordinates(v)))) {
	std::cerr << "Error: vertex = " << v << ", ct_ds_d.locate_point(ct_ds_d.cartesian_coordinates(v)) = " << ct_ds_d.locate_point(ct_ds_d.cartesian_coordinates(v)) << ".\n";
      }
      number_of_vertices++;
    }
    // number of vertices is c.dimension()+1
    if (number_of_vertices != c.dimension()+1) {
      std::cerr << "Error: cell = " << c << ", n_vertices = " << number_of_vertices << ".\n";
      return 1;
    }
  }
  std::cout << "Success!\n";
  return 0;
}

