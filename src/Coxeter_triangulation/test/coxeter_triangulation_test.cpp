#include <gudhi/Coxeter_triangulation_ds.h>


int main() {
  Gudhi::Coxeter_triangulation_ds ct_ds(3);
  std::cout << "ct_ds.dimension() = " << ct_ds.dimension() << "\n";
  std::cout << "ct_ds.pos_root_count() = " << ct_ds.pos_root_count() << "\n";
  std::cout << "ct_ds.simple_root_matrix() =\n" << ct_ds.simple_root_matrix() << "\n";

  Gudhi::Cell_id c_id_0(0);
  c_id_0.push_back(0, true);
  c_id_0.push_back(1, true);
  c_id_0.push_back(1, true);
  c_id_0.push_back(0, true);
  c_id_0.push_back(1, true);
  c_id_0.push_back(1, true);
  std::cout << "ct_ds.alcove_dimension(" << c_id_0 << ") = " << ct_ds.alcove_dimension(c_id_0) << "\n";

  Gudhi::Cell_id c_id_1(0);
  c_id_1.push_back(0, true);
  c_id_1.push_back(1, false);
  c_id_1.push_back(1, false);
  c_id_1.push_back(0, true);
  c_id_1.push_back(1, false);
  c_id_1.push_back(1, false);
  std::cout << "ct_ds.alcove_dimension(" << c_id_1 << ") = " << ct_ds.alcove_dimension(c_id_1) << "\n";

  Gudhi::Cell_id c_id_2(0);
  c_id_2.push_back(0, true);
  c_id_2.push_back(1, false);
  c_id_2.push_back(1, false);
  c_id_2.push_back(0, false);
  c_id_2.push_back(1, false);
  c_id_2.push_back(1, false);
  std::cout << "ct_ds.alcove_dimension(" << c_id_2 << ") = " << ct_ds.alcove_dimension(c_id_2) << "\n";

  Gudhi::Cell_id c_id_3(0);
  c_id_3.push_back(0, false);
  c_id_3.push_back(1, false);
  c_id_3.push_back(1, false);
  c_id_3.push_back(0, false);
  c_id_3.push_back(1, false);
  c_id_3.push_back(1, false);
  std::cout << "ct_ds.alcove_dimension(" << c_id_3 << ") = " << ct_ds.alcove_dimension(c_id_3) << "\n";

  std::vector<double> p1 {1.5, 0.2, -1.7};
  std::cout << "ct_ds.locate_point(p1) = " << ct_ds.locate_point(p1) << "\n";

  std::vector<double> p0 {0, 0, 0};
  std::cout << "ct_ds.locate_point(0) = " << ct_ds.locate_point(p0) << "\n";

  return 0;
}

