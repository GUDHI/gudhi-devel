#include <gudhi/Freudenthal_triangulation.h>

int main() {
  Gudhi::coxeter_triangulation::Freudenthal_triangulation<> tr(3);
  std::vector<double> point({3,-1,0});
  std::cout << tr.locate_point(point) << "\n";
  std::vector<double> point2({3.5,-1.5,0.5});
  std::cout << tr.locate_point(point2) << "\n";
  std::vector<double> point3({3.5,-1.8,0.5});
  std::cout << tr.locate_point(point3) << "\n";
  std::vector<double> point4({3.5,-1.8,0.3});
  std::cout << tr.locate_point(point4) << "\n";

  std::cout << "Triangulation's dimension is " << tr.dimension() << "\n";
  std::cout << "Triangulation's matrix is\n" << tr.matrix() << "\n";
  std::cout << "Triangulation's offset is\n" << tr.offset() << "\n";
  Eigen::MatrixXd new_matrix(3,3);
  new_matrix << 1, 0, 0, 0, 1, 0, -1, 0, 1;
  std::cout << "Determinant of the new matrix is " << new_matrix.determinant() << "\n";
  Eigen::Vector3d new_offset(1.5, 1, 0.5);
  
  tr.change_matrix(new_matrix);
  tr.change_offset(new_offset);
  std::cout << "Triangulation's dimension is " << tr.dimension() << "\n";
  std::cout << "Triangulation's matrix is\n" << tr.matrix() << "\n";
  std::cout << "Triangulation's offset is\n" << tr.offset() << "\n";
  return 0;
}
