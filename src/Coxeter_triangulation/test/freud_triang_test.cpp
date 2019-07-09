#include <gudhi/Freudenthal_triangulation.h>
#include <gudhi/Coxeter_triangulation.h>

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
  new_matrix << 1, 0, 0, -1, 1, 0, -1, 0, 1;
  std::cout << "Determinant of the new matrix is " << new_matrix.determinant() << "\n";
  Eigen::Vector3d new_offset(1.5, 1, 0.5);
  
  tr.change_matrix(new_matrix);
  tr.change_offset(new_offset);
  std::cout << "Triangulation's new matrix is\n" << tr.matrix() << "\n";
  std::cout << "Triangulation's new offset is\n" << tr.offset() << "\n";

  std::cout << "The simplex containing point is " << tr.locate_point(point) << "\n";
  std::cout << "The barycenter of the simplex containing point is:\n"
	    << tr.barycenter(tr.locate_point(point)) << "\n";
  std::cout << "The simplex containing point (with scaling) is " << tr.locate_point(point, 2) << "\n";
  std::cout << "The barycenter of the simplex containing point (with scaling) is:\n"
	    << tr.barycenter(tr.locate_point(point, 2), 2) << "\n";

  Gudhi::coxeter_triangulation::Coxeter_triangulation<> cox_tr(3);
  std::cout << "Coxeter triangulation's matrix is\n" << cox_tr.matrix() << "\n";
  std::cout << "Coxeter triangulation's offset is\n" << cox_tr.offset() << "\n";

  return 0;
}
