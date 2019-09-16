#include <gudhi/Freudenthal_triangulation.h>
#include <gudhi/Coxeter_triangulation.h>

int main(int argc, char** argv) {
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

  {
    std::size_t d = atoi(argv[1]);
    double level = atof(argv[2]);
    Gudhi::coxeter_triangulation::Coxeter_triangulation<> cox_tr(d);
    cox_tr.change_matrix(level * cox_tr.matrix());
    cox_tr.change_offset(Eigen::VectorXd::Random(d));
    double diam = 0;
    auto s = cox_tr.locate_point(Eigen::VectorXd::Zero(d));
    for (auto v1: s.vertex_range())
      for (auto v2: s.vertex_range()) {
	Eigen::VectorXd p1 = cox_tr.cartesian_coordinates(v1);
	Eigen::VectorXd p2 = cox_tr.cartesian_coordinates(v2);
	double dist = (p1-p2).norm();
	if (dist > diam)
	  diam = dist;
      }
    std::cout << "The diameter of a simplex in a " << d << "-dimensional Coxeter triangulation with level " << level << " = " << diam << "\n";
  }
  
  return 0;
}
