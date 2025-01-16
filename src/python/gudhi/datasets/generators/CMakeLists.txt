nanobind_add_module(_points _points.cc)
target_link_libraries(_points PRIVATE CGAL::CGAL Eigen3::Eigen)
