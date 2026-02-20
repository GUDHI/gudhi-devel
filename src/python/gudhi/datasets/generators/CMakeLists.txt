nanobind_add_module_with_targets(_points_ext _points.cc TARGETS Python::NumPy Eigen3::Eigen CGAL::CGAL random)
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_dependencies(_points_ext random)
endif()
