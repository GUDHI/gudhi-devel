nanobind_add_module_with_targets(_points_ext _points.cc TARGETS Python::NumPy Eigen3::Eigen CGAL::CGAL _random_ext)
# For Windows, where points is a client of Random_generator here (_random_ext is the provider)
add_dependencies(_points_ext _random_ext)
