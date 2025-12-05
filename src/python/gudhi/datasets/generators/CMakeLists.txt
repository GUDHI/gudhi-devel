nanobind_add_module_with_targets(_points_ext _points.cc TARGETS Python::NumPy Eigen3::Eigen CGAL::CGAL)
# For Windows, where points is a client of Random_generator here (_random_ext is the provider)
target_link_libraries(_points_ext INTERFACE _random_ext)
