nanobind_add_module_with_targets(_points_ext _points.cc TARGETS Python::NumPy Eigen3::Eigen CGAL::CGAL)
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  target_compile_definitions(_points_ext PUBLIC RANDOM_DLL_IMPORT)
endif()
