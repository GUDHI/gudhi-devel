if (TARGET CGAL::CGAL)
  nanobind_add_module_with_targets(_points_ext _points.cc TARGETS Eigen3::Eigen CGAL::CGAL)
  target_link_libraries(_points_ext PRIVATE _default_cgal_random)
  add_dependencies(_points_ext _default_cgal_random)
  # CMake translates "$ORIGIN/../.." to "@loader_path/../.." on MacOS. Potentially unreliable on MacOS < 11
  # Does not work on Windows - cf. datasets/generators/points.py on how PATH is modified before importing the libs
  set_target_properties(_points_ext PROPERTIES BUILD_RPATH "$ORIGIN/../.." INSTALL_RPATH "$ORIGIN/../..")
endif()
