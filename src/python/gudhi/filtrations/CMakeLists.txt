# Sparse Rips uses choose_n_farthest that uses Random
nanobind_add_module_with_targets(_rips_complex_ext _rips_complex.cc TARGETS TBB::tbb)
link_with_default_random(_rips_complex_ext)

# delaunay complex is split in several files to avoid to much memory consumption during compilation
nanobind_add_module_with_targets(_delaunay_complex_ext _delaunay_complex.cc _delaunay_complex_2d.cc _delaunay_complex_3d.cc _delaunay_complex_dyn.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
