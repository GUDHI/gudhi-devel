nanobind_add_module_with_targets(_tangential_complex_ext _tangential_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)

nanobind_add_module_with_targets(_witness_complex_ext _witness_complex.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_strong_witness_complex_ext _strong_witness_complex.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_euclidean_witness_complex_ext _euclidean_witness_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
nanobind_add_module_with_targets(_euclidean_strong_witness_complex_ext _euclidean_strong_witness_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
