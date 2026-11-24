nanobind_add_module_with_targets(_hera_ext _hera.cc)
nanobind_add_module_with_targets(_bottleneck_ext _bottleneck.cc TARGETS CGAL::CGAL TBB::tbb)
