nanobind_add_module(_pers_cub_low_dim _pers_cub_low_dim.cc)

nanobind_add_module(_cubical_complex_ext _cubical_complex.cc)

nanobind_add_module(bottleneck bottleneck.cc)
nanobind_add_module(hera hera.cc)
nanobind_add_module(simplex_tree_ext simplex_tree.cc)

nanobind_add_module(delaunay_complex_ext delaunay_complex.cc)
target_link_libraries(delaunay_complex_ext PRIVATE Eigen3::Eigen)
target_link_libraries(delaunay_complex_ext PRIVATE CGAL::CGAL)

if(TARGET TBB::tbb)
    target_link_libraries(_pers_cub_low_dim PRIVATE TBB::tbb)
    target_link_libraries(_cubical_complex_ext PRIVATE TBB::tbb)
    target_link_libraries(bottleneck PRIVATE TBB::tbb)
endif()
