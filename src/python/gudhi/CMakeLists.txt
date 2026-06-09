add_library(_default_random SHARED _default_random.cc)
# Required by Windows wheel
install(TARGETS _default_random DESTINATION ${CMAKE_INSTALL_PREFIX}/gudhi)
# _default_random needs to be the on that EXPORT, and when you link with it you need to IMPORT
target_compile_definitions(_default_random
                           PRIVATE GUDHI_DEFAULT_RANDOM_DLL_EXPORT
                           INTERFACE GUDHI_DEFAULT_RANDOM_DLL_IMPORT)

function( link_with_default_random THE_TARGET )
  target_link_libraries(${THE_TARGET} PRIVATE _default_random)
  add_dependencies(${THE_TARGET} _default_random)
  if(NOT WIN32)
    # CMake translates "$ORIGIN" to "@loader_path" on MacOS. Windows searchs in local path, then in PATH
    set_target_properties(${THE_TARGET} PROPERTIES BUILD_RPATH "$ORIGIN" INSTALL_RPATH "$ORIGIN")
  endif()
endfunction( link_with_default_random )

if (TARGET CGAL::CGAL)
  add_library(_default_cgal_random SHARED _default_cgal_random.cc)
  # Required by Windows wheel
  install(TARGETS _default_cgal_random DESTINATION ${CMAKE_INSTALL_PREFIX}/gudhi)
  # _default_cgal_random needs to be the on that EXPORT, and when you link with it you need to IMPORT
  target_compile_definitions(_default_cgal_random
                             PRIVATE GUDHI_DEFAULT_CGAL_RANDOM_DLL_EXPORT
                             INTERFACE GUDHI_DEFAULT_CGAL_RANDOM_DLL_IMPORT)
endif()

nanobind_add_module_with_targets(_pers_cub_low_dim_ext _pers_cub_low_dim.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_cubical_complex_ext _cubical_complex.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_ripser_ext _ripser.cc)
nanobind_add_module_with_targets(_edge_collapse_ext _edge_collapse.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_simplex_tree_ext _simplex_tree.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_reader_utils_ext _reader_utils.cc)
nanobind_add_module_with_targets(_witness_complex_ext _witness_complex.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_strong_witness_complex_ext _strong_witness_complex.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_hera_ext _hera.cc)
nanobind_add_module_with_targets(_bottleneck_ext _bottleneck.cc TARGETS CGAL::CGAL TBB::tbb)
# delaunay complex is split in several files to avoid to much memory consumption during compilation
nanobind_add_module_with_targets(_delaunay_complex_ext _delaunay_complex.cc _delaunay_complex_2d.cc _delaunay_complex_3d.cc _delaunay_complex_dyn.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
nanobind_add_module_with_targets(_tangential_complex_ext _tangential_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
nanobind_add_module_with_targets(_euclidean_witness_complex_ext _euclidean_witness_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
nanobind_add_module_with_targets(_euclidean_strong_witness_complex_ext _euclidean_strong_witness_complex.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)

nanobind_add_module_with_targets(_nerve_gic_ext _nerve_gic.cc TARGETS TBB::tbb)
link_with_default_random(_nerve_gic_ext)

# Sparse Rips uses choose_n_farthest that uses Random
nanobind_add_module_with_targets(_rips_complex_ext _rips_complex.cc TARGETS TBB::tbb)
link_with_default_random(_rips_complex_ext)

nanobind_add_module_with_targets(_subsampling_ext _subsampling.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
link_with_default_random(_subsampling_ext)

nanobind_add_module_with_targets(_default_random_seed_ext _default_random_seed.cc)
link_with_default_random(_default_random_seed_ext)

if (TARGET CGAL::CGAL)
  nanobind_add_module_with_targets(_default_cgal_random_seed_ext _default_cgal_random_seed.cc TARGETS CGAL::CGAL)
  target_link_libraries(_default_cgal_random_seed_ext PRIVATE _default_cgal_random)
  add_dependencies(_default_cgal_random_seed_ext _default_cgal_random)
  # CMake translates "$ORIGIN" to "@loader_path" on MacOS. Windows searchs in local path, then in PATH
  set_target_properties(_default_cgal_random_seed_ext PROPERTIES BUILD_RPATH "$ORIGIN" INSTALL_RPATH "$ORIGIN")
endif()
