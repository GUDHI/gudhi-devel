add_library(_default_random SHARED _default_random.cc)
# Required by Windows wheel
install(TARGETS _default_random DESTINATION ${CMAKE_INSTALL_PREFIX}/gudhi)
# _default_random needs to define EXPORT at compilation, and when you link with it you need to IMPORT
target_compile_definitions(_default_random
                           PRIVATE GUDHI_DEFAULT_RANDOM_DLL_EXPORT
                           INTERFACE GUDHI_DEFAULT_RANDOM_DLL_IMPORT)

# This function assumes that THE_TARGET is directly under gudhi/ and not in a subdirectory
function( link_with_default_random THE_TARGET )
  target_link_libraries(${THE_TARGET} PRIVATE _default_random)
  add_dependencies(${THE_TARGET} _default_random)
  if(APPLE)
    set_target_properties(${THE_TARGET} PROPERTIES BUILD_RPATH "@loader_path" INSTALL_RPATH "@loader_path")
  elseif(UNIX)
    set_target_properties(${THE_TARGET} PROPERTIES BUILD_RPATH "$ORIGIN" INSTALL_RPATH "$ORIGIN")
  endif()
endfunction( link_with_default_random )

if (TARGET CGAL::CGAL)
  add_library(_default_cgal_random SHARED _default_cgal_random.cc)
  # Required by Windows wheel
  install(TARGETS _default_cgal_random DESTINATION ${CMAKE_INSTALL_PREFIX}/gudhi)
  # _default_cgal_random needs to define EXPORT at compilation, and when you link with it you need to IMPORT
  target_compile_definitions(_default_cgal_random
                             PRIVATE GUDHI_DEFAULT_CGAL_RANDOM_DLL_EXPORT
                             INTERFACE GUDHI_DEFAULT_CGAL_RANDOM_DLL_IMPORT)
endif()

nanobind_add_module_with_targets(_pers_cub_low_dim_ext _pers_cub_low_dim.cc TARGETS TBB::tbb)
nanobind_add_module_with_targets(_ripser_ext _ripser.cc)
nanobind_add_module_with_targets(_reader_utils_ext _reader_utils.cc)
nanobind_add_module_with_targets(_hera_ext _hera.cc)
nanobind_add_module_with_targets(_bottleneck_ext _bottleneck.cc TARGETS CGAL::CGAL TBB::tbb)

nanobind_add_module_with_targets(_nerve_gic_ext _nerve_gic.cc TARGETS TBB::tbb)
link_with_default_random(_nerve_gic_ext)

nanobind_add_module_with_targets(_subsampling_ext _subsampling.cc TARGETS Eigen3::Eigen CGAL::CGAL TBB::tbb)
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  link_with_default_random(_subsampling_ext)
endif()

nanobind_add_module_with_targets(_default_random_seed_ext _default_random_seed.cc)
link_with_default_random(_default_random_seed_ext)

if (TARGET CGAL::CGAL)
  nanobind_add_module_with_targets(_default_cgal_random_seed_ext _default_cgal_random_seed.cc TARGETS CGAL::CGAL)
  target_link_libraries(_default_cgal_random_seed_ext PRIVATE _default_cgal_random)
  add_dependencies(_default_cgal_random_seed_ext _default_cgal_random)
  if(APPLE)
    set_target_properties(_default_cgal_random_seed_ext PROPERTIES BUILD_RPATH "@loader_path" INSTALL_RPATH "@loader_path")
  elseif(UNIX)
    set_target_properties(_default_cgal_random_seed_ext PROPERTIES BUILD_RPATH "$ORIGIN" INSTALL_RPATH "$ORIGIN")
  endif()
endif()
