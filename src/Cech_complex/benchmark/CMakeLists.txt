project(Cech_complex_benchmark)

# Cech complex requires CGAL >= 5.0.1
if(TARGET CGAL::CGAL AND EIGEN3_FOUND AND TARGET Boost::filesystem AND NOT CGAL_VERSION VERSION_LESS 5.0.1)
  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_executable(cech_complex_benchmark cech_complex_benchmark.cpp)
  target_link_libraries(cech_complex_benchmark CGAL::CGAL Boost::filesystem)
  
  if(TARGET TBB::tbb)
    target_link_libraries(cech_complex_benchmark TBB::tbb)
  endif()
endif()
