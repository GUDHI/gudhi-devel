project(Cech_complex_benchmark)

if (CGAL_FOUND AND EIGEN3_FOUND)
  if(TARGET Boost::filesystem)
    # Do not forget to copy test files in current binary dir
    file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

    add_executable(cech_complex_benchmark cech_complex_benchmark.cpp)
    target_link_libraries(cech_complex_benchmark CGAL::CGAL Eigen3::Eigen Boost::filesystem)
  
    if(TARGET TBB::tbb)
      target_link_libraries(cech_complex_benchmark TBB::tbb)
    endif()
  endif()
endif()
