project(Alpha_complex_benchmark)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable(Alpha_complex_3d_benchmark Alpha_complex_3d_benchmark.cpp)
  target_link_libraries(Alpha_complex_3d_benchmark ${CGAL_LIBRARY})
  if(TARGET TBB::tbb)
    target_link_libraries(Alpha_complex_3d_benchmark TBB::tbb)
  endif()
endif ()
