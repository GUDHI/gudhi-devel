project(Tangential_complex_benchmark)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable(Tangential_complex_benchmark benchmark_tc.cpp)
  target_link_libraries(Tangential_complex_benchmark CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Tangential_complex_benchmark TBB::tbb)
  endif()
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
