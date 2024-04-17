project(Tangential_complex_benchmark)

if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable(Tangential_complex_benchmark benchmark_tc.cpp)
  target_link_libraries(Tangential_complex_benchmark CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Tangential_complex_benchmark TBB::tbb)
  endif()
endif(CGAL_FOUND AND EIGEN3_FOUND)
