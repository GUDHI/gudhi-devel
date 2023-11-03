project(Tangential_complex_benchmark)

if(CGAL_FOUND AND TARGET Eigen3::Eigen)
  add_executable(Tangential_complex_benchmark benchmark_tc.cpp)
  target_link_libraries(Tangential_complex_benchmark CGAL::CGAL Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(Tangential_complex_benchmark TBB::tbb)
  endif()
endif()
