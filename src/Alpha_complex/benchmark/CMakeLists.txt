project(Alpha_complex_benchmark)

if(CGAL_FOUND AND TARGET Eigen3::Eigen)
  add_executable(Alpha_complex_3d_benchmark Alpha_complex_3d_benchmark.cpp)
  target_link_libraries(Alpha_complex_3d_benchmark CGAL::CGAL Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(Alpha_complex_3d_benchmark TBB::tbb)
  endif()
endif ()
