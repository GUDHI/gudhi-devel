project(Cech_complex_examples)

# Cech complex requires CGAL >= 5.0.1
if(CGAL_FOUND AND NOT CGAL_VERSION VERSION_LESS 5.0.1 AND TARGET Eigen3::Eigen)
  add_executable ( Cech_complex_example_from_points cech_complex_example_from_points.cpp)
  target_link_libraries(Cech_complex_example_from_points CGAL::CGAL Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(Cech_complex_example_from_points TBB::tbb)
  endif()
  add_test(NAME Cech_complex_example_from_points COMMAND $<TARGET_FILE:Cech_complex_example_from_points>)
endif()
