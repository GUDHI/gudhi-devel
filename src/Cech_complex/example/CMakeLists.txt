project(Cech_complex_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.0.1)
  add_executable ( Cech_complex_example_from_points cech_complex_example_from_points.cpp)
  if (TBB_FOUND)
    target_link_libraries(Cech_complex_example_from_points ${TBB_LIBRARIES})
  endif()
  add_test(NAME Cech_complex_example_from_points COMMAND $<TARGET_FILE:Cech_complex_example_from_points>)
endif()
