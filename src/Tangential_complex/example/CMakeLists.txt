project(Tangential_complex_examples)

if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable( Tangential_complex_example_basic example_basic.cpp )
  target_link_libraries(Tangential_complex_example_basic CGAL::CGAL)
  add_executable( Tangential_complex_example_with_perturb example_with_perturb.cpp )
  target_link_libraries(Tangential_complex_example_with_perturb CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Tangential_complex_example_basic TBB::tbb)
    target_link_libraries(Tangential_complex_example_with_perturb TBB::tbb)
  endif()

  add_test(NAME Tangential_complex_example_basic
      COMMAND $<TARGET_FILE:Tangential_complex_example_basic>)
  add_test(NAME Tangential_complex_example_with_perturb
      COMMAND $<TARGET_FILE:Tangential_complex_example_with_perturb>)

endif(CGAL_FOUND AND EIGEN3_FOUND)
