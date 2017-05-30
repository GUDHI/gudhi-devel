cmake_minimum_required(VERSION 2.6)
project(Tangential_complex_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable( Tangential_complex_example_basic example_basic.cpp )
  target_link_libraries(Tangential_complex_example_basic ${CGAL_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
  add_executable( Tangential_complex_example_with_perturb example_with_perturb.cpp )
  target_link_libraries(Tangential_complex_example_with_perturb ${CGAL_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Tangential_complex_example_basic ${TBB_LIBRARIES})
    target_link_libraries(Tangential_complex_example_with_perturb ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(NAME Tangential_complex_example_basic
      COMMAND $<TARGET_FILE:Tangential_complex_example_basic>)
  add_test(NAME Tangential_complex_example_with_perturb
      COMMAND $<TARGET_FILE:Tangential_complex_example_with_perturb>)

  install(TARGETS Tangential_complex_example_basic DESTINATION bin)
  install(TARGETS Tangential_complex_example_with_perturb DESTINATION bin)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
