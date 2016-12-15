cmake_minimum_required(VERSION 2.6)
project(Tangential_complex_examples)

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)
      add_executable( Tangential_complex_example_basic example_basic.cpp )
      target_link_libraries(Tangential_complex_example_basic ${CGAL_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
      add_executable( Tangential_complex_example_with_perturb example_with_perturb.cpp )
      target_link_libraries(Tangential_complex_example_with_perturb ${CGAL_LIBRARY} ${Boost_DATE_TIME_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(Tangential_complex_example_basic ${TBB_LIBRARIES})
        target_link_libraries(Tangential_complex_example_with_perturb ${TBB_LIBRARIES})
      endif(TBB_FOUND)

      add_test(Tangential_complex_example_basic
               ${CMAKE_CURRENT_BINARY_DIR}/Tangential_complex_example_basic)

      add_test(Tangential_complex_example_with_perturb
               ${CMAKE_CURRENT_BINARY_DIR}/Tangential_complex_example_with_perturb)
    endif(EIGEN3_FOUND)
  endif(NOT CGAL_VERSION VERSION_LESS 4.8.0)
endif(CGAL_FOUND)
