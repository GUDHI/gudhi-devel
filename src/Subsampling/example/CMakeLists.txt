cmake_minimum_required(VERSION 2.6)
project(Subsampling_examples)

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_executable( Subsampling_example_pick_n_random_points example_pick_n_random_points.cpp )

      add_executable( Subsampling_example_choose_n_farthest_points example_choose_n_farthest_points.cpp )
      #target_link_libraries(Subsampling_example_choose_n_farthest_points ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      
      add_executable(Subsampling_example_sparsify_point_set example_sparsify_point_set.cpp)
      target_link_libraries(Subsampling_example_sparsify_point_set ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Subsampling feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Subsampling examples. Version 4.8.0 is required.")
  endif ()
else()
  message(WARNING "CGAL not found. It is required for the Subsampling examples.")
endif()
