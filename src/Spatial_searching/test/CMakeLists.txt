project(Spatial_searching_tests)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_test_coverage)

  add_executable( Spatial_searching_test_Kd_tree_search test_Kd_tree_search.cpp )
  target_link_libraries(Spatial_searching_test_Kd_tree_search
    ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  gudhi_add_coverage_test(Spatial_searching_test_Kd_tree_search)
endif ()
